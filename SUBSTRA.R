library(pROC)
library(clusteval)
library(foreach)
library(doParallel)


# NOTICE: The current version assumes BINARY 'data' and 'labels'. The next versions will be more general and applicable for categorical and continuous values.


# This function splits the dataset (columns) into 'nfolds' parts with equal ratio of positive to negative 'labels'.
# The output is a list with two elements:
#   The first element is a list of parts of 'data'
#   The second element is a list of parts of 'labels'
#   (the first element of the first element corresponds to the first element of the second element, the second element of the first element corresponds to the second element of the second element, etc.)

SUBSTRA.split_dataset = function(data, labels, nfolds){
  minv=min(labels)
  maxv=max(labels)
  cases=which(labels==maxv)
  controls=which(labels==minv)
  casestops=round(c(0,c(1:nfolds)*length(cases)/nfolds))
  controlstops=round(c(0,c(1:nfolds)*length(controls)/nfolds))
  dataparts=list()
  labelparts=list()
  for(i in 1:nfolds){
    dataparts[[i]]=cbind.data.frame(data[,cases[(casestops[i]+1):casestops[i+1]]],
                                    data[,controls[(controlstops[i]+1):controlstops[i+1]]])
    labelparts[[i]]=c(labels[cases[(casestops[i]+1):casestops[i+1]]],
                      labels[controls[(controlstops[i]+1):controlstops[i+1]]])
  }
  
  output=list()
  output[[1]]=dataparts
  output[[2]]=labelparts
  output
}



# This function performs supervised biclustering assuming the the columns of 'data' has labels' (the target features): length(labels)==ncol(data)
# The function assings weights of given 'magnitude' to the rows of 'data', which indicates their relevance to the 'labels'.
# The Output is a list containing the following elements:
#   column.clusters: the cluster indexes for columns of 'data'
#   row.clusters: the cluster indexes for rows of 'data'
#   bicluster.counts: the numbers of 1's and 2's (or 0's and 1's) inside each bicluster (the first matrix correponds to 1's and the second to 2's)
#   counts.for.rows: the numbers of 1's and 2's inside each column cluster associated with each row (two matrices corresponding to the numbers of 1's and 2's)
#   counts.for.columns: the numbers of 1's and 2's inside each row cluster associated with each column (two matrices corresponding to the numbers of 1's and 2's)
#   counts.for.labels: a matrix with each column corresponding to one of the values of 'labels' and each row corresponding to one column clusters
#   Pi: estimated value of the model parameter for column clusters
#   Theta: estimated value of the model parameter for biclusters (rows correspond to row clusters and columns correspond to column clusters)
#   Psi: estimated value of the model parameter for labels inside each column cluster
#   row.weights: values of weights assigned to each row/feature

SUBSTRA.train <- function(data, labels,
                          magnitude=1,       # magnitude of feature/row weights
                          phase_1_ite=100,   # iterations with equal feature weights
                          phase_2_ite=50,    # iterations for learning the weights
                          verbose=T,
                          alpha_l = 1e-100, # hyper-parameter for labels (small value if mislabeling is low and vice versa)
                          alpha_p = 1,       # hyper-parameter for column clustering
                          alpha_t = 1,       # hyper-parameter for row clustering
                          alpha_e = 1){      # hyper-parameter for bicluster parameters
  
  data=as.matrix(data)
  data = data + (1-min(data))
  m=max(data)
  data[data==m]=2
  labels = labels + (1-min(labels))
  Nv = length(unique(labels))
  Np = ncol(data)
  Kp=Np  # maximum number of column clusters
  Nt = nrow(data)
  Kt=Nt  # maximum number of row clusters
  logDBLMax=log(.Machine$double.xmax)
  # Clustering Variables
  kp = c()       # column cluster indexes
  kt = c()       # row cluster indexes
  npINkp = c()   # sizes of column clusters
  ntINkt = c()   # sizes of row clusters
  # Bicluster Counters
  nesINkpkt = array(0,dim = c(Kp, Kt, 2))   # counts inside each bicluster
  neINpFORkt = array(0,dim = c(Np, Kt, 2))  # counts for each column inside each row cluster
  neINtFORkp = array(0,dim = c(Nt, Kp, 2))  # counts for each column inside each row cluster
  row.weights=rep(magnitude,Nt)
  # Labels Counters
  nvsINKp =  array(0,dim = c(Kp, Nv))       # counts of labels inside each column cluster
  

  ########### Phase 0: Random Initialization ############
  
  for(i in 1:Np)
    kp[i] = sample(seq(labels[i],round(sqrt(Np)),by=2),1)
  kt = sample(c(1:round(sqrt(Nt))),Nt,replace = T)
  for(i in 1:Kp){
    kpinds=which(kp==i)
    npINkp[i]=length(kpinds)
    if(npINkp[i]!=0){
      for(v in 1:Nv) nvsINKp[i,v]=length(which(labels[kpinds]==v))
      
      for(j in 1:Kt){
        ktinds=which(kt==j)
        if(length(ktinds)==0) next
        nesINkpkt[i,j,1]=sum(2-data[ktinds,kpinds])
        nesINkpkt[i,j,2]=sum(data[ktinds,kpinds]-1)
      }
      if(npINkp[i]==1){
        neINtFORkp[,i,1]= (2-data[,kpinds])
        neINtFORkp[,i,2]= (data[,kpinds]-1)
      }
      else{
        neINtFORkp[,i,1]= rowSums(2-data[,kpinds])
        neINtFORkp[,i,2]= rowSums(data[,kpinds]-1)
      }
    }
  }
  for(j in 1:Kt){
    ktinds=which(kt==j)
    ntINkt[j]=length(ktinds)
    if(ntINkt[j]!=0){
      if(ntINkt[j]==1){
        neINpFORkt[,j,1]=(2-data[ktinds,])
        neINpFORkt[,j,2]=(data[ktinds,]-1)  
      }else{
        neINpFORkt[,j,1]=colSums(2-data[ktinds,])
        neINpFORkt[,j,2]=colSums(data[ktinds,]-1)
      }
    }
  }
  
  #Initial Status
  occupied.kp=unique(kp)
  occupied.kt=unique(kt)
  empty.kp=c()
  empty.kt=c()
  if(length(occupied.kp) < Kp)
    empty.kp=min(which(npINkp==0))
  if(length(occupied.kt) < Kt)
    empty.kt=min(which(ntINkt==0))
  
  Pi=npINkp[occupied.kp]/Np
  Theta=(nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e/2)/
    (nesINkpkt[occupied.kp,occupied.kt,1]+nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e)
  Psi=array(nvsINKp[occupied.kp,]+alpha_l/Nv,dim=c(length(occupied.kp),2))/array(rep(npINkp[occupied.kp]+alpha_l,2),dim=c(length(occupied.kp),2))
  best.model=list(column.clusters=match(kp,occupied.kp), row.clusters=match(kt,occupied.kt),
                  bicluster.counts=nesINkpkt[occupied.kp,occupied.kt,], counts.for.rows=neINtFORkp[,occupied.kp,], counts.for.pats=neINpFORkt[,occupied.kt,], counts.for.labels=nvsINKp[occupied.kp,],
                  Pi=Pi,Theta=t(Theta),Psi=Psi,row.weights=row.weights)
  preds=SUBSTRA.predict(best.model,data)
  r=roc(as.factor(labels),preds[,2],direction = "<")
  bestauc=r$auc
  besterr=mean((labels-1-preds[,2])^2)
  if(verbose){
    print("Initial Status:")
    print(paste0("Train Set AUC: ", bestauc))
    print(paste0("Train Set MSE: ", besterr))
  }
  
  # Sampling
  
  ################### PHASE I ###################
  
  for(ite in 1:phase_1_ite){
    model.prob = 0
    prev.kp=kp
    prev.kt=kt
    
    ## Row Clusters
    t.diff = 0
    for(t in 1:Nt){
      ### take out
      tempkt=kt[t]
      ntINkt[tempkt]=ntINkt[tempkt]-1
      if(ntINkt[tempkt]==0){
        occupied.kt=setdiff(occupied.kt,tempkt)
        empty.kt=tempkt
      }
      nesINkpkt[occupied.kp,tempkt,] = nesINkpkt[occupied.kp,tempkt,] - neINtFORkp[t,occupied.kp,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] - (2-data[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] - (data[t,]-1)
      
      sig.kp=c(occupied.kp,empty.kp)
      sig.kt=c(occupied.kt,empty.kt)
      ### compute probability
      Theta=(nesINkpkt[sig.kp,sig.kt,2]+alpha_e/2)/(nesINkpkt[sig.kp,sig.kt,1]+nesINkpkt[sig.kp,sig.kt,2]+alpha_e)
      t.prob=colSums(log(1-Theta)*matrix(rep(neINtFORkp[t,sig.kp,1],length(sig.kt)),nrow = length(sig.kp),ncol = length(sig.kt),byrow = F))+
        colSums(log(Theta)*matrix(rep(neINtFORkp[t,sig.kp,2],length(sig.kt)),nrow = length(sig.kp),ncol = length(sig.kt),byrow = F))
      
      t.prob=t.prob+log(c(ntINkt[occupied.kt],alpha_t))

      
      ### put back
      t.prob.max=max(t.prob)
      t.prob.max = t.prob.max - (logDBLMax - log(length(sig.kt)))
      t.prob = exp(t.prob - t.prob.max)
      total.prob = sum(t.prob)
      tempkt=sample(sig.kt,1,prob = t.prob/total.prob)
      if (kt[t] != tempkt)
        t.diff=t.diff+1
      model.prob = model.prob + log(t.prob[tempkt] / total.prob)
      
      kt[t]=tempkt
      ntINkt[tempkt]=ntINkt[tempkt]+1
      nesINkpkt[occupied.kp,tempkt,] = nesINkpkt[occupied.kp,tempkt,] + neINtFORkp[t,occupied.kp,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] + (2-data[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] + (data[t,]-1)
      
      if(tempkt==empty.kt){
        occupied.kt=unique(kt)
        empty.kt=c()
        if(length(occupied.kt) < Kt)
          empty.kt=min(which(ntINkt==0))
      }
    }#t
    

    ## Column Clusters
    p.diff = 0
    for(p in 1:Np){
      ### take out
      tempkp=kp[p]
      npINkp[tempkp]=npINkp[tempkp]-1
      if(npINkp[tempkp]==0){
        occupied.kp=setdiff(occupied.kp,tempkp)
        empty.kp=tempkp
      }
      nesINkpkt[tempkp,occupied.kt,] = nesINkpkt[tempkp,occupied.kt,] - neINpFORkt[p,occupied.kt,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] - (2-data[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] - (data[,p]-1)
      
      nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]-1
      
      ### compute probability
      sig.kp=c(occupied.kp,empty.kp)
      sig.kt=c(occupied.kt,empty.kt)
      reverseData=matrix(rep(2-data[,p],length(sig.kp)),nrow = length(sig.kp),ncol = Nt,byrow = T)
      Theta=(nesINkpkt[sig.kp,sig.kt,2]+alpha_e/2)/(nesINkpkt[sig.kp,sig.kt,1]+nesINkpkt[sig.kp,sig.kt,2]+alpha_e)
      p.prob=rowSums(log(abs(reverseData-Theta[,match(kt,sig.kt)]))*
                       matrix(rep(row.weights,length(sig.kp)),nrow = length(sig.kp),ncol = Nt,byrow = T))+
        log((nvsINKp[sig.kp,labels[p]]+alpha_l/Nv)/(npINkp[sig.kp]+alpha_l))
      
      p.prob=p.prob+log(c(npINkp[occupied.kp],alpha_p))
      
      ### put back
      p.prob.max=max(p.prob)
      p.prob.max = p.prob.max - (logDBLMax - log(length(sig.kp)))
      p.prob = exp(p.prob - p.prob.max)
      total.prob = sum(p.prob)
      tempkp=sample(sig.kp,1,prob = p.prob/total.prob)
      if (kp[p] != tempkp)
        p.diff=p.diff+1
      model.prob = model.prob + log(p.prob[tempkp] / total.prob)
      
      kp[p]=tempkp
      npINkp[tempkp]=npINkp[tempkp]+1
      nesINkpkt[tempkp,occupied.kt,] = nesINkpkt[tempkp,occupied.kt,] + neINpFORkt[p,occupied.kt,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] + (2-data[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] + (data[,p]-1)
      
      nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]+1
      
      if(tempkp==empty.kp){
        occupied.kp=unique(kp)
        empty.kp=c()
        if(length(occupied.kp) < Kp)
          empty.kp=min(which(npINkp==0))
      }
    }#p
    
    Pi=npINkp[occupied.kp]/Np
    Theta=(nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e/2)/
      (nesINkpkt[occupied.kp,occupied.kt,1]+nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e)
    Psi=(nvsINKp[occupied.kp,]+alpha_l/Nv)/array(rep(npINkp[occupied.kp]+alpha_l,2),dim=c(length(occupied.kp),2))
    model=list(column.clusters=match(kp,occupied.kp), row.clusters=match(kt,occupied.kt),
               bicluster.counts=nesINkpkt[occupied.kp,occupied.kt,], counts.for.rows=neINtFORkp[,occupied.kp,], counts.for.pats=neINpFORkt[,occupied.kt,], counts.for.labels=nvsINKp[occupied.kp,],
               Pi=Pi,Theta=t(Theta), Psi=Psi,row.weights=row.weights)
    
    preds=SUBSTRA.predict(model,data)
    r=roc(as.factor(labels),preds[,2],direction = "<")
    err=mean((labels-1-preds[,2])^2)
    besterr=err
    bestauc=r$auc
    best.model=model
    
    rkp=cluster_similarity(kp,prev.kp)
    rkt=cluster_similarity(kt,prev.kt)
    if(verbose){
      print("")
      print("***************************************************************")
      print(paste0("Phase I - Iteration: ",ite))
      print("Changes from the previous iteration:")
      print(paste0("Column Clusters: ",p.diff," - ",length(occupied.kp)," clusters occupied."))
      print(paste0("Row Clusters: ",t.diff," - ",length(occupied.kt)," clusters occupied."))
      print("")
      print("Similarity with the previous iteration (Rand Indexes):")
      print(paste0("Column Clustering:    ",rkp))
      print(paste0("Row Clustering: ",rkt))
      print("")
      print(paste0("Current Weights (Min - Average - Max): ",min(row.weights)," - ",mean(row.weights)," - ",max(row.weights)))
      print("")
      print(paste0("Train Set AUC: ", r$auc))
      print(paste0("Train Set MSE: ",err))
    }
    
    if(rkp > 0.95 && rkt > 0.95){
      if(verbose){
        print("")
        print("Converged!!! Going to Phase II!!!")
      }
      break
    }
  }#phase_1_ite
  
  
  
  ################### PHASE II ###################
  
  for(ite in 1:phase_2_ite){
    model.prob = 0
    prev.kp=kp
    prev.kt=kt
    
    ## Row Clusters
    t.diff = 0
    for(t in 1:Nt){
      ### take out
      tempkt=kt[t]
      ntINkt[tempkt]=ntINkt[tempkt]-1
      if(ntINkt[tempkt]==0){
        occupied.kt=setdiff(occupied.kt,tempkt)
        empty.kt=tempkt
      }
      nesINkpkt[occupied.kp,tempkt,] = nesINkpkt[occupied.kp,tempkt,] - neINtFORkp[t,occupied.kp,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] - (2-data[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] - (data[t,]-1)
      
      if(length(which(nesINkpkt<0))>0)
        print(t)
      
      sig.kp=c(occupied.kp,empty.kp)
      sig.kt=c(occupied.kt,empty.kt)
      ### compute probability
      Theta=(nesINkpkt[sig.kp,sig.kt,2]+alpha_e/2)/(nesINkpkt[sig.kp,sig.kt,1]+nesINkpkt[sig.kp,sig.kt,2]+alpha_e)
      t.prob=colSums(log(1-Theta)*matrix(rep(neINtFORkp[t,sig.kp,1],length(sig.kt)),nrow = length(sig.kp),ncol = length(sig.kt),byrow = F))+
        colSums(log(Theta)*matrix(rep(neINtFORkp[t,sig.kp,2],length(sig.kt)),nrow = length(sig.kp),ncol = length(sig.kt),byrow = F))
   
      t.prob=t.prob+log(c(ntINkt[occupied.kt],alpha_t))
      
      ### put back
      t.prob.max=max(t.prob)
      t.prob.max = t.prob.max - (logDBLMax - log(length(sig.kt)))
      t.prob = exp(t.prob - t.prob.max)
      total.prob = sum(t.prob)
      tempkt=sample(sig.kt,1,prob = t.prob/total.prob)
      if (kt[t] != tempkt)
        t.diff=t.diff+1
      model.prob = model.prob + log(t.prob[tempkt] / total.prob)
      
      kt[t]=tempkt
      ntINkt[tempkt]=ntINkt[tempkt]+1
      nesINkpkt[occupied.kp,tempkt,] = nesINkpkt[occupied.kp,tempkt,] + neINtFORkp[t,occupied.kp,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] + (2-data[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] + (data[t,]-1)
      
      if(tempkt==empty.kt){
        occupied.kt=unique(kt)
        empty.kt=c()
        if(length(occupied.kt) < Kt)
          empty.kt=min(which(ntINkt==0))
      }
    }#t
    
    ## Column Clusters
    p.diff = 0
    prev.weights=row.weights
    for(p in 1:Np){
      ### take out
      tempkp=kp[p]
      npINkp[tempkp]=npINkp[tempkp]-1
      if(npINkp[tempkp]==0){
        occupied.kp=setdiff(occupied.kp,tempkp)
        empty.kp=tempkp
      }
      nesINkpkt[tempkp,occupied.kt,] = nesINkpkt[tempkp,occupied.kt,] - neINpFORkt[p,occupied.kt,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] - (2-data[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] - (data[,p]-1)
      
      nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]-1
      
      ### compute weights
      Theta=(nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e/2)/
        (nesINkpkt[occupied.kp,occupied.kt,1]+nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e)
      Pi=npINkp[occupied.kp]/Np
      Psi=(nvsINKp[occupied.kp,]+alpha_l/Nv)/array(rep(npINkp[occupied.kp]+alpha_l,2),dim=c(length(occupied.kp),2))
      
      # Error before adjusting the weights
      BO.model=list(column.clusters=match(kp,occupied.kp), row.clusters=match(kt,occupied.kt),
                    bicluster.counts=nesINkpkt[occupied.kp,occupied.kt,], counts.for.rows=neINtFORkp[,occupied.kp,], counts.for.pats=neINpFORkt[,occupied.kt,], counts.for.labels=nvsINKp[occupied.kp,],
                    Pi=Pi,Theta=t(Theta), Psi=Psi,row.weights=row.weights)
      preds=SUBSTRA.predict(BO.model,as.matrix(data[,p]))
      BOerr=(labels[p]-1-preds[2])^2
      
      # Derivative
      reverseData=matrix(rep(2-data[,p],length(occupied.kp)),nrow = length(occupied.kp),ncol = Nt,byrow = T)
      logTheta=log(abs(reverseData-Theta[,match(kt,occupied.kt)]))
      temp.weights=row.weights
      temp.weights.matrix=matrix(rep(temp.weights,length(occupied.kp)),nrow = length(occupied.kp),ncol = Nt,byrow = T)
      logs=rowSums(logTheta*temp.weights.matrix)+log(Pi)
      logs=logs+log(.Machine$double.xmax)-max(logs)-log(length(logs))
      values=exp(logs)
      values=values/sum(values)
      A=sum(values*Psi[,2])
      B=1
      dA=colSums(matrix(rep(values*Psi[,2],Nt),byrow = F,nrow = length(occupied.kp),ncol = Nt)*logTheta)
      dB=colSums(matrix(rep(values,Nt),byrow = F,nrow = length(occupied.kp),ncol = Nt)*logTheta)
      f=A/B
      df=(dA*B-dB*A)/(B*B)
      # Updating the weights
      temp.weights=temp.weights+magnitude*2*(labels[p]-1-f)*df
      temp.weights[which(temp.weights<0)]=0
      
      # Checking if the weights improve the prediction accuracy
      AO.model=list(column.clusters=match(kp,occupied.kp), row.clusters=match(kt,occupied.kt),
                    bicluster.counts=nesINkpkt[occupied.kp,occupied.kt,], counts.for.rows=neINtFORkp[,occupied.kp,], counts.for.pats=neINpFORkt[,occupied.kt,], counts.for.labels=nvsINKp[occupied.kp,],
                    Pi=Pi,Theta=t(Theta), Psi=Psi,row.weights=temp.weights)
      preds=SUBSTRA.predict(AO.model,as.matrix(data[,p]))
      AOerr=(labels[p]-1-preds[2])^2
      
      if(AOerr < BOerr){
        row.weights=temp.weights
      }
      
      ### compute probability
      sig.kp=c(occupied.kp,empty.kp)
      sig.kt=c(occupied.kt,empty.kt)
      reverseData=matrix(rep(2-data[,p],length(sig.kp)),nrow = length(sig.kp),ncol = Nt,byrow = T)
      Theta=(nesINkpkt[sig.kp,sig.kt,2]+alpha_e/2)/(nesINkpkt[sig.kp,sig.kt,1]+nesINkpkt[sig.kp,sig.kt,2]+alpha_e)
      p.prob=rowSums(log(abs(reverseData-Theta[,match(kt,sig.kt)]))*
                       matrix(rep(row.weights,length(sig.kp)),nrow = length(sig.kp),ncol = Nt,byrow = T))+
        log((nvsINKp[sig.kp,labels[p]]+alpha_l/Nv)/(npINkp[sig.kp]+alpha_l))
      
      p.prob=p.prob+log(c(npINkp[occupied.kp],alpha_p))
      
      ### put back
      p.prob.max=max(p.prob)
      p.prob.max = p.prob.max - (logDBLMax - log(length(sig.kp)))
      p.prob = exp(p.prob - p.prob.max)
      total.prob = sum(p.prob)
      tempkp=sample(sig.kp,1,prob = p.prob/total.prob)
      if (kp[p] != tempkp)
        p.diff=p.diff+1
      model.prob = model.prob + log(p.prob[tempkp] / total.prob)
      
      kp[p]=tempkp
      npINkp[tempkp]=npINkp[tempkp]+1
      nesINkpkt[tempkp,occupied.kt,] = nesINkpkt[tempkp,occupied.kt,] + neINpFORkt[p,occupied.kt,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] + (2-data[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] + (data[,p]-1)
      
      nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]+1
      
      if(tempkp==empty.kp){
        occupied.kp=unique(kp)
        empty.kp=c()
        if(length(occupied.kp) < Kp)
          empty.kp=min(which(npINkp==0))
      }
    }#p
    
    w.diff=length(which(row.weights!=prev.weights))
    
    Pi=npINkp[occupied.kp]/Np
    Theta=(nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e/2)/
      (nesINkpkt[occupied.kp,occupied.kt,1]+nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e)
    Psi=(nvsINKp[occupied.kp,]+alpha_l/Nv)/array(rep(npINkp[occupied.kp]+alpha_l,2),dim=c(length(occupied.kp),2))
    model=list(column.clusters=match(kp,occupied.kp), row.clusters=match(kt,occupied.kt),
               bicluster.counts=nesINkpkt[occupied.kp,occupied.kt,], counts.for.rows=neINtFORkp[,occupied.kp,], counts.for.columns=neINpFORkt[,occupied.kt,], counts.for.labels=nvsINKp[occupied.kp,],
               Pi=Pi,Theta=t(Theta), Psi=Psi,row.weights=row.weights)
    
    preds=SUBSTRA.predict(model,data)
    r=roc(as.factor(labels),preds[,2],direction = "<")
    err=mean((labels-1-preds[,2])^2)
    
    rkp=cluster_similarity(kp,prev.kp)
    rkt=cluster_similarity(kt,prev.kt)
    if(verbose){
      print("")
      print("***************************************************************")
      print(paste0("Phase II - Iteration: ",ite))
      print("Changes from the previous iteration:")
      print(paste0("Column Clusters: ",p.diff," - ",length(occupied.kp)," clusters occupied."))
      print(paste0("Row Clusters: ",t.diff," - ",length(occupied.kt)," clusters occupied."))
      print(paste0("Row Weights: ",w.diff))
      print("")
      print("Similarity with the previous iteration (Rand Indexes):")
      print(paste0("Column Clustering:    ",rkp))
      print(paste0("Row Clustering: ",rkt))
      print("")
      print(paste0("Current Weights (Min - Average - Max): ",min(row.weights)," - ",mean(row.weights)," - ",max(row.weights)))
      print("")
      print(paste0("Train Set AUC: ", r$auc))
      print(paste0("Train Set MSE: ",err))
    }
    
    if(r$auc > bestauc){
      best.model=model
      besterr=err
      bestauc=r$auc
      if(verbose){
        print("")
        print("Best Model Updated!!!")
      }
    }else if(r$auc == bestauc){
      if(err <= besterr){
        best.model=model
        besterr=err
        bestauc=r$auc
        if(verbose){
          print("")
          print("Best Model Updated!!!")
        }
      }
    }
  }#phase_2_ite
  
  return(best.model)
}



# This function predicts the labels for 'data' columns using the SUBSTRA 'model'.
# The output is a matrix with each column corresponding to a possible output value and each row corresponding to a column of 'data'.
# The elements of that matrix are probabilities of the corresponding output value for the corresponding column of 'data'.

SUBSTRA.predict <- function(model, data){
  data=as.matrix(data)
  data = data + (1-min(data))
  
  probs=matrix(0,nrow = ncol(data),ncol = ncol(model$Psi))
  for(i in 1:ncol(data)){
    reverseData=matrix(rep(2-data[,i],length(model$Pi)),nrow = length(model$Pi),ncol = nrow(data),byrow = T)
    logs=rowSums(log(abs(reverseData-t(model$Theta[model$row.clusters,])))*
                   matrix(rep(model$row.weights,length(model$Pi)),nrow=length(model$Pi),ncol=nrow(data),byrow=T))+
      log(model$Pi)
    logs=logs+log(.Machine$double.xmax)-max(logs)-log(length(logs))
    values=exp(logs)
    clusterProbs=values/sum(values)
    probs[i,] = colSums(matrix(rep(clusterProbs,2),nrow=nrow(model$Psi),ncol=ncol(model$Psi))*model$Psi)
  }
  return(round(probs,6))
}



# This function tunes the feature weight 'magnitude' for a dataset through 'nfold'-fold cross-validation
# The output is the 'magnitude' value with the best cross-validation AUC (Area Under the ROC Curve).
# This value should be given as input to the 'SUBSTRA.train' function when training with 'data'

SUBSTRA.tune <- function(data, labels,
                         pathToSUBSTRA,     # the address of SUBSTRA.R file
                         magnitudes=c(0.001,0.01,0.1,1,10),       # the set of magnitudes of feature weights to be checked
                         nfolds=3,          # number of cross-validation folds for parameter tuning
                         phase_1_ite=80,    # iterations with equal feature weights
                         phase_2_ite=30,    # iterations for learning the weights
                         parallel=T,        # whether the cross-validation iterations should be executed in parallel
                         ncores=3,          # number of CPU cores available to this function (used if parellel=TRUE)
                         verbose=T,
                         trainVerbose=F,    # wether the printings of 'SUBSTRA.train' function should be displayed
                         alpha_l = 1e-100, # hyper-parameter for labels (small value if mislabeling is low and vice versa)
                         alpha_p = 1,       # hyper-parameter for column clustering
                         alpha_t = 1,       # hyper-parameter for row clustering
                         alpha_e = 1){      # hyper-parameter for bicluster parameters
  
  data=as.matrix(data)
  dataparts=SUBSTRA.split_dataset(data,labels,nfolds)
  bestAUC=0
  bestMagnitude=magnitudes[1]
  for(magnitude in magnitudes){
    if(verbose)
      print(paste0("Checking magnitude=",magnitude))
    if(parallel){
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      results=foreach(i=1:nfolds, .combine = rbind, .packages = c("pROC","clusteval","e1071")) %dopar% {
        source(pathToSUBSTRA)
        trainData=cbind.data.frame(dataparts[[1]][-i])
        trainLabs=unlist(dataparts[[2]][-i])
        testData=as.data.frame(dataparts[[1]][i])
        testLabs=unlist(dataparts[[2]][i])
        
        model=SUBSTRA.train(trainData,trainLabs,phase_1_ite = phase_1_ite,phase_2_ite = phase_2_ite,verbose = trainVerbose,
                            magnitude = magnitude,alpha_p = alpha_p,alpha_t = alpha_t,alpha_l = alpha_l,alpha_e = alpha_e)
        ps=SUBSTRA.predict(model,testData)[,2]
        r=roc(as.factor(testLabs), ps, direction = "<")
        length(testLabs)*r$auc
      }
      stopCluster(cl)
      registerDoSEQ()
      AUC=sum(results)/length(labels)
    }else{
      AUC=0
      for(i in 1:nfolds){
        trainData=cbind.data.frame(dataparts[[1]][-i])
        trainLabs=unlist(dataparts[[2]][-i])
        testData=as.data.frame(dataparts[[1]][i])
        testLabs=unlist(dataparts[[2]][i])
        
        model=SUBSTRA.train(trainData,trainLabs,phase_1_ite = phase_1_ite,phase_2_ite = phase_2_ite,verbose = trainVerbose,
                            magnitude = magnitude,alpha_p = alpha_p,alpha_t = alpha_t,alpha_l = alpha_l,alpha_e = alpha_e)
        ps=SUBSTRA.predict(model,testData)[,2]
        r=roc(as.factor(testLabs), ps, direction = "<")
        AUC=AUC+length(testLabs)*r$auc
      }
      AUC=AUC/length(labels)
    }
    
    if(verbose)
      print(paste0("Achieved AUC = ",AUC))
    
    if(AUC > bestAUC){
      bestAUC=AUC
      bestMagnitude=magnitude
    }
  }
  
  if(verbose)
    print(paste0("Selected 'magnitude': ",bestMagnitude))
  
  bestMagnitude
}



# This function trains 'ntimes' models and aggregates them into one model through consensus clustering
# The Output is a list containing the following elements:
#   column.clusters: the cluster indexes for columns of 'data'
#   row.clusters: the cluster indexes for rows of 'data'
#   row.weights: values of weights assigned to each row/feature

SUBSTRA.train_consensus <- function(data, labels,
                                    ntimes=5,          # the number of models to be trained and aggregated
                                    parallel=T,        # whether the models should be trained in parallel
                                    ncores=5,          # number of CPU cores available to this function (used if parellel=TRUE)
                                    verbose=T,
                                    ## Tuning Arguments
                                    pathToSUBSTRA,     # the address of SUBSTRA.R file
                                    magnitudes=c(0.001,0.01,0.1,1,10),     # the set of magnitudes of feature weights to be checked
                                    nfolds_tune=3,     # number of cross-validation folds for parameter tuning
                                    phase_1_ite_tune=80,      # Phase I iterations for SUBSTRA.tune
                                    phase_2_ite_tune=30,      # Phase II iterations for SUBSTRA.tune
                                    parallel_tune=T,   # whether the cross-validation iterations of 'SUBSTRA.tune' should be executed in parallel
                                    verbose_tune=F,    # wether the printings of 'SUBSTRA.tune' function should be displayed
                                    trainVerbose_tune=F,      # wether the printings of 'SUBSTRA.train' inside the SUBSTRA.tune' function should be displayed
                                    ## Training Arguments
                                    phase_1_ite_train=100,    # Phase I iterations for SUBSTRA.train
                                    phase_2_ite_train=50,     # Phase II iterations for SUBSTRA.train
                                    verbose_train=F,  # wether the printings of 'SUBSTRA.train' function should be displayed
                                    ## Common Arguments
                                    alpha_l = 1e-100, # hyper-parameter for labels (small value if mislabeling is low and vice versa)
                                    alpha_p = 1,       # hyper-parameter for column clustering
                                    alpha_t = 1,       # hyper-parameter for row clustering
                                    alpha_e = 1){      # hyper-parameter for bicluster parameters
  
  # Tuning
  if(verbose)
    print("Tuning the 'magnitude'...")
  
  magnitude=SUBSTRA.tune(data,labels,pathToSUBSTRA,magnitudes,nfolds_tune,phase_1_ite_tune,phase_1_ite_tune,parallel_tune,
                         min(ncores,nfolds_tune),verbose_tune,trainVerbose_tune,alpha_l,alpha_p,alpha_t,alpha_e)
  if(verbose)
    print(paste0("Selected 'magnitude': ",magnitude))
  
  # Training
  if(verbose)
    print("Training the individual models...")
  
  results.kp=matrix(0,nrow = ntimes,ncol = ncol(data))
  results.kt=results.w=matrix(0,nrow = ntimes,ncol = nrow(data))
  if(parallel){
    cl=makeCluster(ncores)
    registerDoParallel(cl)
    results=foreach(i=1:ntimes, .combine = rbind, .packages = c("pROC","clusteval")) %dopar% {
      source(pathToSUBSTRA)
      model=SUBSTRA.train(data,labels,magnitude,phase_1_ite_train,phase_2_ite_train,verbose_train,alpha_l,alpha_p,alpha_t,alpha_e)
      rbind(model$column.clusters,model$row.clusters,rep(model$row.weights,ceiling(ncol(data)/nrow(data))))
    }
    stopCluster(cl)
    registerDoSEQ()
    results.kp=results[seq(1,ntimes*3,3),1:ncol(data)]
    results.kt=results[seq(2,ntimes*3,3),1:nrow(data)]
    results.w=results[seq(3,ntimes*3,3),1:nrow(data)]
  }else{
    for(i in 1:ntimes){
      model=SUBSTRA.train(data,labels,magnitude,phase_1_ite_train,phase_2_ite_train,verbose_train,alpha_l,alpha_p,alpha_t,alpha_e)
      results.kp[i,]=model$column.clusters
      results.kt[i,]=model$row.clusters
      results.w[i,]=model$row.weights
    }
  }
  
  # Consensus
  if(verbose)
    print("Aggregating the models...")
  
  SUBSTRA.aggregate(results.kp,results.kt,results.w)
}
  


# This function aggregates the results from several models through consensus clustering (based on hierarchical clustering)
# The output is a list containing aggregate column clusters, aggregate, row clusters, and aggregate row weights

SUBSTRA.aggregate <- function(results.kp,  # column.clusters from different models
                              results.kt,  # row.clusters from different models
                              results.w){  # row.weights from different models
  
  ntimes=nrow(results.kp)
  ## Column Clusters
  results.kp=apply(results.kp,2,as.numeric)
  kpcooc=matrix(0,ncol = ncol(results.kp),nrow = ncol(results.kp))
  for(i1 in 1:ncol(results.kp))
    for(i2 in 1:ncol(results.kp)){
      kpcooc[i1,i2]=length(which(results.kp[,i1]==results.kp[,i2]))
    }
  numcl=round(max(apply(results.kp,1,function(x) length(unique(x)))))
  for(thr in seq(0,ntimes,ntimes/100)){
    kphc=hclust(as.dist(ntimes-kpcooc),method = 'average')
    kp=cutree(kphc,h=thr)
    if(max(kp)<numcl)
      break
  }
  binarymat=matrix(0,nrow = length(kp),ncol = length(kp))
  for(i in 1:numcl){
    binarymat[which(kp==i),which(kp==i)]=1
  }
  print("Cophenetic Coefficient for Column Clustering:")
  print(abs(cor(as.vector(binarymat),as.vector(kpcooc))))
  
  ## Row Clusters
  results.kt=apply(results.kt,2,as.numeric)
  ktcooc=matrix(0,ncol = ncol(results.kt),nrow = ncol(results.kt))
  for(i1 in 1:ncol(results.kt))
    for(i2 in 1:ncol(results.kt)){
      ktcooc[i1,i2]=length(which(results.kt[,i1]==results.kt[,i2]))
    }
  numcl=round(max(apply(results.kt,1,function(x) length(unique(x)))))
  for(thr in seq(0,ntimes,ntimes/100)){
    kthc=hclust(as.dist(ntimes-ktcooc),method = 'average')
    kt=cutree(kthc,h=thr)
    if(max(kt)<numcl)
      break
  }
  binarymat=matrix(0,nrow = length(kt),ncol = length(kt))
  for(i in 1:numcl){
    binarymat[which(kt==i),which(kt==i)]=1
  }
  print("Cophenetic Coefficient for Row Clustering:")
  print(abs(cor(as.vector(binarymat),as.vector(ktcooc))))
  
  ## Row Weights
  weights=apply(results.w,2,mean)
  
  list(column.clusters=kp,row.clusters=kt,row.weights=weights)
}



# This function draws a heatmap for a given biclustering
# The output is a heatmap, a visualization of the Theta and other model parameters, with rows corresponding to row clusters and columns to column clusters
# The heatmap consists of the following elements:
#    LR: Label Ratio--the ratio of '1' labels inside the column clusters;
#    RCSize/CCSize: the size of row/column cluster
#    AvgW: average weight for row cluster
# Row clusters are sorted by 'AvgW' in descending order from top to bottom
# The color inside each cell indicates the proportion of 1 data points for the elements inside the corresponding bicluster, with red meaning high proportion and blue indicating low proportion

SUBSTRA.heatmap <- function(data, labels,
                            model,      # the model for which the heatmap is being drawn
                            pathToPDF,  # the address to store the heatmap
                            title){     # the title of the heatmap
  
  column.clusters=model$column.clusters
  row.clusters=model$row.clusters
  row.weights=model$row.weights
  
  ncc=max(column.clusters)
  nrc=max(row.clusters)
  data=as.matrix(data)
  maxv=max(data)
  Theta=matrix(0,nrow=nrc,ncol=ncc)
  for(i in 1:nrc){
    tinds=which(row.clusters==i)
    for(j in 1:ncc){
      pinds=which(column.clusters==j)
      Theta[i,j]=(length(which(data[tinds,pinds]==maxv)))/(length(pinds)*length(tinds))
    }
  }

  colnames(Theta)=c(1:ncc)
  rownames(Theta)=c(1:nrc)

  df=data.frame(RCl=row.clusters,W=row.weights)
  temp=aggregate(df$W,by=list(RCl=df$RCl),FUN=function(x) c(AvgW=mean(x),Size=length(x)))
  rowanot=data.frame(RCSize=temp$x[,2],
                     AvgW=temp$x[,1])
  rownames(rowanot)=temp$RCl
  roworder=order(rowanot$AvgW,decreasing = T)

  df=data.frame(CCl=column.clusters,Labels=labels)
  temp=aggregate(df$Labels,by=list(CCl=df$CCl),FUN=function(x) c(mean(x),length(x)))
  colanot=data.frame(CCSize=temp$x[,2],LR=temp$x[,1])
  rownames(colanot)=temp$CCl
  colorder=order(colanot$LR)


  pheatmap(Theta[match(rownames(rowanot[roworder,]),rownames(Theta)),match(rownames(colanot[colorder,]),colnames(Theta))],
           annotation_row = rowanot,annotation_col = colanot,#fontsize_row = 20,fontsize_col = 20,
           fontsize = 13,show_colnames = T,show_rownames = T, filename = pathToPDF,width = 11.69,height = 8.27,
           cluster_cols = F,cluster_rows = F,main = title,legend_breaks = seq(0,1,by = 0.1),gaps_col = length(which(colanot$LR==0)))
}
