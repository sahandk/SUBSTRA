library(pROC)
library(clusteval)



SUBSTRA.split_dataset = function(data, labels, folds){
  cases=which(labels==1)
  controls=which(labels==0)
  casestops=round(c(0,c(1:folds)*length(cases)/folds))
  controlstops=round(c(0,c(1:folds)*length(controls)/folds))
  dataparts=list()
  labelparts=list()
  for(i in 1:folds){
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



SUBSTRA.train <- function(expresData, labels, sideData=NULL,Kp, Kt, phase_1_ite, phase_2_ite=50,
                         alpha_p = 1, alpha_t = 1, alpha_ph = 1e-100, alpha_e = 1, alpha_sd = 1, magp=1, return.best=T){
  if(min(min(expresData))==0) expresData = expresData + 1
  withLabel=!is.null(labels)
  if(withLabel && min(labels)==0) labels = labels + 1
  withSideData=!is.null(sideData)
  if(withSideData && min(min(sideData))==0) sideData = sideData + 1
  Nv = length(unique(labels))
  Np = ncol(expresData)
  Nt = nrow(expresData)
  logKp=log(Kp)
  logKt=log(Kt)
  logDBLMax=log(.Machine$double.xmax)
  # Clustering
  kp = c()
  kt = c()
  npINkp = c()
  ntINkt = c()
  # Bicluster Counters & Parameters
  nesINkpkt = array(0,dim = c(Kp, Kt, 2))
  neINpFORkt = array(0,dim = c(Np, Kt, 2))
  neINtFORkp = array(0,dim = c(Nt, Kp, 2))
  trans.weights=array(magp,dim = c(1, Nt))
  # Side expresData Counters
  nvsINktkt = array(0,dim = c(Kt, Kt, 2))
  nvINtFORkt = array(0,dim = c(Nt, Kt, 2))
  # Labels Counters
  nvsINKp =  array(0,dim = c(Kp, Nv))
  
  # Random Initialization
  for(i in 1:Np)
    kp[i] = sample(seq(labels[i],Kp,by=2),1)
  kt = sample(c(1:Kt),Nt,replace = T)
  for(i in 1:Kp){
    kpinds=which(kp==i)
    npINkp[i]=length(kpinds)
    if(npINkp[i]!=0){
      if(withLabel){
        for(v in 1:Nv) nvsINKp[i,v]=length(which(labels[kpinds]==v))
      }
      
      for(j in 1:Kt){
        ktinds=which(kt==j)
        if(length(ktinds)==0) next
        nesINkpkt[i,j,1]=sum(2-expresData[ktinds,kpinds])
        nesINkpkt[i,j,2]=sum(expresData[ktinds,kpinds]-1)
      }
      if(npINkp[i]==1){
        neINtFORkp[,i,1]= (2-expresData[,kpinds])
        neINtFORkp[,i,2]= (expresData[,kpinds]-1)
      }
      else{
        neINtFORkp[,i,1]= rowSums(2-expresData[,kpinds])
        neINtFORkp[,i,2]= rowSums(expresData[,kpinds]-1)
      }
    }
  }
  for(j in 1:Kt){
    ktinds=which(kt==j)
    ntINkt[j]=length(ktinds)
    if(ntINkt[j]!=0){
      if(withSideData){
        for(jj in 1:Kt){
          ntinkt2=which(kt==jj)
          nvsINktkt[j,jj,1]=sum(2-sideData[ktinds,ktinds2])
          nvsINktkt[j,jj,2]=sum(sideData[ktinds,ktinds2]-1)
        }
        
        if(ntINkt[j]==1){
          nvINtFORkt[,j,1]=(2-sideData[,ktinds])
          nvINtFORkt[,j,2]=(sideData[,ktinds]-1)
        }else{
          nvINtFORkt[,j,1]=rowSums(2-sideData[,ktinds])
          nvINtFORkt[,j,2]=rowSums(sideData[,ktinds]-1)
        }
      }
      if(ntINkt[j]==1){
        neINpFORkt[,j,1]=(2-expresData[ktinds,])
        neINpFORkt[,j,2]=(expresData[ktinds,]-1)  
      }else{
        neINpFORkt[,j,1]=colSums(2-expresData[ktinds,])
        neINpFORkt[,j,2]=colSums(expresData[ktinds,]-1)
      }
    }
  }
  
  if(length(which(nesINkpkt<0))>0)
    print("initialization")
  
  Pi=(npINkp+alpha_p/Kp)/(Np+alpha_p)
  Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
  Psi=(nvsINKp+alpha_ph/Nv)/array(rep(npINkp+alpha_ph,2),dim=dim(nvsINKp))
  # Psi=nvsINKp/array(rep(rowSums(nvsINKp),2),dim=dim(nvsINKp))
  # Psi[is.nan(Psi)]=0.5
  best.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
                  bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
                  Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights,
                  sd.bicluster.counts=nvsINktkt, sd.counts.for.trans=nvINtFORkt)
  preds=SUBSTRA.predict(best.model,expresData)
  # print(preds)
  r=roc(as.factor(labels),preds[,2],direction = "<")
  bestauc=r$auc
  besterr=mean((labels-1-preds[,2])^2)
  print("Initial Status:")
  print(bestauc)
  print(besterr)
  
  # Sampling
  iterations=phase_1_ite+phase_2_ite
  for(ite in 1:iterations){
    model.prob = 0
    prev.kp=kp
    prev.kt=kt
    
    ## Transcript Clusters
    t.diff = 0
    for(t in 1:Nt){
      ### take out
      tempkt=kt[t]
      ntINkt[tempkt]=ntINkt[tempkt]-1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] - neINtFORkp[t,,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] - (2-expresData[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] - (expresData[t,]-1)
      # for(p in 1:Np){
      #   neINpFORkt[p,tempkt,expresData[t,p]]=neINpFORkt[p,tempkt,expresData[t,p]]-1
      # }
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]-nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]-nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]-ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]-ifelse(sideData[tt,t]==2,1,0)
        }
      }
      
      if(length(which(nesINkpkt<0))>0)
        print(t)
      
      ### compute probability
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      t.prob=log(ntINkt + alpha_t / Kt)+
        colSums(log(Theta[,,1])*matrix(rep(neINtFORkp[t,,1],Kt),nrow = Kp,ncol = Kt,byrow = F))+
        colSums(log(Theta[,,2])*matrix(rep(neINtFORkp[t,,2],Kt),nrow = Kp,ncol = Kt,byrow = F))
      
      # t.prob.max=-1*.Machine$double.xmax
      # t.prob=c()
      # for(cl in 1:Kt){
      #   t.prob[cl]=log(ntINkt[cl] + alpha_t / Kt)
      #   t.prob[cl]=t.prob[cl]+sum(log((nesINkpkt[,cl,1]+alpha_e/2)/(ntINkt[cl]*npINkp+alpha_e))*neINtFORkp[t,,1])+
      #     sum(log((nesINkpkt[,cl,2]+alpha_e/2)/(ntINkt[cl]*npINkp+alpha_e))*neINtFORkp[t,,2])
      #   if(withSideData){
      #     t.prob[cl]=t.prob[cl]+sum(log((nesINktkt[,cl,1]+alpha_sd/2)/(ntINkt[cl]*ntINkt+alpha_sd))*nvINtFORkt[t,,1])+
      #       sum(log((nesINktkt[,cl,2]+alpha_sd/2)/(ntINkt[cl]*ntINkt+alpha_sd))*nvINtFORkt[t,,2])
      #   }
      #   t.prob[cl]=t.prob[cl]*magt
      #   if(t.prob[cl]>t.prob.max)
      #     t.prob.max=t.prob[cl]
      # }
      
      ### put back
      t.prob.max=max(t.prob)
      t.prob.max = t.prob.max - (logDBLMax - logKt)
      t.prob = exp(t.prob - t.prob.max)
      total.prob = sum(t.prob)
      tempkt=sample(1:Kt,1,prob = t.prob/total.prob)
      if (kt[t] != tempkt)
        t.diff=t.diff+1
      model.prob = model.prob + log(t.prob[tempkt] / total.prob)
      
      kt[t]=tempkt
      ntINkt[tempkt]=ntINkt[tempkt]+1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] + neINtFORkp[t,,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] + (2-expresData[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] + (expresData[t,]-1)
      # for(p in 1:Np){
      #   neINpFORkt[p,tempkt,expresData[t,p]]=neINpFORkt[p,tempkt,expresData[t,p]]+1
      # }
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]+nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]+nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]+ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]+ifelse(sideData[tt,t]==2,1,0)
        }
      }
    }#t
    
    if(length(which(nesINkpkt<0))>0)
      print("transcripts")
    
    # trans.weights=rep(magp,Nt)
    ## Patient Clusters
    p.diff = 0
    counter=0
    for(p in 1:Np){
      ### take out
      tempkp=kp[p]
      npINkp[tempkp]=npINkp[tempkp]-1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] - neINpFORkt[p,,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] - (2-expresData[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] - (expresData[,p]-1)
      # for(t in 1:Nt)
      #   neINtFORkp[t,tempkp,expresData[t,p]]=neINtFORkp[t,tempkp,expresData[t,p]]-1
      if (withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]-1
      
      ### compute weights
      if(ite > phase_1_ite){
        flag=F
        best.weights=trans.weights
        Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
        Pi=(npINkp+alpha_p/Kp)/(Np-1+alpha_p)
        Psi=(nvsINKp+alpha_ph/Nv)/array(rep(npINkp+alpha_ph,2),dim=dim(nvsINKp))
        
        BO.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
                      bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
                      Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights,
                      sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
        preds=SUBSTRA.predict(BO.model,as.matrix(expresData[,p]))
        AOerr=(labels[p]-1-preds[2])^2
        
        reverseExp=matrix(rep(2-expresData[,p],Kp),nrow = Kp,ncol = Nt,byrow = T)
        logTheta=log(abs(reverseExp-Theta[,kt,2]))
        
        for(curnu in c(magp)){
          # print(sprintf("rep = %d",rep))
          temp.weights=trans.weights
          # temp.weights=rep(magp,Nt)/ntINkt[kt]
          # temp.weights=rep(1,Nt)
          temp.weights.matrix=matrix(rep(temp.weights,Kp),nrow = Kp,ncol = Nt,byrow = T)
          # temp.weights=rep(1,Nt)
          # featno=2000
          
          # inds=c(1:featno)
          # for(lii in 1:ceiling(Nt/featno)){
          #   ti=inds+featno*(lii-1)
          #   if(ti[featno]>Nt)
          #     ti=c(ti[1]:Nt)
          remaining=c(1:Nt)
          while(length(remaining) > 0){
            # ti=sample(remaining,min(featno,length(remaining)))
            ti=remaining
            remaining=setdiff(remaining,ti)
            # reverseExp=matrix(rep(2-expresData[ti,p],Kp),nrow = Kp,ncol = length(ti),byrow = T)
            # trans.weights.matrix=matrix(rep(trans.weights[ti],Kp),nrow = Kp,ncol = length(ti),byrow = T)
            
            logs=rowSums(logTheta[,ti]*temp.weights.matrix[,ti])+log(Pi)
            logs=logs+log(.Machine$double.xmax)-max(logs)-log(length(logs))
            values=exp(logs)
            values=values/sum(values)
            # print(paste0("values: ", values))
            A=sum(values*Psi[,2])
            # print(paste0("A: ",A))
            B=1
            # dA=colSums(matrix(rep(values*Psi[,2],Nt),byrow = F,nrow = Kp,ncol = Nt)*logTheta)
            dA=colSums(matrix(rep(values*Psi[,2],length(ti)),byrow = F,nrow = Kp,ncol = length(ti))*logTheta[,ti])
            # print(paste0("dA: ",min(dA)," - ",max(dA)))
            # dB=colSums(matrix(rep(values,Nt),byrow = F,nrow = Kp,ncol = Nt)*logTheta)
            dB=colSums(matrix(rep(values,length(ti)),byrow = F,nrow = Kp,ncol = length(ti))*logTheta[,ti])
            # print(paste0("dB: ",min(dB)," - ",max(dB)))
            f=A/B
            df=(dA*B-dB*A)/(B*B)
            # print(paste0("df: ",min(df)," - ",max(df)))
            # trans.weights=trans.weights+nu*2*(labels[p]-1-f)*df
            temp.weights[ti]=temp.weights[ti]+curnu*2*(labels[p]-1-f)*df
            # print("trans.weights:")
            # print(paste0(min(trans.weights)," - ",max(trans.weights)))
            
          }#lii
          temp.weights[which(temp.weights<0)]=0
          
          AO.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
                        bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
                        Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=temp.weights,
                        sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
          preds=SUBSTRA.predict(AO.model,as.matrix(expresData[,p]))
          curerr=(labels[p]-1-preds[2])^2
          
          # preds=SUBSTRA.predict(AO.model,as.matrix(expresData))
          # curerr=mean((labels-1-preds[,2])^2)
          # print(sprintf("%d -> %f", p, curerr))
          
          # if(is.nan(curerr)){
          #   next
          #   # print("not changed")
          # }else 
          if(curerr < AOerr){
            best.weights=temp.weights
            AOerr=curerr
            flag=T
          }
        }#curnu
        
        # print(sprintf("%d -> err: %f to %f - w in [%f, %f] - distance: %g",p,BOerr,AOerr,min(best.weights),max(best.weights),mean(abs(trans.weights-best.weights))))
        if(flag)
          counter=counter+1
        
        trans.weights=best.weights
      }
      
      # if(ite > burn_in){
      #   # ran=rank(best.weights,ties.method = "random")
      #   # best.weights[ran>Nt-100]=1
      #   # best.weights[ran<=Nt-100]=0
      #   best.weights[best.weights<0]=0
      #   trans.weights=trans.weights+best.weights
      # }
      
      ### compute probability
      reverseExp=matrix(rep(2-expresData[,p],Kp),nrow = Kp,ncol = Nt,byrow = T)
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      p.prob=rowSums(log(abs(reverseExp-Theta[,kt,2]))*
                       matrix(rep(trans.weights,Kp),nrow = Kp,ncol = Nt,byrow = T))+
        log(npINkp + alpha_p / Kp)
      # p.prob=log(npINkp + alpha_p / Kp)+
      #   rowSums(log(Theta[,,1])*matrix(rep(neINpFORkt[p,,1],Kp),nrow = Kp,ncol = Kt,byrow = T))+
      #   rowSums(log(Theta[,,2])*matrix(rep(neINpFORkt[p,,2],Kp),nrow = Kp,ncol = Kt,byrow = T))
      p.prob.max=max(p.prob[seq(labels[p],Kp,by=2)])
      
      ### put back
      p.prob.max = p.prob.max - (logDBLMax - logKp)
      p.prob = exp(p.prob - p.prob.max)
      # print(sprintf("label: %d",labels[p]))
      p.prob[seq(3-labels[p],Kp,by=2)]=0
      # print(p.prob/sum(p.prob))
      total.prob = sum(p.prob)
      tempkp=sample(1:Kp,1,prob = p.prob/total.prob)
      if (kp[p] != tempkp)
        p.diff=p.diff+1
      model.prob = model.prob + log(p.prob[tempkp] / total.prob)
      
      kp[p]=tempkp
      npINkp[tempkp]=npINkp[tempkp]+1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] + neINpFORkt[p,,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] + (2-expresData[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] + (expresData[,p]-1)
      # for(t in 1:Nt)
      #   neINtFORkp[t,tempkp,expresData[t,p]]=neINtFORkp[t,tempkp,expresData[t,p]]+1
      if (withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]+1
    }#p
    
    if(length(which(nesINkpkt<0))>0)
      print("patients")
    
    print(paste0("Iteration: ",ite))
    print("Difference:")
    print(paste0("p: ",p.diff," - ",length(which(npINkp>0))," of ",Kp," clusters occupied."))
    print(paste0("t: ",t.diff," - ",length(which(ntINkt>0))," of ",Kt," clusters occupied."))
    print(paste0("w: ",counter))
    rkp=cluster_similarity(kp,prev.kp)
    rkt=cluster_similarity(kt,prev.kt)
    print("Rand Indexes:")
    print(paste0("Patient Clustering:    ",rkp))
    print(paste0("Transcript Clustering: ",rkt))
    
    if(rkp > 0.90 && rkt > 0.90)
      ite=phase_1_ite
    
    # trans.weights=trans.weights/Np
    # trans.weights[which(trans.weights<1)]=0
    # trans.weights=Nt/(Kt*ntINkt[kt])
    # trans.weights[which(kt %in% which(kt.weights==0))]=0  
    # trans.weights=rep(magp,Nt)/ntINkt[kt]
    print(paste0(min(trans.weights), " - ", max(trans.weights), " : ", mean(trans.weights)))
    Pi=(npINkp+alpha_p/Kp)/(Np+alpha_p)
    Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
    Psi=(nvsINKp+alpha_ph/Nv)/array(rep(npINkp+alpha_ph,2),dim=dim(nvsINKp))
    # Psi=nvsINKp/array(rep(rowSums(nvsINKp),2),dim=dim(nvsINKp))
    # Psi[is.nan(Psi)]=0.5
    model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
               bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
               Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights,
               sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
    preds=SUBSTRA.predict(model,expresData)
    r=roc(as.factor(labels),preds[,2],direction = "<")
    print(r$auc)
    err=mean((labels-1-preds[,2])^2)
    print(err)
    
    # if(err <= besterr && r$auc >= bestauc){
    if(r$auc > bestauc || ite == phase_1_ite){
      best.model=model
      besterr=err
      bestauc=r$auc
    }else if(r$auc == bestauc){
      if(err <= besterr){
        best.model=model
        besterr=err
        bestauc=r$auc
      }
    }
  }#ite
  
  if(return.best)
    return(best.model)
  else
    return(model)
}



SUBSTRA.train_np <- function(expresData, labels, sideData=NULL, phase_1_ite, phase_2_ite=50, ubKp=1000,ubKt=2000,
                            alpha_p = 1, alpha_t = 1, alpha_ph = 1e-100, alpha_e = 1, alpha_sd = 1, magnitude=1, return.best=T){
  if(min(min(expresData))==0) expresData = expresData + 1
  withLabel=!is.null(labels)
  if(withLabel && min(labels)==0) labels = labels + 1
  withSideData=!is.null(sideData)
  if(withSideData && min(min(sideData))==0) sideData = sideData + 1
  Nv = length(unique(labels))
  Np = ncol(expresData)
  Kp=min(Np,ubKp)
  Nt = nrow(expresData)
  Kt=min(Nt,ubKt)
  logKp=log(Kp)
  logKt=log(Kt)
  logDBLMax=log(.Machine$double.xmax)
  # Clustering
  kp = c()
  kt = c()
  npINkp = c()
  ntINkt = c()
  # Bicluster Counters & Parameters
  nesINkpkt = array(0,dim = c(Kp, Kt, 2))
  neINpFORkt = array(0,dim = c(Np, Kt, 2))
  neINtFORkp = array(0,dim = c(Nt, Kp, 2))
  trans.weights=array(magnitude,dim = c(1, Nt))
  # Side expresData Counters
  nvsINktkt = array(0,dim = c(Kt, Kt, 2))
  nvINtFORkt = array(0,dim = c(Nt, Kt, 2))
  # Labels Counters
  nvsINKp =  array(0,dim = c(Kp, Nv))
  
  # Random Initialization
  for(i in 1:Np)
    kp[i] = sample(seq(labels[i],Np/5,by=2),1)
  kt = sample(c(1:(Nt/50)),Nt,replace = T)
  for(i in 1:Kp){
    kpinds=which(kp==i)
    npINkp[i]=length(kpinds)
    if(npINkp[i]!=0){
      if(withLabel){
        for(v in 1:Nv) nvsINKp[i,v]=length(which(labels[kpinds]==v))
      }
      
      for(j in 1:Kt){
        ktinds=which(kt==j)
        if(length(ktinds)==0) next
        nesINkpkt[i,j,1]=sum(2-expresData[ktinds,kpinds])
        nesINkpkt[i,j,2]=sum(expresData[ktinds,kpinds]-1)
      }
      if(npINkp[i]==1){
        neINtFORkp[,i,1]= (2-expresData[,kpinds])
        neINtFORkp[,i,2]= (expresData[,kpinds]-1)
      }
      else{
        neINtFORkp[,i,1]= rowSums(2-expresData[,kpinds])
        neINtFORkp[,i,2]= rowSums(expresData[,kpinds]-1)
      }
    }
  }
  for(j in 1:Kt){
    ktinds=which(kt==j)
    ntINkt[j]=length(ktinds)
    if(ntINkt[j]!=0){
      if(withSideData){
        for(jj in 1:Kt){
          ntinkt2=which(kt==jj)
          nvsINktkt[j,jj,1]=sum(2-sideData[ktinds,ktinds2])
          nvsINktkt[j,jj,2]=sum(sideData[ktinds,ktinds2]-1)
        }
        
        if(ntINkt[j]==1){
          nvINtFORkt[,j,1]=(2-sideData[,ktinds])
          nvINtFORkt[,j,2]=(sideData[,ktinds]-1)
        }else{
          nvINtFORkt[,j,1]=rowSums(2-sideData[,ktinds])
          nvINtFORkt[,j,2]=rowSums(sideData[,ktinds]-1)
        }
      }
      if(ntINkt[j]==1){
        neINpFORkt[,j,1]=(2-expresData[ktinds,])
        neINpFORkt[,j,2]=(expresData[ktinds,]-1)  
      }else{
        neINpFORkt[,j,1]=colSums(2-expresData[ktinds,])
        neINpFORkt[,j,2]=colSums(expresData[ktinds,]-1)
      }
    }
  }
  
  if(length(which(nesINkpkt<0))>0)
    print("initialization")
  
  occupied.kp=unique(kp)
  occupied.kt=unique(kt)
  empty.kp=c()
  empty.kt=c()
  if(length(occupied.kp) < Kp)
    empty.kp=min(which(npINkp==0))
  if(length(occupied.kt) < Kt)
    empty.kt=min(which(ntINkt==0))
  
  Pi=npINkp[occupied.kp]/Np
  Theta=(nesINkpkt[occupied.kp,occupied.kt,]+alpha_e/2)/array(rep(nesINkpkt[occupied.kp,occupied.kt,1]+nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e,2),dim = c(length(occupied.kp),length(occupied.kt),2))
  Psi=(nvsINKp[occupied.kp,]+alpha_ph/Nv)/array(rep(npINkp[occupied.kp]+alpha_ph,2),dim=c(length(occupied.kp),2))
  best.model=list(patient.clusters=kp, transcript.clusters=match(kt,occupied.kt),
                  bicluster.content=nesINkpkt[occupied.kp,occupied.kt,], counts.for.trans=neINtFORkp[,occupied.kp,], counts.for.pats=neINpFORkt[,occupied.kt,], counts.for.pheno=nvsINKp[occupied.kp,],
                  Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights)
  preds=SUBSTRA.predict(best.model,expresData)
  # print(preds)
  r=roc(as.factor(labels),preds[,2],direction = "<")
  bestauc=r$auc
  besterr=mean((labels-1-preds[,2])^2)
  print("Initial Status:")
  print(bestauc)
  print(besterr)
  
  # Sampling
  
  ################### PHASE I ###################
  
  for(ite in 1:phase_1_ite){
    model.prob = 0
    prev.kp=kp
    prev.kt=kt
    
    ## Transcript Clusters
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
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] - (2-expresData[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] - (expresData[t,]-1)
      # for(p in 1:Np){
      #   neINpFORkt[p,tempkt,expresData[t,p]]=neINpFORkt[p,tempkt,expresData[t,p]]-1
      # }
      # if (withSideData){
      #   nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]-nvINtFORkt[tempkt,]
      #   nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]-nvINtFORkt[tempkt,]
      #   nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
      #   for(tt in 1:Nt)
      #   {
      #     nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]-ifelse(sideData[tt,t]==1,1,0)
      #     nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]-ifelse(sideData[tt,t]==2,1,0)
      #   }
      # }
      
      if(length(which(nesINkpkt<0))>0)
        print(t)
      
      sig.kp=c(occupied.kp,empty.kp)
      sig.kt=c(occupied.kt,empty.kt)
      ### compute probability
      Theta=(nesINkpkt[sig.kp,sig.kt,]+alpha_e/2)/array(rep(nesINkpkt[sig.kp,sig.kt,1]+nesINkpkt[sig.kp,sig.kt,2]+alpha_e,2),dim = c(length(sig.kp),length(sig.kt),2))
      t.prob=colSums(log(Theta[,,1])*matrix(rep(neINtFORkp[t,sig.kp,1],length(sig.kt)),nrow = length(sig.kp),ncol = length(sig.kt),byrow = F))+
        colSums(log(Theta[,,2])*matrix(rep(neINtFORkp[t,sig.kp,2],length(sig.kt)),nrow = length(sig.kp),ncol = length(sig.kt),byrow = F))
      
      if(is.null(empty.kt))
        t.prob=t.prob+log(ntINkt[occupied.kt])
      else
        t.prob=t.prob+log(c(ntINkt[occupied.kt],alpha_t))
      # t.prob.max=-1*.Machine$double.xmax
      # t.prob=c()
      # for(cl in 1:Kt){
      #   t.prob[cl]=log(ntINkt[cl] + alpha_t / Kt)
      #   t.prob[cl]=t.prob[cl]+sum(log((nesINkpkt[,cl,1]+alpha_e/2)/(ntINkt[cl]*npINkp+alpha_e))*neINtFORkp[t,,1])+
      #     sum(log((nesINkpkt[,cl,2]+alpha_e/2)/(ntINkt[cl]*npINkp+alpha_e))*neINtFORkp[t,,2])
      #   if(withSideData){
      #     t.prob[cl]=t.prob[cl]+sum(log((nesINktkt[,cl,1]+alpha_sd/2)/(ntINkt[cl]*ntINkt+alpha_sd))*nvINtFORkt[t,,1])+
      #       sum(log((nesINktkt[,cl,2]+alpha_sd/2)/(ntINkt[cl]*ntINkt+alpha_sd))*nvINtFORkt[t,,2])
      #   }
      #   t.prob[cl]=t.prob[cl]*magt
      #   if(t.prob[cl]>t.prob.max)
      #     t.prob.max=t.prob[cl]
      # }
      
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
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] + (2-expresData[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] + (expresData[t,]-1)
      
      if(!is.null(empty.kt) && tempkt==empty.kt){
        occupied.kt=unique(kt)
        empty.kt=c()
        if(length(occupied.kt) < Kt)
          empty.kt=min(which(ntINkt==0))
      }
      # for(p in 1:Np){
      #   neINpFORkt[p,tempkt,expresData[t,p]]=neINpFORkt[p,tempkt,expresData[t,p]]+1
      # }
      # if (withSideData){
      #   nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]+nvINtFORkt[tempkt,]
      #   nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]+nvINtFORkt[tempkt,]
      #   nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
      #   for(tt in 1:Nt)
      #   {
      #     nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]+ifelse(sideData[tt,t]==1,1,0)
      #     nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]+ifelse(sideData[tt,t]==2,1,0)
      #   }
      # }
    }#t
    
    if(length(which(nesINkpkt<0))>0)
      print("transcripts")
    
    # trans.weights=rep(magnitude,Nt)
    ## Patient Clusters
    p.diff = 0
    counter=0
    for(p in 1:Np){
      ### take out
      tempkp=kp[p]
      npINkp[tempkp]=npINkp[tempkp]-1
      if(npINkp[tempkp]==0){
        occupied.kp=setdiff(occupied.kp,tempkp)
        empty.kp=tempkp
      }
      nesINkpkt[tempkp,occupied.kt,] = nesINkpkt[tempkp,occupied.kt,] - neINpFORkt[p,occupied.kt,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] - (2-expresData[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] - (expresData[,p]-1)
      
      if(withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]-1
      
      ### compute probability
      sig.kp=c(occupied.kp,empty.kp)
      sig.kt=c(occupied.kt,empty.kt)
      reverseExp=matrix(rep(2-expresData[,p],length(sig.kp)),nrow = length(sig.kp),ncol = Nt,byrow = T)
      Theta=(nesINkpkt[sig.kp,sig.kt,]+alpha_e/2)/array(rep(nesINkpkt[sig.kp,sig.kt,1]+nesINkpkt[sig.kp,sig.kt,2]+alpha_e,2),dim = c(length(sig.kp),length(sig.kt),2))
      p.prob=rowSums(log(abs(reverseExp-Theta[,match(kt,sig.kt),2]))*
                       matrix(rep(trans.weights,length(sig.kp)),nrow = length(sig.kp),ncol = Nt,byrow = T))+
        log((nvsINKp[sig.kp,labels[p]]+alpha_ph/Nv)/(npINkp[sig.kp]+alpha_ph))
      
      if(is.null(empty.kp))
        p.prob=p.prob+log(npINkp[occupied.kp])
      else
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
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] + (2-expresData[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] + (expresData[,p]-1)
      
      if (withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]+1
      
      if(!is.null(empty.kp) && tempkp==empty.kp){
        occupied.kp=unique(kp)
        empty.kp=c()
        if(length(occupied.kp) < Kp)
          empty.kp=min(which(npINkp==0))
      }
    }#p
    
    if(length(which(nesINkpkt<0))>0)
      print("patients")
    
    print(paste0("Phase I - Iteration: ",ite))
    print("Difference:")
    print(paste0("p: ",p.diff," - ",length(occupied.kp)," clusters occupied."))
    print(paste0("t: ",t.diff," - ",length(occupied.kt)," clusters occupied."))
    print(paste0("w: ",counter))
    rkp=cluster_similarity(kp,prev.kp)
    rkt=cluster_similarity(kt,prev.kt)
    print("Rand Indexes:")
    print(paste0("Patient Clustering:    ",rkp))
    print(paste0("Transcript Clustering: ",rkt))
    
    print(paste0(min(trans.weights), " - ", max(trans.weights), " : ", mean(trans.weights)))
    Pi=npINkp[occupied.kp]/Np
    Theta=(nesINkpkt[occupied.kp,occupied.kt,]+alpha_e/2)/array(rep(nesINkpkt[occupied.kp,occupied.kt,1]+nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e,2),dim = c(length(occupied.kp),length(occupied.kt),2))
    Psi=(nvsINKp[occupied.kp,]+alpha_ph/Nv)/array(rep(npINkp[occupied.kp]+alpha_ph,2),dim=c(length(occupied.kp),2))
    model=list(patient.clusters=kp, transcript.clusters=match(kt,occupied.kt),
               bicluster.content=nesINkpkt[occupied.kp,occupied.kt,], counts.for.trans=neINtFORkp[,occupied.kp,], counts.for.pats=neINpFORkt[,occupied.kt,], counts.for.pheno=nvsINKp[occupied.kp,],
               Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights)
    
    preds=SUBSTRA.predict(model,expresData)
    r=roc(as.factor(labels),preds[,2],direction = "<")
    print(r$auc)
    err=mean((labels-1-preds[,2])^2)
    print(err)
    
    if(rkp > 0.95 && rkt > 0.95)
      break
    
    # # if(err <= besterr && r$auc >= bestauc){
    # if(r$auc > bestauc || ite == phase_1_ite){
    #   best.model=model
    #   besterr=err
    #   bestauc=r$auc
    # }else if(r$auc == bestauc){
    #   if(err <= besterr){
    #     best.model=model
    #     besterr=err
    #     bestauc=r$auc
    #   }
    # }
  }#phase_1_ite
  
  
  
  ################### PHASE II ###################
  
  for(ite in 1:phase_2_ite){
    model.prob = 0
    prev.kp=kp
    prev.kt=kt
    
    ## Transcript Clusters
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
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] - (2-expresData[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] - (expresData[t,]-1)
      
      # if (withSideData){
      #   nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]-nvINtFORkt[tempkt,]
      #   nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]-nvINtFORkt[tempkt,]
      #   nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
      #   for(tt in 1:Nt)
      #   {
      #     nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]-ifelse(sideData[tt,t]==1,1,0)
      #     nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]-ifelse(sideData[tt,t]==2,1,0)
      #   }
      # }
      
      if(length(which(nesINkpkt<0))>0)
        print(t)
      
      sig.kp=c(occupied.kp,empty.kp)
      sig.kt=c(occupied.kt,empty.kt)
      ### compute probability
      Theta=(nesINkpkt[sig.kp,sig.kt,]+alpha_e/2)/array(rep(nesINkpkt[sig.kp,sig.kt,1]+nesINkpkt[sig.kp,sig.kt,2]+alpha_e,2),dim = c(length(sig.kp),length(sig.kt),2))
      t.prob=colSums(log(Theta[,,1])*matrix(rep(neINtFORkp[t,sig.kp,1],length(sig.kt)),nrow = length(sig.kp),ncol = length(sig.kt),byrow = F))+
        colSums(log(Theta[,,2])*matrix(rep(neINtFORkp[t,sig.kp,2],length(sig.kt)),nrow = length(sig.kp),ncol = length(sig.kt),byrow = F))
      
      if(is.null(empty.kt))
        t.prob=t.prob+log(ntINkt[occupied.kt])
      else
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
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] + (2-expresData[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] + (expresData[t,]-1)
      
      if(!is.null(empty.kt) && tempkt==empty.kt){
        occupied.kt=unique(kt)
        empty.kt=c()
        if(length(occupied.kt) < Kt)
          empty.kt=min(which(ntINkt==0))
      }
      
      # if (withSideData){
      #   nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]+nvINtFORkt[tempkt,]
      #   nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]+nvINtFORkt[tempkt,]
      #   nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
      #   for(tt in 1:Nt)
      #   {
      #     nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]+ifelse(sideData[tt,t]==1,1,0)
      #     nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]+ifelse(sideData[tt,t]==2,1,0)
      #   }
      # }
    }#t
    
    if(length(which(nesINkpkt<0))>0)
      print("transcripts")
    
    # trans.weights=rep(magnitude,Nt)
    ## Patient Clusters
    p.diff = 0
    counter=0
    for(p in 1:Np){
      ### take out
      tempkp=kp[p]
      npINkp[tempkp]=npINkp[tempkp]-1
      if(npINkp[tempkp]==0){
        occupied.kp=setdiff(occupied.kp,tempkp)
        empty.kp=tempkp
      }
      nesINkpkt[tempkp,occupied.kt,] = nesINkpkt[tempkp,occupied.kt,] - neINpFORkt[p,occupied.kt,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] - (2-expresData[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] - (expresData[,p]-1)
      
      if(withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]-1
      
      ### compute weights
      flag=F
      best.weights=trans.weights
      Theta=(nesINkpkt[occupied.kp,occupied.kt,]+alpha_e/2)/array(rep(nesINkpkt[occupied.kp,occupied.kt,1]+nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e,2),dim = c(length(occupied.kp),length(occupied.kt),2))
      Pi=npINkp[occupied.kp]/Np
      Psi=(nvsINKp[occupied.kp,]+alpha_ph/Nv)/array(rep(npINkp[occupied.kp]+alpha_ph,2),dim=c(length(occupied.kp),2))
      
      BO.model=list(patient.clusters=kp, transcript.clusters=match(kt,occupied.kt),
                    bicluster.content=nesINkpkt[occupied.kp,occupied.kt,], counts.for.trans=neINtFORkp[,occupied.kp,], counts.for.pats=neINpFORkt[,occupied.kt,], counts.for.pheno=nvsINKp[occupied.kp,],
                    Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights)
      preds=SUBSTRA.predict(BO.model,as.matrix(expresData[,p]))
      AOerr=(labels[p]-1-preds[2])^2
      
      reverseExp=matrix(rep(2-expresData[,p],length(occupied.kp)),nrow = length(occupied.kp),ncol = Nt,byrow = T)
      logTheta=log(abs(reverseExp-Theta[,match(kt,occupied.kt),2]))
      
      for(curnu in c(magnitude)){
        # print(sprintf("rep = %d",rep))
        temp.weights=trans.weights
        # temp.weights=rep(magnitude,Nt)/ntINkt[kt]
        # temp.weights=rep(1,Nt)
        temp.weights.matrix=matrix(rep(temp.weights,length(occupied.kp)),nrow = length(occupied.kp),ncol = Nt,byrow = T)
        # temp.weights=rep(1,Nt)
        # featno=2000
        
        # inds=c(1:featno)
        # for(lii in 1:ceiling(Nt/featno)){
        #   ti=inds+featno*(lii-1)
        #   if(ti[featno]>Nt)
        #     ti=c(ti[1]:Nt)
        remaining=c(1:Nt)
        while(length(remaining) > 0){
          # ti=sample(remaining,min(featno,length(remaining)))
          ti=remaining
          remaining=setdiff(remaining,ti)
          # reverseExp=matrix(rep(2-expresData[ti,p],Kp),nrow = Kp,ncol = length(ti),byrow = T)
          # trans.weights.matrix=matrix(rep(trans.weights[ti],Kp),nrow = Kp,ncol = length(ti),byrow = T)
          
          logs=rowSums(logTheta[,ti]*temp.weights.matrix[,ti])+log(Pi)
          logs=logs+log(.Machine$double.xmax)-max(logs)-log(length(logs))
          values=exp(logs)
          values=values/sum(values)
          # print(paste0("values: ", values))
          A=sum(values*Psi[,2])
          # print(paste0("A: ",A))
          B=1
          # dA=colSums(matrix(rep(values*Psi[,2],Nt),byrow = F,nrow = Kp,ncol = Nt)*logTheta)
          dA=colSums(matrix(rep(values*Psi[,2],length(ti)),byrow = F,nrow = length(occupied.kp),ncol = length(ti))*logTheta[,ti])
          # print(paste0("dA: ",min(dA)," - ",max(dA)))
          # dB=colSums(matrix(rep(values,Nt),byrow = F,nrow = Kp,ncol = Nt)*logTheta)
          dB=colSums(matrix(rep(values,length(ti)),byrow = F,nrow = length(occupied.kp),ncol = length(ti))*logTheta[,ti])
          # print(paste0("dB: ",min(dB)," - ",max(dB)))
          f=A/B
          df=(dA*B-dB*A)/(B*B)
          # print(paste0("df: ",min(df)," - ",max(df)))
          # trans.weights=trans.weights+nu*2*(labels[p]-1-f)*df
          temp.weights[ti]=temp.weights[ti]+curnu*2*(labels[p]-1-f)*df
          # print("trans.weights:")
          # print(paste0(min(trans.weights)," - ",max(trans.weights)))
          
        }#lii
        temp.weights[which(temp.weights<0)]=0
        
        AO.model=list(patient.clusters=kp, transcript.clusters=match(kt,occupied.kt),
                      bicluster.content=nesINkpkt[occupied.kp,occupied.kt,], counts.for.trans=neINtFORkp[,occupied.kp,], counts.for.pats=neINpFORkt[,occupied.kt,], counts.for.pheno=nvsINKp[occupied.kp,],
                      Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=temp.weights)
        preds=SUBSTRA.predict(AO.model,as.matrix(expresData[,p]))
        curerr=(labels[p]-1-preds[2])^2
        
        # preds=SUBSTRA.predict(AO.model,as.matrix(expresData))
        # curerr=mean((labels-1-preds[,2])^2)
        # print(sprintf("%d -> %f", p, curerr))
        
        # if(is.nan(curerr)){
        #   next
        #   # print("not changed")
        # }else 
        if(curerr < AOerr){
          best.weights=temp.weights
          AOerr=curerr
          flag=T
        }
      }#curnu
      
      # print(sprintf("%d -> err: %f to %f - w in [%f, %f] - distance: %g",p,BOerr,AOerr,min(best.weights),max(best.weights),mean(abs(trans.weights-best.weights))))
      if(flag)
        counter=counter+1
      
      trans.weights=best.weights
      
      ### compute probability
      sig.kp=c(occupied.kp,empty.kp)
      sig.kt=c(occupied.kt,empty.kt)
      reverseExp=matrix(rep(2-expresData[,p],length(sig.kp)),nrow = length(sig.kp),ncol = Nt,byrow = T)
      Theta=(nesINkpkt[sig.kp,sig.kt,]+alpha_e/2)/array(rep(nesINkpkt[sig.kp,sig.kt,1]+nesINkpkt[sig.kp,sig.kt,2]+alpha_e,2),dim = c(length(sig.kp),length(sig.kt),2))
      p.prob=rowSums(log(abs(reverseExp-Theta[,match(kt,sig.kt),2]))*
                       matrix(rep(trans.weights,length(sig.kp)),nrow = length(sig.kp),ncol = Nt,byrow = T))+
        log((nvsINKp[sig.kp,labels[p]]+alpha_ph/Nv)/(npINkp[sig.kp]+alpha_ph))
      
      if(is.null(empty.kp))
        p.prob=p.prob+log(npINkp[occupied.kp])
      else
        p.prob=p.prob+log(c(npINkp[occupied.kp],alpha_p))
      
      # p.prob=log(npINkp + alpha_p / Kp)+
      #   rowSums(log(Theta[,,1])*matrix(rep(neINpFORkt[p,,1],Kp),nrow = Kp,ncol = Kt,byrow = T))+
      #   rowSums(log(Theta[,,2])*matrix(rep(neINpFORkt[p,,2],Kp),nrow = Kp,ncol = Kt,byrow = T))
      
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
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] + (2-expresData[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] + (expresData[,p]-1)
      
      if (withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]+1
      
      if(!is.null(empty.kp) && tempkp==empty.kp){
        occupied.kp=unique(kp)
        empty.kp=c()
        if(length(occupied.kp) < Kp)
          empty.kp=min(which(npINkp==0))
      }
    }#p
    
    if(length(which(nesINkpkt<0))>0)
      print("patients")
    
    print(paste0("Phase II - Iteration: ",ite))
    print("Difference:")
    print(paste0("p: ",p.diff," - ",length(occupied.kp)," clusters occupied."))
    print(paste0("t: ",t.diff," - ",length(occupied.kt)," clusters occupied."))
    print(paste0("w: ",counter))
    rkp=cluster_similarity(kp,prev.kp)
    rkt=cluster_similarity(kt,prev.kt)
    print("Rand Indexes:")
    print(paste0("Patient Clustering:    ",rkp))
    print(paste0("Transcript Clustering: ",rkt))
    
    print(paste0(min(trans.weights), " - ", max(trans.weights), " : ", mean(trans.weights)))
    Pi=npINkp[occupied.kp]/Np
    Theta=(nesINkpkt[occupied.kp,occupied.kt,]+alpha_e/2)/array(rep(nesINkpkt[occupied.kp,occupied.kt,1]+nesINkpkt[occupied.kp,occupied.kt,2]+alpha_e,2),dim = c(length(occupied.kp),length(occupied.kt),2))
    Psi=(nvsINKp[occupied.kp,]+alpha_ph/Nv)/array(rep(npINkp[occupied.kp]+alpha_ph,2),dim=c(length(occupied.kp),2))
    # Psi=nvsINKp/array(rep(rowSums(nvsINKp),2),dim=dim(nvsINKp))
    # Psi[is.nan(Psi)]=0.5
    model=list(patient.clusters=kp, transcript.clusters=match(kt,occupied.kt),
               bicluster.content=nesINkpkt[occupied.kp,occupied.kt,], counts.for.trans=neINtFORkp[,occupied.kp,], counts.for.pats=neINpFORkt[,occupied.kt,], counts.for.pheno=nvsINKp[occupied.kp,],
               Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights)
    
    preds=SUBSTRA.predict(model,expresData)
    r=roc(as.factor(labels),preds[,2],direction = "<")
    print(r$auc)
    err=mean((labels-1-preds[,2])^2)
    print(err)
    
    if(r$auc > bestauc || ite == 1){
      best.model=model
      besterr=err
      bestauc=r$auc
    }else if(r$auc == bestauc){
      if(err <= besterr){
        best.model=model
        besterr=err
        bestauc=r$auc
      }
    }
  }#phase_2_ite
  
  if(return.best)
    return(best.model)
  else
    return(model)
}



SUBSTRA.tune_parameters <- function(data, labels, Kp, Kt, nfolds=3, phase_1_ite, phase_2_ite,
                                   alpha_p = 1, alpha_t = 1, alpha_ph = 1e-100, alpha_e = 1, magnitudes=c(0.01,0.1,0.5,1,2), return.best=T){
  
  dataparts=SUBSTRA.split_dataset(data,labels,nfolds)
  bestAUC=0
  bestMagnitude=1
  for(magnitude in magnitudes){
    AUC=0
    for(i in 1:nfolds){
      trainData=cbind.data.frame(dataparts[[1]][-i])
      trainLabs=unlist(dataparts[[2]][-i])
      testData=as.data.frame(dataparts[[1]][i])
      testLabs=unlist(dataparts[[2]][i])
      
      model=SUBSTRA.train_np(as.matrix(trainData),trainLabs,phase_1_ite = phase_1_ite,phase_2_ite = phase_2_ite,
                            magnitude = magnitude,alpha_p = alpha_p,alpha_t = alpha_t,alpha_ph = alpha_ph,alpha_e = alpha_e,return.best = return.best)
      ps=SUBSTRA.predict(model,testData)[,2]
      r=roc(as.factor(testLabs), ps, direction = "<")
      AUC=AUC+length(testLabs)*r$auc
    }
    AUC=AUC/length(labels)
    
    if(AUC > bestAUC){
      bestAUC=AUC
      bestMagnitude=magnitude
    }
  }
  
  bestMagnitude
}



SUBSTRA.train_new <- function(expresData, labels=NULL, sideData=NULL,Kp, Kt, iterations, burn_in,
                             alpha_p = 1, alpha_t = 1, alpha_ph = 1e-100, alpha_e = 1, alpha_sd = 1,
                             magp=1, stdv=0.1, stdv0=2, return.best=T){
  # alpha_ph=alpha_ph^magp
  if(min(min(expresData))==0) expresData = expresData + 1
  withLabel=!is.null(labels)
  if(withLabel && min(labels)==0) labels = labels + 1
  withSideData=!is.null(sideData)
  if(withSideData && min(min(sideData))==0) sideData = sideData + 1
  Nv = length(unique(labels))
  Np = ncol(expresData)
  Nt = nrow(expresData)
  logKp=log(Kp)
  logKt=log(Kt)
  logDBLMax=log(.Machine$double.xmax)
  # Clustering
  kp = c()
  kt = c()
  npINkp = c()
  ntINkt = c()
  # Bicluster Counters & Parameters
  nesINkpkt = array(0,dim = c(Kp, Kt, 2))
  neINpFORkt = array(0,dim = c(Np, Kt, 2))
  neINtFORkp = array(0,dim = c(Nt, Kp, 2))
  trans.weights=array(magp,dim = c(1, Nt))
  # Side expresData Counters
  nvsINktkt = array(0,dim = c(Kt, Kt, 2))
  nvINtFORkt = array(0,dim = c(Nt, Kt, 2))
  # Labels Counters
  nvsINKp =  array(0,dim = c(Kp, Nv))
  
  # Random Initialization
  for(i in 1:Np)
    kp[i] = sample(seq(labels[i],Kp,by=2),1)
  kt = sample(c(1:Kt),Nt,replace = T)
  for(i in 1:Kp){
    npINkp[i]=length(which(kp==i))
    if(npINkp[i]!=0){
      if(withLabel)
        for(v in 1:Nv) nvsINKp[i,v]=length(which(labels[which(kp==i)]==v))
        for(j in 1:Kt){
          if(length(which(kt==j))==0) next
          ktinds=which(kt==j)
          kpinds=which(kp==i)
          nesINkpkt[i,j,1]=abs(sum((expresData[ktinds,kpinds]-2)*
                                     matrix(rep(trans.weights[ktinds],length(kpinds)),nrow = length(ktinds),ncol = length(kpinds),byrow = F)))
          nesINkpkt[i,j,2]=sum((expresData[ktinds,kpinds]-1)*
                                 matrix(rep(trans.weights[ktinds],length(kpinds)),nrow = length(ktinds),ncol = length(kpinds),byrow = F))
        }
        for(t in 1:Nt){
          neINtFORkp[t,i,1]= length(which(expresData[t,which(kp==i)]==1))*trans.weights[t]
          neINtFORkp[t,i,2]= length(which(expresData[t,which(kp==i)]==2))*trans.weights[t]
        }
    }
  }
  for(j in 1:Kt){
    ntINkt[j]=length(which(kt==j))
    if(ntINkt[j]!=0){
      if(withSideData){
        for(jj in 1:Kt){
          nvsINktkt[j,jj,1]=length(which(sideData[which(kt==j),which(kt==jj)]==1))
          nvsINktkt[j,jj,2]=length(which(sideData[which(kt==j),which(kt==jj)]==2))
        }
        for(t in 1:Nt){
          nvINtFORkt[t,j,1]=length(which(sideData[t,which(kt==j)]==1))
          nvINtFORkt[t,j,2]=length(which(sideData[t,which(kt==j)]==2))
        }
      }
      for(p in 1:Np){
        neINpFORkt[p,j,1]=abs(sum((expresData[kt==j,p]-2)*trans.weights[kt==j]))
        neINpFORkt[p,j,2]=sum((expresData[kt==j,p]-1)*trans.weights[kt==j])
      }
    }
  }
  
  
  Pi=(npINkp+alpha_p/Kp)/(Np+alpha_p)
  Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
  Psi=(nvsINKp+alpha_ph/Nv)/array(rep(rowSums(nvsINKp)+alpha_ph,2),dim=dim(nvsINKp))
  best.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
                  bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
                  Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights,
                  sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
  preds=SUBSTRA.predict(best.model,expresData)
  # print(preds)
  r=roc(as.factor(labels),preds[,2],direction = "<")
  bestauc=r$auc
  besterr=mean((labels-1-preds[,2])^2)
  print("Initial Status:")
  print(bestauc)
  print(besterr)
  
  # Sampling
  for(ite in 1:iterations){
    model.prob = 0
    
    ## Transcript Clusters
    t.diff = 0
    for(t in 1:Nt){
      ### take out
      tempkt=kt[t]
      ntINkt[tempkt]=ntINkt[tempkt]-1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] - neINtFORkp[t,,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1]-(2-expresData[t,])*trans.weights[t]
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2]-(expresData[t,]-1)*trans.weights[t]
      # for(p in 1:Np){
      #   neINpFORkt[p,tempkt,expresData[t,p]]=neINpFORkt[p,tempkt,expresData[t,p]]-trans.weights[t]
      # }
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]-nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]-nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]-ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]-ifelse(sideData[tt,t]==2,1,0)
        }
      }
      ### compute probability
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      t.prob=log(ntINkt + alpha_t / Kt)+
        colSums(log(Theta[,,1])*matrix(rep(neINtFORkp[t,,1]/magp),nrow = Kp,ncol = Kt,byrow = F))+
        colSums(log(Theta[,,2])*matrix(rep(neINtFORkp[t,,2]/magp),nrow = Kp,ncol = Kt,byrow = F))
      # colSums(log(Theta[,,1])*matrix(rep(neINtFORkp[t,,1]/ifelse(trans.weights[t]==0,1,trans.weights[t])*
      #                                    magt,Kt),nrow = Kp,ncol = Kt,byrow = F))+1
      # colSums(log(Theta[,,2])*matrix(rep(neINtFORkp[t,,2]/ifelse(trans.weights[t]==0,1,trans.weights[t])*
      #                                      magt,Kt),nrow = Kp,ncol = Kt,byrow = F))
      
      # aa=ntINkt/(2*stdv*stdv)+1/(2*stdv0*stdv0)
      # aap=(ntINkt+1)/(2*stdv*stdv)+1/(2*stdv0*stdv0)
      # tempdf=cbind.data.frame(cl=as.factor(kt),w=as.numeric(trans.weights))
      # ag=aggregate(tempdf$w,by=list(tempdf$cl),FUN=sum)
      # ag[ag[,1]==kt[t],2]=ag[ag[,1]==kt[t],2]-trans.weights[t]
      # trans.weights.sum=rep(0,Kt)
      # trans.weights.sum[ag[,1]]=ag[,2]
      # bb=-trans.weights.sum/(2*stdv*stdv)-magp/(2*stdv0*stdv0)
      # bbp=-(trans.weights.sum+trans.weights[t])/(2*stdv*stdv)-magp/(2*stdv0*stdv0)
      # 
      # t.prob=t.prob+(bbp*bbp/aap-bb*bb/aa)
      
      # t.prob=c()
      # for(cl in 1:Kt){
      #   aa=ntINkt[cl]/(2*stdv*stdv)+1/(2*stdv0*stdv0)
      #   aap=(ntINkt[cl]+1)/(2*stdv*stdv)+1/(2*stdv0*stdv0)
      #   trans.weights.sum=sum(trans.weights[kt==cl])
      #   if(cl==tempkt)
      #     trans.weights.sum=trans.weights.sum - trans.weights[t]
      #   bb=-trans.weights.sum/(2*stdv*stdv)-magp/(2*stdv0*stdv0)
      #   bbp=-(trans.weights.sum+trans.weights[t])/(2*stdv*stdv)-magp/(2*stdv0*stdv0)
      #   # t.prob[cl]=log(ntINkt[cl] + alpha_t / Kt)
      #   # t.prob[cl]=t.prob[cl]+sum(lgamma(nesINkpkt[,cl,1]+neINtFORkp[t,,1]+alpha_e/2)-lgamma(nesINkpkt[,cl,1]+alpha_e/2)+
      #   #                             lgamma(nesINkpkt[,cl,2]+neINtFORkp[t,,2]+alpha_e/2)-lgamma(nesINkpkt[,cl,2]+alpha_e/2)-
      #   #                             lgamma(nesINkpkt[,cl,1]+neINtFORkp[t,,1]+nesINkpkt[,cl,2]+neINtFORkp[t,,2]+alpha_e)+
      #   #                             lgamma(nesINkpkt[,cl,1]+nesINkpkt[,cl,2]+alpha_e))+
      #   #   (bbp*bbp/aap-bb*bb/aa)
      #   t.prob[cl]=log(ntINkt[cl] + alpha_t / Kt)
      #   # t.prob[cl]=t.prob[cl]+trans.weights[t]*(sum(log((nesINkpkt[,cl,1]+alpha_e/2)/(ntINkt[cl]*npINkp+alpha_e))*neINtFORkp[t,,1])+
      #   #   sum(log((nesINkpkt[,cl,2]+alpha_e/2)/(ntINkt[cl]*npINkp+alpha_e))*neINtFORkp[t,,2]))+
      #   #   (bbp*bbp/aap-bb*bb/aa)
      #   t.prob[cl]=t.prob[cl]+sum(log((nesINkpkt[,cl,1]+alpha_e/2)/(trans.weights.sum*npINkp+alpha_e))*neINtFORkp[t,,1])+
      #                            sum(log((nesINkpkt[,cl,2]+alpha_e/2)/(trans.weights.sum*npINkp+alpha_e))*neINtFORkp[t,,2])
      #     # (bbp*bbp/aap-bb*bb/aa)
      #   if(withSideData){
      #     t.prob[cl]=t.prob[cl]+sum(log((nesINktkt[,cl,1]+alpha_sd/2)/(ntINkt[cl]*ntINkt+alpha_sd))*nvINtFORkt[t,,1])+
      #       sum(log((nesINktkt[,cl,2]+alpha_sd/2)/(ntINkt[cl]*ntINkt+alpha_sd))*nvINtFORkt[t,,2])
      #   }
      # }
      
      # t.prob=t.prob*trans.weights[t]
      
      ### put back
      t.prob.max=max(t.prob)
      t.prob.max = t.prob.max - (logDBLMax - logKt)
      t.prob = exp(t.prob - t.prob.max)
      total.prob = sum(t.prob)
      tempkt=sample(1:Kt,1,prob = t.prob/total.prob)
      if (kt[t] != tempkt)
        t.diff=t.diff+1
      model.prob = model.prob + log(t.prob[tempkt] / total.prob)
      
      kt[t]=tempkt
      ntINkt[tempkt]=ntINkt[tempkt]+1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] + neINtFORkp[t,,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1]+(2-expresData[t,])*trans.weights[t]
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2]+(expresData[t,]-1)*trans.weights[t]
      # for(p in 1:Np){
      #   neINpFORkt[p,tempkt,expresData[t,p]]=neINpFORkt[p,tempkt,expresData[t,p]]+trans.weights[t]
      # }
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]+nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]+nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]+ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]+ifelse(sideData[tt,t]==2,1,0)
        }
      }
    }#t
    
    print("Transcripts Done!")
    # print(Np*sum(trans.weights))
    # print(sum(nesINkpkt))
    # print(sum(neINpFORkt))
    # print(sum(neINtFORkp))
    print(length(which(ntINkt>0)))
    
    # trans.weights=rep(magp,Nt)
    ## Patient Clusters
    p.diff = 0
    counter=0
    for(p in 1:Np){
      ### take out
      tempkp=kp[p]
      prevkp=kp[p]
      kp[p]=-1
      npINkp[tempkp]=npINkp[tempkp]-1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] - neINpFORkt[p,,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] - (2-expresData[,p])*trans.weights
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] - (expresData[,p]-1)*trans.weights
      # for(t in 1:Nt)
      #   neINtFORkp[t,tempkp,expresData[t,p]]=neINtFORkp[t,tempkp,expresData[t,p]]-trans.weights[t]
      if (withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]-1
      
      
      if(ite > burn_in){
        ### compute weights
        flag=F
        best.weights=trans.weights
        Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
        Pi=(npINkp+alpha_p/Kp)/(Np-1+alpha_p)
        Psi=(nvsINKp+alpha_ph/Nv)/(npINkp+alpha_ph)
        
        BO.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
                      bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
                      Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights,
                      sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
        preds=SUBSTRA.predict(BO.model,as.matrix(expresData[,p]))
        AOerr=(labels[p]-1-preds[2])^2
        
        reverseExp=matrix(rep(2-expresData[,p],Kp),nrow = Kp,ncol = Nt,byrow = T)
        logTheta=log(abs(reverseExp-Theta[,kt,2]))
        
        for(curnu in c(magp)){
          temp.weights=trans.weights
          temp.weights.matrix=matrix(rep(temp.weights,Kp),nrow = Kp,ncol = Nt,byrow = T)
          remaining=c(1:Nt)
          while(length(remaining) > 0){
            ti=remaining
            remaining=setdiff(remaining,ti)
            
            logs=rowSums(logTheta[,ti]*temp.weights.matrix[,ti])+log(Pi)
            logs=logs+log(.Machine$double.xmax)-max(logs)-log(length(logs))
            values=exp(logs)
            values=values/sum(values)
            # print(paste0("values: ", values))
            # print(max(values))
            A=sum(values*Psi[,2])
            # print(paste0("A: ",A))
            B=1
            # dA=colSums(matrix(rep(values*Psi[,2],Nt),byrow = F,nrow = Kp,ncol = Nt)*logTheta)
            dA=colSums(matrix(rep(values*Psi[,2],length(ti)),byrow = F,nrow = Kp,ncol = length(ti))*logTheta[,ti])
            # print(paste0("dA: ",min(dA)," - ",max(dA)))
            # dB=colSums(matrix(rep(values,Nt),byrow = F,nrow = Kp,ncol = Nt)*logTheta)
            dB=colSums(matrix(rep(values,length(ti)),byrow = F,nrow = Kp,ncol = length(ti))*logTheta[,ti])
            # print(paste0("dB: ",min(dB)," - ",max(dB)))
            f=A/B
            df=(dA*B-dB*A)/(B*B)
            # print(paste0("df: ",min(df)," - ",max(df)))
            # trans.weights=trans.weights+nu*2*(labels[p]-1-f)*df
            temp.weights[ti]=temp.weights[ti]+curnu*2*(labels[p]-1-f)*df
            # print("trans.weights:")
            # print(paste0(min(trans.weights)," - ",max(trans.weights)))
          }#lii
          
          # epsilon=1e-300*magp
          temp.weights[which(temp.weights<0)]=0
          
          AO.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
                        bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
                        Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=temp.weights,
                        sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
          preds=SUBSTRA.predict(AO.model,as.matrix(expresData[,p]))
          curerr=(labels[p]-1-preds[2])^2
          
          if(curerr < AOerr){
            best.weights=temp.weights
            AOerr=curerr
            flag=T
          }
        }#curnu
        
        # print(sprintf("%d -> err: %f to %f - w in [%f, %f] - distance: %g",p,BOerr,AOerr,min(best.weights),max(best.weights),mean(abs(trans.weights-best.weights))))
        if(flag)
          counter=counter+1
        
        ### Updating the Parameters
        # print("before:")
        # print(nesINkpkt)
        nesINkpkt[,,]=0
        neINpFORkt[,,]=0
        neINtFORkp[,,]=0
        
        for(pcl in 1:Kp){
          if(npINkp[pcl]==0)
            next
          kpinds=which(kp==pcl)
          multiplier=matrix(rep(best.weights,npINkp[pcl]),nrow = Nt,ncol = npINkp[pcl],byrow = F)
          neINtFORkp[,pcl,1]=rowSums((2-expresData[,kpinds])*multiplier)
          neINtFORkp[,pcl,2]=rowSums((expresData[,kpinds]-1)*multiplier)
        }
        for(tcl in 1:Kt){
          if(ntINkt[tcl]==0)
            next
          ktinds=which(kt==tcl)
          multiplier=matrix(rep(best.weights[ktinds],Np),nrow = ntINkt[tcl],ncol = Np,byrow = F)
          neINpFORkt[,tcl,1]=colSums((2-expresData[ktinds,])*multiplier)
          neINpFORkt[,tcl,2]=colSums((expresData[ktinds,]-1)*multiplier)
          
          if(ntINkt[tcl]==1){
            nesINkpkt[,tcl,1]=neINtFORkp[ktinds,,1]
            nesINkpkt[,tcl,2]=neINtFORkp[ktinds,,2]
          }else{
            nesINkpkt[,tcl,1]=colSums(as.matrix(neINtFORkp[ktinds,,1]))
            nesINkpkt[,tcl,2]=colSums(as.matrix(neINtFORkp[ktinds,,2]))
          }
        }
        
        # print(Np*sum(trans.weights))
        # print(sum(nesINkpkt))
        # print(sum(neINpFORkt))
        # print(sum(neINtFORkp))
        
        w.diff=Nt-length(which(trans.weights == best.weights))
        
        trans.weights=best.weights
      } # if(ite)
      
      ### compute probability
      reverseExp=matrix(rep(2-expresData[,p],Kp),nrow = Kp,ncol = Nt,byrow = T)
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      # p.prob=rowSums(log(abs(reverseExp-Theta[,kt,2]))*
      #                  matrix(rep(trans.weights,Kp),nrow = Kp,ncol = Nt,byrow = T))+
      #   log(npINkp + alpha_p / Kp)+
      #   log((nvsINKp[,labels[p]]+alpha_ph/Nv)/(npINkp+alpha_ph))
      
      p.prob=log(npINkp + alpha_p / Kp)+
        log((nvsINKp[,labels[p]]+alpha_ph/Nv)/(npINkp+alpha_ph))+
        rowSums(log(Theta[,,1])*matrix(rep(neINpFORkt[p,,1],Kp),nrow = Kp,ncol = Kt,byrow = T))+
        rowSums(log(Theta[,,2])*matrix(rep(neINpFORkt[p,,2],Kp),nrow = Kp,ncol = Kt,byrow = T))
      
      # print(rowSums(log(abs(reverseExp-Theta[,kt,2]))*
      # matrix(rep(trans.weights,Kp),nrow = Kp,ncol = Nt,byrow = T)))
      # print(log(npINkp + alpha_p / Kp))
      # print(log((nvsINKp[,labels[p]]+alpha_ph/Nv)/(npINkp+alpha_ph)))
      
      # p.prob.max=max(p.prob[seq(labels[p],Kp,by=2)])
      p.prob.max=max(p.prob)
      # print(p.prob)
      ### put back
      p.prob.max = p.prob.max - (logDBLMax - logKp)
      p.prob = exp(p.prob - p.prob.max)
      # print(sprintf("label: %d",labels[p]))
      # p.prob[seq(3-labels[p],Kp,by=2)]=0
      # print(p.prob/sum(p.prob))
      total.prob = sum(p.prob)
      # print(round(p.prob/total.prob,2))
      tempkp=sample(1:Kp,1,prob = p.prob/total.prob)
      if (prevkp != tempkp)
        p.diff=p.diff+1
      model.prob = model.prob + log(p.prob[tempkp] / total.prob)
      
      kp[p]=tempkp
      npINkp[tempkp]=npINkp[tempkp]+1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] + neINpFORkt[p,,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] + (2-expresData[,p])*trans.weights
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] + (expresData[,p]-1)*trans.weights
      # for(t in 1:Nt)
      #   neINtFORkp[t,tempkp,expresData[t,p]]=neINtFORkp[t,tempkp,expresData[t,p]]+trans.weights[t]
      if (withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]+1
    }#p
    
    print("Patients Done!")
    # print(Np*sum(trans.weights))
    # print(sum(nesINkpkt))
    # print(sum(neINpFORkt))
    # print(sum(neINtFORkp))
    print(length(which(npINkp>0)))
    
    print(paste0("Iteration: ",ite))
    print("Difference:")
    print(p.diff)
    print(t.diff)
    print(counter)
    
    print(paste0(min(trans.weights), " - ", max(trans.weights), " : ", mean(trans.weights)))
    Pi=(npINkp+alpha_p/Kp)/(Np+alpha_p)
    Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
    Psi=(nvsINKp+alpha_ph/Nv)/array(rep(rowSums(nvsINKp)+alpha_ph,2),dim=dim(nvsINKp))
    
    model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
               bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
               Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights,
               sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
    preds=SUBSTRA.predict(model,expresData)
    r=roc(as.factor(labels),preds[,2],direction = "<")
    print(r$auc)
    err=mean((labels-1-preds[,2])^2)
    print(err)
    
    if(r$auc > bestauc || ite == burn_in){
      best.model=model
      besterr=err
      bestauc=r$auc
    }else if(r$auc == bestauc){
      if(err <= besterr){
        best.model=model
        besterr=err
        bestauc=r$auc
      }
    }
  }#ite
  
  if(return.best)
    return(best.model)
  else
    return(model)
}



SUBSTRA.update <- function(oldModel, expresData, labels=NULL, sideData=NULL,Kp, Kt, iterations, burn_in=3,
                          alpha_p = 1, alpha_t = 1, alpha_ph = 1e-4, alpha_e = 1, alpha_sd = 1,mag=1){
  if(min(min(expresData))==0) expresData = expresData + 1
  withLabel=!is.null(labels)
  if(withLabel && min(labels)==0) labels = labels + 1
  withSideData=!is.null(sideData)
  if(withSideData && min(min(sideData))==0) sideData = sideData + 1
  Nv = length(unique(labels))
  Np = ncol(expresData)
  Nt = nrow(expresData)
  logKp=log(Kp)
  logKt=log(Kt)
  logDBLMax=log(.Machine$double.xmax)
  # Clustering
  kp = c()
  kt = oldModel$transcript.clusters
  npINkp = oldModel$pat.cluster.counts
  ntINkt = oldModel$trans.cluster.counts
  # Bicluster Counters & Parameters
  nesINkpkt = oldModel$bicluster.content
  neINpFORkt = array(0,dim = c(Np, Kt, 2))
  neINtFORkp = oldModel$counts.for.trans
  trans.weights=oldModel$trans.weights
  # Side expresData Counters
  nvsINktkt = oldModel$sd.bicluster.counts
  nvINtFORkt = oldModel$sd.counts.for.trans
  # Labels Counters
  nvsINKp =  oldModel$counts.for.pheno
  
  # Random Initialization
  for(i in 1:Np)
    kp[i] = sample(seq(labels[i],Kp,by=2),1)
  # kt = sample(c(1:Kt),Nt,replace = T)
  for(i in 1:Kp){
    npINkp[i]=npINkp[i]+length(which(kp==i))
    if(npINkp[i]!=0){
      if(withLabel)
        for(v in 1:Nv) nvsINKp[i,v]=nvsINKp[i,v]+length(which(labels[which(kp==i)]==v))
        for(j in 1:Kt){
          if(ntINkt[j]==0) next
          nesINkpkt[i,j,1]=nesINkpkt[i,j,1]+length(which(expresData[which(kt==j),which(kp==i)]==1))
          nesINkpkt[i,j,2]=nesINkpkt[i,j,2]+length(which(expresData[which(kt==j),which(kp==i)]==2))
        }
        for(t in 1:Nt){
          neINtFORkp[t,i,1]=neINtFORkp[t,i,1]+length(which(expresData[t,which(kp==i)]==1))
          neINtFORkp[t,i,2]=neINtFORkp[t,i,2]+length(which(expresData[t,which(kp==i)]==2))
        }
    }
  }
  for(j in 1:Kt){
    if(ntINkt[j]!=0){
      for(p in 1:Np){
        neINpFORkt[p,j,1]=length(which(expresData[which(kt==j),p]==1))
        neINpFORkt[p,j,2]=length(which(expresData[which(kt==j),p]==2))
      }
    }
  }
  
  # print(min(nesINkpkt))
  # print(min(neINtFORkp))
  # print(min(neINpFORkt))
  
  Pi=(npINkp+alpha_p/Kp)/(sum(npINkp)+alpha_p)
  Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
  # Psi=(nvsINKp+alpha_ph/Nv)/array(rep(rowSums(nvsINKp)+alpha_ph,2),dim=dim(nvsINKp))
  Psi=nvsINKp/array(rep(rowSums(nvsINKp),2),dim=dim(nvsINKp))
  Psi[is.nan(Psi)]=0.5
  best.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
                  bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
                  Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights,
                  sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
  preds=SUBSTRA.predict(best.model,expresData,mag)
  # print(preds)
  r=roc(as.factor(labels),preds[,2],direction = "<")
  bestauc=r$auc
  besterr=mean((labels-1-preds[,2])^2)
  print("Initial Status:")
  print(bestauc)
  print(besterr)
  
  # Sampling
  for(ite in 1:iterations){
    model.prob = 0
    # Pi=(npINkp+alpha_p/Kp)/(Np+alpha_p)
    # Theta=(nesINkpkt+alpha_e)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e[,,1]+alpha_e[,,2],2),dim = dim(nesINkpkt))
    # Psi=(nvsINKp+alpha_ph/Nv)/array(rep(rowSums(nvsINKp)+alpha_ph,2),dim=dim(nvsINKp))
    # 
    # print(paste0("Iteration: ",ite-1))
    # print(paste0(min(trans.weights)," - ",max(trans.weights)))
    # BS.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
    #               bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
    #               Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights,
    #               sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
    # preds=SUBSTRA.predictWeighted2(BS.model,expresData)
    # # print(preds)
    # r=roc(as.factor(labels),preds[,2],direction = "<")
    # BSauc=r$auc
    # print(r$auc)
    # print(mean((labels-1-preds[,2])^2))
    
    ## Transcript Clusters
    t.diff = 0
    for(t in 1:Nt){
      ### take out
      tempkt=kt[t]
      ntINkt[tempkt]=ntINkt[tempkt]-1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] - neINtFORkp[t,,]
      for(p in 1:Np){
        neINpFORkt[p,tempkt,expresData[t,p]]=neINpFORkt[p,tempkt,expresData[t,p]]-1
      }
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]-nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]-nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]-ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]-ifelse(sideData[tt,t]==2,1,0)
        }
      }
      ### compute probability
      t.prob.max=-1*.Machine$double.xmax
      t.prob=c()
      for(cl in 1:Kt){
        t.prob[cl]=log(ntINkt[cl] + alpha_t / Kt)
        # print(nesINkpkt[,cl,2])
        # print(ntINkt[cl]*npINkp)
        # print((nesINkpkt[,cl,2]+alpha_e/2)/(ntINkt[cl]*npINkp+alpha_e))
        t.prob[cl]=t.prob[cl]+sum(log((nesINkpkt[,cl,1]+alpha_e/2)/(ntINkt[cl]*npINkp+alpha_e))*neINtFORkp[t,,1])+
          sum(log((nesINkpkt[,cl,2]+alpha_e/2)/(ntINkt[cl]*npINkp+alpha_e))*neINtFORkp[t,,2])
        if(withSideData){
          t.prob[cl]=t.prob[cl]+sum(log((nesINktkt[,cl,1]+alpha_sd/2)/(ntINkt[cl]*ntINkt+alpha_sd))*nvINtFORkt[t,,1])+
            sum(log((nesINktkt[,cl,2]+alpha_sd/2)/(ntINkt[cl]*ntINkt+alpha_sd))*nvINtFORkt[t,,2])
        }
        
        if(t.prob[cl]>t.prob.max)
          t.prob.max=t.prob[cl]
      }
      
      ### put back
      t.prob.max = t.prob.max - (logDBLMax - logKt)
      t.prob = exp(t.prob - t.prob.max)
      total.prob = sum(t.prob)
      tempkt=sample(1:Kt,1,prob = t.prob/total.prob)
      if (kt[t] != tempkt)
        t.diff=t.diff+1
      model.prob = model.prob + log(t.prob[tempkt] / total.prob)
      
      kt[t]=tempkt
      ntINkt[tempkt]=ntINkt[tempkt]+1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] + neINtFORkp[t,,]
      for(p in 1:Np){
        neINpFORkt[p,tempkt,expresData[t,p]]=neINpFORkt[p,tempkt,expresData[t,p]]+1
      }
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]+nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]+nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]+ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]+ifelse(sideData[tt,t]==2,1,0)
        }
      }
    }#t
    
    
    # trans.weights=rep(0,Nt)
    ## Patient Clusters
    p.diff = 0
    counter=0
    for(p in 1:Np){
      ### take out
      tempkp=kp[p]
      npINkp[tempkp]=npINkp[tempkp]-1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] - neINpFORkt[p,,]
      for(t in 1:Nt)
        neINtFORkp[t,tempkp,expresData[t,p]]=neINtFORkp[t,tempkp,expresData[t,p]]-1
      if (withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]-1
      
      ### compute weights
      reverseExp=matrix(rep(2-expresData[,p],Kp),nrow = Kp,ncol = Nt,byrow = T)
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      Pi=(npINkp+alpha_p/Kp)/(sum(npINkp)+alpha_p)
      # Psi=(nvsINKp+alpha_ph/Nv)/array(rep(rowSums(nvsINKp)+alpha_ph,2),dim=dim(nvsINKp))
      Psi=nvsINKp/array(rep(rowSums(nvsINKp),2),dim=dim(nvsINKp))
      Psi[is.nan(Psi)]=0.5
      
      best.weights=trans.weights
      
      # print(sprintf("Before Optimization: %d",p))
      BO.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
                    bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
                    Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=best.weights,
                    sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
      preds=SUBSTRA.predict(BO.model,as.matrix(expresData[,p]),mag)
      # r=roc(as.factor(labels),preds[,2],direction = "<")
      BOerr=(labels[p]-1-preds[2])^2
      # print(BOerr)
      
      curmodel=NULL
      curauc=0
      AOerr=BOerr
      flag=F
      
      logTheta=log(abs(reverseExp-Theta[,kt,2]))
      if(ite > burn_in){
        curnu=mag
        temp.weights=trans.weights
        temp.weights.matrix=matrix(rep(temp.weights,Kp),nrow = Kp,ncol = Nt,byrow = T)
        ti=c(1:Nt)
        logs=rowSums(logTheta[,ti]*temp.weights.matrix[,ti])+mag*log(Pi)
        logs=logs+log(.Machine$double.xmax)-max(logs)-log(length(logs))
        values=exp(logs)
        values=values/sum(values)
        # print(paste0("values: ", values))
        A=sum(values*Psi[,2])
        # print(paste0("A: ",A))
        B=1
        dA=colSums(matrix(rep(values*Psi[,2],length(ti)),byrow = F,nrow = Kp,ncol = length(ti))*logTheta[,ti])
        # print(paste0("dA: ",min(dA)," - ",max(dA)))
        dB=colSums(matrix(rep(values,length(ti)),byrow = F,nrow = Kp,ncol = length(ti))*logTheta[,ti])
        # print(paste0("dB: ",min(dB)," - ",max(dB)))
        f=A/B
        df=(dA*B-dB*A)/(B*B)
        # print(paste0("df: ",min(df)," - ",max(df)))
        temp.weights[ti]=temp.weights[ti]+curnu*2*(labels[p]-1-f)*df
        # print("trans.weights:")
        # print(paste0(min(trans.weights)," - ",max(trans.weights)))
        
        AO.model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
                      bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
                      Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=temp.weights,
                      sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
        preds=SUBSTRA.predict(AO.model,as.matrix(expresData[,p]),mag)
        curerr=(labels[p]-1-preds[2])^2
        # print(sprintf("%d -> %f",p,curerr))
        
        if(is.nan(curerr)){
          next
          # print("not changed")
        }else if(curerr < AOerr){
          best.weights=temp.weights
          AOerr=curerr
          flag=T
        }
      }#if
      
      # print(sprintf("%d -> err: %f to %f - w in [%f, %f]",p,BOerr,AOerr,min(best.weights),max(best.weights)))
      if(flag)
        counter=counter+1
      
      trans.weights=best.weights
      
      ### compute probability
      reverseExp=matrix(rep(2-expresData[,p],Kp),nrow = Kp,ncol = Nt,byrow = T)
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      p.prob=rowSums(log(abs(reverseExp-Theta[,kt,2]))*
                       matrix(rep(best.weights,Kp),nrow = Kp,ncol = Nt,byrow = T))+
        mag*log(npINkp + alpha_p / Kp)
      p.prob.max=max(p.prob[seq(labels[p],Kp,by=2)])
      
      ### put back
      p.prob.max = p.prob.max - (logDBLMax - logKp)
      p.prob = exp(p.prob - p.prob.max)
      # print(sprintf("label: %d",labels[p]))
      p.prob[seq(3-labels[p],Kp,by=2)]=0
      # print(p.prob/sum(p.prob))
      total.prob = sum(p.prob)
      tempkp=sample(1:Kp,1,prob = p.prob/total.prob)
      if (kp[p] != tempkp)
        p.diff=p.diff+1
      model.prob = model.prob + log(p.prob[tempkp] / total.prob)
      
      kp[p]=tempkp
      npINkp[tempkp]=npINkp[tempkp]+1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] + neINpFORkt[p,,]
      for(t in 1:Nt)
        neINtFORkp[t,tempkp,expresData[t,p]]=neINtFORkp[t,tempkp,expresData[t,p]]+1
      if (withLabel)
        nvsINKp[tempkp,labels[p]]=nvsINKp[tempkp,labels[p]]+1
    }#p
    
    print(counter)
    print(paste0("Iteration: ",ite))
    print("Difference:")
    print(p.diff)
    print(t.diff)
    
    print(paste0(min(trans.weights), " - ", max(trans.weights), " : ", mean(trans.weights)))
    Pi=(npINkp+alpha_p/Kp)/(sum(npINkp)+alpha_p)
    Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
    Psi=(nvsINKp+alpha_ph/Nv)/array(rep(rowSums(nvsINKp)+alpha_ph,2),dim=dim(nvsINKp))
    # Psi=nvsINKp/array(rep(rowSums(nvsINKp),2),dim=dim(nvsINKp))
    # Psi[is.nan(Psi)]=0.5
    model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
               bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt, counts.for.pheno=nvsINKp,
               Pi=Pi,Theta=Theta, Psi=Psi,trans.weights=trans.weights,
               sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
    preds=SUBSTRA.predict(model,expresData,mag)
    r=roc(as.factor(labels),preds[,2],direction = "<")
    print(r$auc)
    err=mean((labels-1-preds[,2])^2)
    print(err)
    
    if(r$auc > bestauc){
      best.model=model
      besterr=err
      bestauc=r$auc
    }else if(r$auc == bestauc){
      if(err <= besterr){
        best.model=model
        besterr=err
        bestauc=r$auc
      }
    }
  }#ite
  
  return(best.model)
}



SUBSTRA.predict <- function(model, expresData){
  if(min(min(expresData))==0) expresData = expresData + 1
  probs=matrix(0,nrow = ncol(expresData),ncol = ncol(model$Psi))
  for(i in 1:ncol(expresData)){
    reverseExp=matrix(rep(2-expresData[,i],length(model$Pi)),nrow = length(model$Pi),ncol = nrow(expresData),byrow = T)
    logs=rowSums(log(abs(reverseExp-model$Theta[,model$transcript.clusters,2]))*
                   matrix(rep(model$trans.weights,length(model$Pi)),nrow=length(model$Pi),ncol=nrow(expresData),byrow=T))+
      log(model$Pi)
    logs=logs+log(.Machine$double.xmax)-max(logs)-log(length(logs))
    values=exp(logs)
    clusterProbs=values/sum(values)
    probs[i,] = colSums(matrix(rep(clusterProbs,2),nrow=nrow(model$Psi),ncol=ncol(model$Psi))*model$Psi)
  }
  return(round(probs,6))
}



B2PS <- function(data, sideData=NULL, Kp, Kt, iterations, alpha_p = 1, alpha_t = 1, alpha_e = 1, alpha_sd = 1){
  if(min(min(data))==0) data = data + 1
  withSideData=!is.null(sideData)
  if(withSideData && min(min(sideData))==0) sideData = sideData + 1
  Np = ncol(data)
  Nt = nrow(data)
  logKp=log(Kp)
  logKt=log(Kt)
  logDBLMax=log(.Machine$double.xmax)
  # Clustering
  kp = c()
  kt = c()
  npINkp = c()
  ntINkt = c()
  # Bicluster Counters & Parameters
  nesINkpkt = array(0,dim = c(Kp, Kt, 2))
  neINpFORkt = array(0,dim = c(Np, Kt, 2))
  neINtFORkp = array(0,dim = c(Nt, Kp, 2))
  # Side Data Counters
  nvsINktkt = array(0,dim = c(Kt, Kt, 2))
  nvINtFORkt = array(0,dim = c(Nt, Kt, 2))
  
  # Random Initialization
  kp = sample(c(1:Kp),Np,replace = T)
  kt = sample(c(1:Kt),Nt,replace = T)
  for(i in 1:Kp){
    npINkp[i]=length(which(kp==i))
    if(npINkp[i]!=0){
      for(j in 1:Kt){
        if(length(which(kt==j))==0) next
        nesINkpkt[i,j,1]=length(which(data[which(kt==j),which(kp==i)]==1))
        nesINkpkt[i,j,2]=length(which(data[which(kt==j),which(kp==i)]==2))
      }
      for(t in 1:Nt){
        neINtFORkp[t,i,1]= length(which(data[t,which(kp==i)]==1))
        neINtFORkp[t,i,2]= length(which(data[t,which(kp==i)]==2))
      }
    }
  }
  for(j in 1:Kt){
    ntINkt[j]=length(which(kt==j))
    if(ntINkt[j]!=0){
      if(withSideData){
        for(jj in 1:Kt){
          nvsINktkt[j,jj,1]=length(which(sideData[which(kt==j),which(kt==jj)]==1))
          nvsINktkt[j,jj,2]=length(which(sideData[which(kt==j),which(kt==jj)]==2))
        }
        for(t in 1:Nt){
          nvINtFORkt[t,j,1]=length(which(sideData[t,which(kt==j)]==1))
          nvINtFORkt[t,j,2]=length(which(sideData[t,which(kt==j)]==2))
        }
      }
      for(p in 1:Np){
        neINpFORkt[p,j,1]=length(which(data[which(kt==j),p]==1))
        neINpFORkt[p,j,2]=length(which(data[which(kt==j),p]==2))
      }
    }
  }
  
  
  # Sampling
  for(ite in 1:iterations){
    print(paste0("Iteration: ",ite))
    model.prob = 0
    ## Transcript Clusters
    t.diff = 0
    for(t in 1:Nt){
      ### take out
      tempkt=kt[t]
      ntINkt[tempkt]=ntINkt[tempkt]-1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] - neINtFORkp[t,,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] - (2-data[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] - (data[t,]-1)
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]-nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]-nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]-ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]-ifelse(sideData[tt,t]==2,1,0)
        }
      }
      
      ### compute probability
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      t.prob=log(ntINkt + alpha_t / Kt)+
        colSums(log(Theta[,,1])*matrix(rep(neINtFORkp[t,,1],Kt),nrow = Kp,ncol = Kt,byrow = F))+
        colSums(log(Theta[,,2])*matrix(rep(neINtFORkp[t,,2],Kt),nrow = Kp,ncol = Kt,byrow = F))
      
      ### put back
      t.prob.max=max(t.prob)
      t.prob.max = t.prob.max - (logDBLMax - logKt)
      t.prob = exp(t.prob - t.prob.max)
      total.prob = sum(t.prob)
      tempkt=sample(1:Kt,1,prob = t.prob/total.prob)
      if (kt[t] != tempkt)
        t.diff=t.diff+1
      model.prob = model.prob + log(t.prob[tempkt] / total.prob)
      
      kt[t]=tempkt
      ntINkt[tempkt]=ntINkt[tempkt]+1
      nesINkpkt[,tempkt,] = nesINkpkt[,tempkt,] + neINtFORkp[t,,]
      neINpFORkt[,tempkt,1]=neINpFORkt[,tempkt,1] + (2-data[t,])
      neINpFORkt[,tempkt,2]=neINpFORkt[,tempkt,2] + (data[t,]-1)
      if (withSideData){
        nvsINktkt[tempkt,,] = nvsINktkt[tempkt,,]+nvINtFORkt[tempkt,]
        nvsINktkt[,tempkt,] = nvsINktkt[,tempkt,]+nvINtFORkt[tempkt,]
        nvsINktkt[tempkt,tempkt,2]=nvsINktkt[tempkt,tempkt,2]+1
        for(tt in 1:Nt)
        {
          nvINtFORkt[tt,tempkt,1] = nvINtFORkt[tt,tempkt,1]+ifelse(sideData[tt,t]==1,1,0)
          nvINtFORkt[tt,tempkt,2] = nvINtFORkt[tt,tempkt,2]+ifelse(sideData[tt,t]==2,1,0)
        }
      }
    }#t
    
    ## Patient Clusters
    p.diff = 0
    for(p in 1:Np){
      ### take out
      tempkp=kp[p]
      npINkp[tempkp]=npINkp[tempkp]-1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] - neINpFORkt[p,,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] - (2-data[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] - (data[,p]-1)
      
      ### compute probability
      reverseExp=matrix(rep(2-data[,p],Kp),nrow = Kp,ncol = Nt,byrow = T)
      Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
      p.prob=rowSums(log(abs(reverseExp-Theta[,kt,2])))+
        log(npINkp + alpha_p / Kp)
      p.prob.max=max(p.prob)
      p.prob.max = p.prob.max - (logDBLMax - logKp)
      p.prob = exp(p.prob - p.prob.max)
      total.prob = sum(p.prob)
      tempkp=sample(1:Kp,1,prob = p.prob/total.prob)
      if (kp[p] != tempkp)
        p.diff=p.diff+1
      model.prob = model.prob + log(p.prob[tempkp] / total.prob)
      
      ### put back
      kp[p]=tempkp
      npINkp[tempkp]=npINkp[tempkp]+1
      nesINkpkt[tempkp,,] = nesINkpkt[tempkp,,] + neINpFORkt[p,,]
      neINtFORkp[,tempkp,1] = neINtFORkp[,tempkp,1] + (2-data[,p])
      neINtFORkp[,tempkp,2] = neINtFORkp[,tempkp,2] + (data[,p]-1)
    }#p
    
    print("Difference:")
    print(paste0("p: ",p.diff))
    print(paste0("t: ",t.diff))
    
  }#ite
  
  Pi=(npINkp+alpha_p/Kp)/(Np+alpha_p)
  Theta=(nesINkpkt+alpha_e/2)/array(rep(nesINkpkt[,,1]+nesINkpkt[,,2]+alpha_e,2),dim = dim(nesINkpkt))
  model=list(patient.clusters=kp, transcript.clusters=kt, pat.cluster.counts=npINkp, trans.cluster.counts=ntINkt,
             bicluster.content=nesINkpkt, counts.for.trans=neINtFORkp, counts.for.pats=neINpFORkt,
             Pi=Pi,Theta=Theta,
             sd.bicluster.counts=nvsINktkt,sd.counts.for.trans=nvINtFORkt)
  return(model)
}
