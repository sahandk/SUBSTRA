pathToSUBSTRA="./SUBSTRA.R"
source(pathToSUBSTRA)

data=read.table("./Inputs/data.txt")
labels=read.table("./Inputs/labels.txt")[,1]


# (1) Training an Individual Model
model=SUBSTRA.train(data = data,labels = labels,magnitude = 1,phase_1_ite = 50,phase_2_ite = 20)
print(model)


# (2) Tuning 'magnitude' and Training an Individual Model
magnitude=SUBSTRA.tune(data = data,labels = labels,pathToSUBSTRA = pathToSUBSTRA,magnitudes = c(0.1,1,10),nfolds = 3,
                       phase_1_ite = 45,phase_2_ite = 15,parallel = T,ncores = 3,verbose = T,trainVerbose = F)
model=SUBSTRA.train(data = data,labels = labels,magnitude = magnitude,phase_1_ite = 50,phase_2_ite = 20,verbose = F)
print(model)


# (3) Concensus Method: Tunes the 'magnitude' and trains several models and aggregates them.
consensus_model=SUBSTRA.train_consensus(data = data,labels = labels,ntimes = 5,parallel = T,ncores = 5,verbose = T,pathToSUBSTRA = pathToSUBSTRA,
                                        magnitudes = c(0.1,1,10),nfolds_tune = 3,phase_1_ite_tune = 45,phase_2_ite_tune = 15,parallel_tune = T,
                                        phase_1_ite_train = 50,phase_2_ite_train = 20)
print(consensus_model)


# (4) Drawing the heatmaps of the outputs of (2) and (3)
SUBSTRA.heatmap(data = data,labels = labels,model = model,pathToPDF = "./model_heatmap.pdf","Model Heatmap")
SUBSTRA.heatmap(data = data,labels = labels,model = consensus_model,pathToPDF = "./consensus_model_heatmap.pdf","Consensus Model Heatmap")
