require(mlr)
library(parallelMap)
setwd('~/git/MIE_Metagenomics')

task = readRDS('data/processed/fcbf_COUNTRY')
task = makeClassifTask(data = task, target = 'target')
print('Removing Constant Features')
task = removeConstantFeatures(task)
print('Normalizing Features')
task = normalizeFeatures(task)
print('Done task!')

# Hyperparameter tuning
ctrl<-makeTuneControlGrid()  #Hacer que se pueda modificar la búsqueda!!!!
inner<-makeResampleDesc("Holdout")
# GLMNET
psglmnet = makeParamSet(
  makeDiscreteParam("lambda", c(0.0001,0.001,0.01,0.1,1)),
  makeDiscreteParam("alpha",c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1))
)
l<-makeLearner("classif.glmnet", predict.type = "prob")
lrn_glmnet<-makeTuneWrapper(l, inner, psglmnet, measures = auc, ctrl, show.info=T)
 
# Random Forest
psrf<-makeParamSet(
  makeDiscreteParam("mtry", values = sqrt(ncol(n[[i]]$env$data))),
  makeDiscreteParam("ntree", values= 1000L),
  makeDiscreteParam("nodesize", values= c(1:3))
)
l<-makeLearner("classif.randomForest", predict.type = "prob")
lrn_rf<-makeTuneWrapper(l,  resampling = inner, par.set = psrf, measures = auc, control=ctrl,  show.info = T)


learners = list(lrn_ksvm, lrn_glmnet, lrn_rf)

# Outer Cross-Validation
outer = makeResampleDesc('RepCV' , reps = 5, folds = 10 , stratify = T)

print('Training the model')

# Benchmarking
parallelStartMulticore(2L , level = 'mlr.tuneParams')
bmr = benchmark(learners, task, outer, measures =  list(acc , auc, mmce), show.info = T, models = T)
parallelStop()


