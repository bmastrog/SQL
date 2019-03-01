## Load example data, from `Hiiragi2013`.
if (FALSE) { # install packages if not available
  source("https://bioconductor.org/biocLite.R")
  biocLite("Hiiragi2013")
  install.packages("mlr")
}
library("knitr")
library("Biobase")
library("Hiiragi2013")
library("glmnet")
library("mlr")
library("ggplot2")
library("limma")
library("kernlab")

#Question 1

#load mouse embryo expression profile data
data( "x", package = "Hiiragi2013" )
table( x$sampleGroup )

#convert expression set x to data frame
rowV <- data.frame( v = rowVars(exprs(x)) )

#non-specific/ class-independent filtering
selectionThreshold <- 10^(0.5)
selectedFeatures <- ( rowV$v > selectionThreshold )
print(sum(selectedFeatures))

embryoSingleCells <- data.frame( t(exprs(x)[selectedFeatures, ]), check.names = TRUE )

#2-class classification: E3.25 v. everything else
#response variable
embryoSingleCells$tg <- factor( ifelse( x$Embryonic.day == "E3.25", "E3.25", "other") )

# randomise 
embryoSingleCells$tg = sample(embryoSingleCells$tg) 

# show counts in each class
with( embryoSingleCells, table( tg ) ) 

#define the task
task <- makeClassifTask( id = "Hiiragi", data = embryoSingleCells, target = "tg" )

#filter genes by significance of fold-change
fv = generateFilterValuesData(task, method = "anova.test")
plotFilterValues(fv) # look at the distribution of the gene p-values

#select just the top 30 genes based on significance of fold change
filtered.task = filterFeatures(task, method = "anova.test", abs=30)

#simple form of neural net classifier
lrn = makeLearner("classif.nnet")

#10-fold cross-validation (repeated 2 times) 
# provides prediction error rate.
rdesc <- makeResampleDesc( method = "RepCV", stratify = TRUE, folds = 10, reps = 2)
r <- resample(learner = lrn, task = filtered.task, resampling = rdesc )
r$aggr # mean error rate

#The mean error rate is 0.1728283

#At first, this result was surprising because the expected error rate for a 
#random sample is 0.5.  After reading the paper by MacLachlan, it became evident 
#that this mean error rate had selection bias.  
#The paper attributes feature selection as a source of selection bias in 
#predictive models of expression data.  Before the models can be generated, 
#genes need to be filtered out to remove genes that have “little or no 
#discriminatory power” (6562).  If the cross-validated error is calculated within 
#this filtered subset of genes, the error has selection bias.  This occurs because 
#the feature selection has already looked at the labels of the training data and 
#has used them, thus the data has already been partially trained.  This training 
#needs to be included in validation.  The expected error rate in this scenario is 0, which 
#the 0.17 error value is close to.  In the above error rate, the feature selection 
#for the top 30 genes with the most significant fold change values was used for the 
#cross-validation (filtered.task was used in CV), so it has selection bias. 
#MacLachlan provides a series of examples that demonstrate that the CV error rate 
#is reduced with selection bias.   To obtain the corrected error rate, we would have 
#to include the feature selection step in the cross-validation.   By ensuring cross 
#validation occurs externally to feature selection, we can correct for selection 
#bias.

#combine filter and classifier
lrn1 = makeFilterWrapper(learner = "classif.nnet", fw.method = "anova.test", fw.abs= 30)
rdesc1 = makeResampleDesc("RepCV", stratify = TRUE, folds = 10, reps = 2)
r1 = resample(learner = lrn1, task = task, resampling = rdesc1, show.info = FALSE, models = TRUE)
r1$aggr

#mmce.test.mean 0.4565152
#This error rate is much closer to the expected .5 value for random data.  
#This value is expected because the calculations corrected for sample bias.  
#In the code describing the learner, the makeFilterWrapper function combines both 
#the learner, neural net classifier, and the feature selection parameters, 30 
#genes with most significant fold change.  This ensures that the CV is done 
#externally to feature selection, thus corrects for selection bias.   

#unshuffle labels
embryoSingleCells$tg <- factor( ifelse( x$Embryonic.day == "E3.25", "E3.25", "other") )

#define the task
task2 <- makeClassifTask( id = "Hiiragi", data = embryoSingleCells, target = "tg" )

#filter genes by significance of fold-change
fv2 = generateFilterValuesData(task, method = "anova.test")

#select just the top 30 genes based on significance of fold change
filtered.task2 = filterFeatures(task2, method = "anova.test", abs=30)

#simple form of neural net classifier
lrn2 = makeLearner("classif.nnet")

#10-fold cross-validation (repeated 2 times) 
# provides prediction error rate
rdesc2 <- makeResampleDesc( method = "RepCV", stratify = TRUE, folds = 10, reps = 2)
r2 <- resample(learner = lrn2, task = filtered.task2, resampling = rdesc2 )
r2$aggr # mean error rate

#mmce.test.mean 0.1683838


#repeat with combined filter and classifier
lrn3 = makeFilterWrapper(learner = "classif.nnet", fw.method = "anova.test", fw.abs= 30)
rdesc3 = makeResampleDesc("RepCV", stratify = TRUE, folds = 10, reps = 2)
r3 = resample(learner = lrn3, task = task2, resampling = rdesc3, show.info = FALSE, models = TRUE)
r3$aggr

#mmce.test.mean 0.160404

#Without the label shuffling, the data is no longer random, so the training is 
#done on an actual subset of the data the model aims to predict.  We would expect 
#the error rate in both of these scenarios to be much lower than the expected 0.5 
#for random data, since the data is no longer randomized.  However, the calculations between 
#the two values are different.  The first performing cross-validation on a 
#feature-selected data, so it has selection bias.  The second calculation has 
#cross-validation done externally to feature selection, so does not have selection
#bias.  The second calculation reflects the correct mean-error rate.

#Question 2

# shuffled labels

# randomise 
embryoSingleCells$tg = sample(embryoSingleCells$tg) 

#execute nested CV

#Define the learner:
lrn7 = "classif.ksvm" 

#set parameters
p7 = makeParamSet(makeNumericParam("C", lower = -1, upper = 1, trafo = function(x) 2^x),
                  makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x))

#generate wrapped learner for nested resampling
tuningLrn7 <- makeTuneWrapper(lrn7, 
                              resampling = makeResampleDesc("CV", iters = 3,  stratify = TRUE), 
                              par.set = p7, 
                              control = makeTuneControlGrid(resolution = 5) )
#define the task
task7 <- makeClassifTask( id = "Hiiragi", data = embryoSingleCells, target = "tg" )

#define resampling strategy
rdesc7 <- makeResampleDesc("RepCV", fold = 10, rep= 3 )

#execute nested CV
r7 <- resample(learner = tuningLrn7, task = task7, resampling = rdesc7 )
r7$aggr #mean error rate

#The mean error rate is 0.495.  This value is much closer to the expected value than 
#the value the neural net model provided, which was 0.45.

#unshuffled labels

#re-order labels
embryoSingleCells$tg <- factor( ifelse( x$Embryonic.day == "E3.25", "E3.25", "other") )

#define the learner
lrn6 = "classif.ksvm" 

#set the parameters
p6 = makeParamSet(makeNumericParam("C", lower = -1, upper = 1, trafo = function(x) 2^x),
                 makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x))
#generate wrapped learner for nested CV
tuningLrn6 <- makeTuneWrapper(lrn6, 
                             resampling = makeResampleDesc("CV", iters = 3,  stratify = TRUE), 
                             par.set = p6, 
                             control = makeTuneControlGrid(resolution = 5) )
#define the task
task6 <- makeClassifTask( id = "Hiiragi", data = embryoSingleCells, target = "tg" )

#define the resampling method
rdesc6 <- makeResampleDesc("RepCV", fold = 10, rep= 3 )

#execute nested CV
r6 <- resample(learner = tuningLrn6, task = task6, resampling = rdesc6 )
r6$aggr #mean error rate 

#The mean error rate is 0.097.  This is much lower than the value obtained for the
#neural net model which was 0.16.
