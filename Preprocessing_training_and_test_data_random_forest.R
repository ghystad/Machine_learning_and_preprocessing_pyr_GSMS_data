## Written by Grethe Hystad, ghystad@pnw.edu
## Reviewed by Anirudh Prabhu, aprabhu@carnegiescience.edu

## Data consist of 134 samples of organic molecular mixtures analyzed by pyrolysis Gass chromatogray-Mass Spectrometry (pyr-GC-MS).
## Each data set has the scan values in the rows and the Mass to Charge ratio (m/z) values in the columns.
## The data values are the raw intensities as measured by pyr-GC-MS.
## Each sample has either an abiotic (A) or biotic (B) origin. 
## Purpose is to build classification rules to predict whether the samples are biotic or abiotic.
## Purpose is also to determine the features, (m/z) and (scan) values, that discriminate between the biotic and the abiotic groups.
## Random Forest is used as the machine learning method.
## The data are preprocessed before we apply random forest.

library(dplyr)         # for dataframe computation
library(MALDIquant)    # for chemometrics processing
library(caret)         # for machine learning
library(mlr3)          # for machine learning
library("mlr3verse")   # for machine learning
library("mlr3learners")# for machine learning
library("mlr3tuning")  # for machine learning
library("data.table")  # for rbindlist
library("ggplot2")     # for plots


###Data preparation
setwd("C:/Users/ghystad/Desktop/Data2022")
species=c(rep("A",5),"B", "A","A","B","B","A","C", rep("B",5),rep("C",4),rep("B",4),"A",rep("B",11),"A","A","B","A","A","B","B","A",rep("B",3),"A","B",rep("A",6),rep("B",3),"A","A","B","B",rep("A",4),"B",rep("A",5),"B",rep("A",10),rep("B",4),rep("A",5),rep("B",3),rep("A",3),rep("B",4),rep("A",8),"B","A","A","B","C","C","A",rep("B",3),"A","B", "B",    "A","B",rep("A",3),"B","A", rep("A",7),"B","A","B")
species2=as.numeric(as.factor(species))
ind=which(species2=="3")
species_n=species2[-ind]

#Reading in the data is based on the code found in: "How to import multiple .csv files simultaneously in R and create a data frame" by Dahna, A., datascience+, 03/09/2019
#https://datascienceplus.com/how-to-import-multiple-csv-files-simultaneously-in-r-and-create-a-data-frame/
fnames <- list.files()
fnames2=fnames[-ind]
z1=lapply(fnames2, read.delim)  #read in all of 134 datasets
z=z1
NN=700 #number of m/z values
mass=seq(50,NN,1) #m/z 
MM=3240 #number of scans
Ndim=length(fnames2)


M=list()
for(i in 1:Ndim){
	colnames(z[[i]])="mass"
	#remove commas
	z[[i]]=data.frame(do.call("rbind", strsplit(as.character(z[[i]]$mass), ",", fixed = TRUE)))
	z[[i]]=data.frame(lapply(z[[i]],as.numeric))
	colnames(z[[i]])=c("scan",as.character(seq(50,NN,1)))
	z[[i]]=z[[i]] %>% slice(1:MM)       #Selects the first MM rows  
	M[[i]]=z[[i]]
	}

################################################################

y=species_n
y=factor(y,labels=c("A","B"))

set.seed(99)

#Stratified random sampling
inTraining = createDataPartition(y, p = 0.7, list = FALSE)

#Training and testing indices
train_set  = as.vector(inTraining)
test_set=  setdiff(seq(1,Ndim),train_set)
y_train=y[train_set]
y_test=y[test_set]

prop.table(table(y_train))
prop.table(table(y_test))

N=length(train_set)
Nt=length(test_set)

################################################################
#Training preprocessing
#Detect the significant peaks as local max above four times the signal to noise ratio

#number of m/z values to use
MZ=151
#Create Chromatograms for each training sample and each m/z value inside each training sample
Sample.list=list()
for (i in 1:N){
     S=list()
     	for(j in 1:MZ){
	     S[[j]] = createMassSpectrum(mass=seq(1,MM,1), intensity=unlist(M[[train_set[i]]][,j+1]),
	  metaData=list(name="Chrom")) 
	 }
   chrom = transformIntensity(S,method="sqrt") #stabilize the variance
   chrom = smoothIntensity(chrom, method="MovingAverage",halfWindowSize=5)#Smooth the data
   chrom = removeBaseline(chrom, method="SNIP") #Remove the baseline

   # Put the processed chromatograms back into a dataframe
   Processed.chrom.list=list()
	for (k in 1:MZ){
	     Processed.chrom.list[[k]] = as.numeric(intensity(chrom[[k]]))
	}

   Processed.Mass_dataframe = as.data.frame(do.call(rbind, Processed.chrom.list))
   Ma=max(Processed.Mass_dataframe)
   Mi=min(Processed.Mass_dataframe)
   #Normalize across sample
   Processed.Mass_dataframe=t((Processed.Mass_dataframe-Mi)/(Ma-Mi))
   #Processed.Mass_dataframe=as.data.frame(Processed.Mass_dataframe)
   S2=list()
	for(t in 1:MZ){
	    S2[[t]] = createMassSpectrum(mass=seq(1,MM,1), intensity=Processed.Mass_dataframe[,t],
	    metaData=list(name="Chrom_normalized"))  
	 }

    peaks = detectPeaks(S2, method="MAD", halfWindowSize=20, SNR=4)
    peak.list=list()
	for(tt in 1:MZ){
      	v=numeric(MM)
        	scan.number=mass(peaks[[tt]])
       	v[scan.number] = intensity(peaks[[tt]])
	  	peak.list[[tt]] = v
	}
 	Processed.peaks = t(as.data.frame(do.call(rbind, peak.list)))
      row.names(Processed.peaks)=c(paste0("R", 1:MM))
      colnames(Processed.peaks)=c(paste0("M", 50:(MZ+50-1)))
	Sample.list[[i]] = Processed.peaks
}

mass.scan.list=list()
for(i in 1:N){
    sm=as.numeric(unlist(Sample.list[[i]]))     #Scan and mass spectrum for each training sample
    mass.scan.list[[i]]=sm
}

#Sample vs mass/scan numbers
data.mass.scan = do.call(rbind, mass.scan.list) #Put the training samples into a dataframe

Bin=as.character(seq(1,3240,1))
MS=as.character(seq(50,(MZ+50-1),1))

colnames(data.mass.scan)=paste(outer(Bin, MS, paste, sep = ';'))


#Next we use hierarchical clustering with a distance of 20 to group the scan numbers for each m/z value.
#The method is based on the paper: "Sample classification from proteing mass spectrometry by peak probability contrasts"
#by Tibshirani, R., et al. (2004), Bioinformatis, 20(17):3034-44, doi:10.1093/bioinformatics/bth357 
ml=50:(MZ+50-1)
mll=length(50:(MZ+50-1))

data.mass.scan_new=as.data.frame(ifelse(data.mass.scan>0,1,0))
MZ_name=c(paste0(";",50:(MZ+50-1)))
MZ_name2=c(paste0(50:(MZ+50-1)))
MZ_name3=c(paste0(".",50:(MZ+50-1)))

L_original=list()   # Scan numbers for each cluster per m/z value selected 
LL_dataframe=list() #The max intensities for the ith m/z value for each bin of scan numbers and each training sample

for (i in 1:MZ){
      data.mass.scan_new2 = data.mass.scan_new %>% select(ends_with(MZ_name[i])) #Select columns that ends with m/z=i
      scan_name=sub(basename(colnames(data.mass.scan_new2)),pattern = MZ_name[i], replacement = "" , fixed = TRUE)
      colnames(data.mass.scan_new2) = scan_name  #Name the columns with the scan numbers
	f=function(x){
      y=sum(x)
      }
      N_elements=as.numeric(apply(as.matrix(data.mass.scan_new2),2,f)) #Count the number of elements for each scan number across the training samples
      vec=as.numeric(rep(colnames(data.mass.scan_new2), times=N_elements)) #Repeat the nonzero scan numbers with its frequency
      hc=hclust(dist(vec),method="complete")  # Hierarchical clustering on the scan numbers
      clusters=cutree(hc,h=20)                # Cluster with a distance of 20
      dataframeN=as.data.frame(data.mass.scan) %>% select(ends_with(c(paste0(MZ_name[i]))))
	 L.original2=list()  #Scan numbers for each cluster for the ith m/z value selected 
       L_dataframe=list()  #The max intensity for the jth cluster in the ith m/z value for each training sample
		 for (j in 1:length(unique(clusters))){
    			indd=which(clusters==unique(clusters)[j])
			#L[[j]]=vec[indd]
			L.original2[[j]]=as.data.frame(unique(vec[indd]))
			colnames(L.original2[[j]])=MZ_name2[i]
			Lt=vec[indd]
			L_mean=round(mean(Lt))   #Mean number of scan number in each cluster
			dataframeN2=dataframeN[,unique(Lt)]
				va=if(length(unique(Lt))>1){
				apply(dataframeN2,1,max) #Find the max intensity over the scan numbers in one cluster for each training sample
				}else {dataframeN2}
			DD=data.frame(matrix(ncol = 1, nrow =  N)) #The max intensity for the jth cluster in the ith m/z value for each training sample
			colnames(DD)= paste(outer(as.character(L_mean),MZ_name2[i], paste, sep = ';'))
			DD[,1]=va
			L_dataframe[[j]]=DD
		   }
	 L_original[[i]]=L.original2
	 LL_dataframe[[i]]=do.call(cbind, L_dataframe) #The max intensities for the ith m/z value for each bin of scan numbers and each training sample
}

data.mass.scan33=do.call(cbind, LL_dataframe) # Dataframe with the training samples as rows and the intertwined m/z values and scan numbers as columns

names2=colnames(data.mass.scan33)

lg=function(x){(length(which(x>0)))/(dim(data.mass.scan33)[2])}
Perc_Nnonzero=apply(data.mass.scan33,1,lg) # Ratio of nonzero feature values for each observation

data.mass.scan33 = data.mass.scan33 %>% mutate(Perc_Nnonzero)

colnames(data.mass.scan33) = make.names(colnames(data.mass.scan33))
names4=colnames(data.mass.scan33)

#Removing variables with near zero variance, nearZeroVar, and 
#correlated predictors, findCorrelation, are based on the book
#"Applied Predictive Modeling" by Kuhn, M. and Johnson, K., Springer, New York, 2013 
near.zero.variance = nearZeroVar(data.mass.scan33) #Variable with near zero variance

data.mass.scan33 = data.mass.scan33[, -near.zero.variance] #Remove features with near zero variance
dim(data.mass.scan33) #new dimension 

#Omit correlated predictors from the training set
corr.ma=cor(data.mass.scan33)
Corr.h = findCorrelation(corr.ma, 0.85)
data.mass.scan33 = data.mass.scan33[, -Corr.h]
dim(data.mass.scan33)  #new dimension 

#Contains training samples with the selected features and corresponding y-values (A or B)
training_transformed=data.frame(data.mass.scan33,y=y_train)

#################################################################################
#Test preprocessing

#number of m/z values to use
MZ=151
#Create Chromatograms for each sample and each m/z value inside each test sample
Sample.list2=list()
for (i in 1:Nt){
  	S=list()
		for(j in 1:MZ){
	  		S[[j]] = createMassSpectrum(mass=seq(1,MM,1), intensity=unlist(M[[test_set[i]]][,j+1]),
	  		metaData=list(name="Chrom"))  
	 	}
    chrom = transformIntensity(S,method="sqrt") #stabilize the variance
    chrom = smoothIntensity(chrom, method="MovingAverage",halfWindowSize=5)
    chrom = removeBaseline(chrom, method="SNIP")

    # Put the processed chromatograms back into a dataframe
    Processed.chrom.list=list()
		for (k in 1:MZ){
			Processed.chrom.list[[k]] = as.numeric(intensity(chrom[[k]]))
		}
    Processed.Mass_dataframe = as.data.frame(do.call(rbind, Processed.chrom.list))
    Ma=max(Processed.Mass_dataframe)
    Mi=min(Processed.Mass_dataframe)
    #Normalize across sample
    Processed.Mass_dataframe=t((Processed.Mass_dataframe-Mi)/(Ma-Mi))
    #Processed.Mass_dataframe=as.data.frame(Processed.Mass_dataframe)
    S2=list()
		for(t in 1:MZ){
	  		S2[[t]] = createMassSpectrum(mass=seq(1,MM,1), intensity=Processed.Mass_dataframe[,t],
	  		metaData=list(name="Chrom_normalized"))  
		 }
	peaks = detectPeaks(S2, method="MAD", halfWindowSize=20, SNR=4)
    	peak.list=list()
		for (tt in 1:MZ){
     			 v=numeric(MM)
      		 scan.number=mass(peaks[[tt]])
     			 v[scan.number] = intensity(peaks[[tt]])
			 peak.list[[tt]] = v
		}
 	Processed.peaks = t(as.data.frame(do.call(rbind, peak.list)))
      row.names(Processed.peaks)=c(paste0("R", 1:MM))
     	colnames(Processed.peaks)=c(paste0("M", 50:(MZ+50-1)))
    	Sample.list2[[i]] = Processed.peaks
}

mass.scan.list_test=list()
for(i in 1:Nt){
  	sm=as.numeric(unlist(Sample.list2[[i]]))    #Scan and mass spectrum for each test set
 	mass.scan.list_test[[i]]=sm
}

#Sample vs mass/scan numbers
data.mass.scan_test = do.call(rbind, mass.scan.list_test) #Put the test data into a dataframe
data.mass.scan_test=as.data.frame(data.mass.scan_test)
colnames(data.mass.scan_test)=paste(outer(Bin, MS, paste, sep = ';'))
colnames(data.mass.scan_test) = make.names(colnames(data.mass.scan_test))


LL_dataframeT=list()
for (i in 1:MZ){
	dataframeNT=as.data.frame(data.mass.scan_test) %>% select(ends_with(c(paste0(MZ_name3[i]))))
	#LT=list()
	#L.original2T=list()
	L_dataframeT=list()
		for (j in 1: length(L_original[[i]])){
			L.ind2=as.numeric(unlist(L_original[[i]][j]))
			L.ind=min(L.ind2):max(L.ind2)
			dataframeNT2=dataframeNT[,L.ind]
				vaT=if(length(L.ind)>1){
				apply(dataframeNT2,1,max) #Find the max intensity over the scan numbers in one cluster for each test sample
				}else {dataframeNT2}
			DDT=data.frame(matrix(ncol = 1, nrow =  Nt)) #The max intensity for the jth cluster in the ith m/z value for each test sample
			DDT[,1]=vaT
			L_dataframeT[[j]]=DDT
		}
	LL_dataframeT[[i]]=do.call(cbind, L_dataframeT) #The max intensities for the ith m/z value over each binned scan number for each test sample
}
data.mass.scan33T=do.call(cbind, LL_dataframeT) # Dataframe with the test samples as rows and the intertwined m/z values and scan numbers as columns
colnames(data.mass.scan33T)=names2


lgT=function(x){(length(which(x>0)))/(dim(data.mass.scan33T)[2])}
Perc_NnonzeroT=apply(data.mass.scan33T,1,lgT) # Ratio of nonzero feature values for each observation

data.mass.scan33T = data.mass.scan33T %>% mutate(Perc_NnonzeroT)
colnames(data.mass.scan33T)=names4
names3=colnames(data.mass.scan33)

data.mass.scan33T=data.mass.scan33T  %>% select(all_of(names3))

#Contains test samples
testing_transformed=data.frame(data.mass.scan33T,y=y_test)

#Contains both training and test samples
data_transformed=rbind(training_transformed,testing_transformed)

########################################################
#Machine learning
#The machine learning R-code is based on the mlr3 library and the book https://mlr3book.mlr-org.com/, #Flexible and Robust Machine Learning Using mlr3 in R#
#by Lang, M. et al

#Create a task 
task=as_task_classif(training_transformed, target = "y", positive="B")

#Using stratified random sampling for each resampling
task$col_roles$stratum = "y"

#To check the strata
#prop.table(table(task$data(rows = r$test_set(1), cols = "y")))

###############################
#Ranger
set.seed(99)
resampling_cv10 = rsmp("cv", folds = 10)  #For 10-fold cross validation
resampling_cv10$instantiate(task)

#Create the learner
learner = lrn("classif.ranger", id = "rf",predict_type="prob",importance = "impurity")

graph = po("removeconstants")  %>>% learner
plot(graph)

graph_learner = as_learner(graph)
graph_learner$param_set$ids()

#rf
graph_learner$param_set$values$rf.mtry = to_tune(p_int(1,2000))
graph_learner$param_set$values$rf.num.trees = to_tune(p_int(1,8000))

graph_learner$id = "graph_learner"


set.seed(99)
future::plan("multisession")

#using stratified random sampling inside each fold
instance = tune(
  method = "random_search",
  task =  task,
  learner = graph_learner,
  resampling = resampling_cv10,
  measure = msr("classif.auc"),
  term_evals = 20
)


instance$result_y
#as.data.table(instance$archive)
#instance$archive$benchmark_result
instance$result
##################################################################

#Ranger
#Use the results from above

learner2 = lrn("classif.ranger", id = "rf",predict_type="prob",importance = "impurity")

graph2 = po("removeconstants")  %>>% learner2

graph2$param_set$ids()

graph_learner2 = as_learner(graph2)

graph_learner2$id = "graph2"

graph_learner2$param_set$values$rf.mtry = 388
graph_learner2$param_set$values$rf.num.trees =  5096  

graph_learner2$train(task)  

#Prediction on the testing set
prediction_rf = graph_learner2$predict_newdata(newdata = testing_transformed)

prediction_rf$confusion
measure = msr("classif.auc")
measure$predict_type
prediction_rf$score(measure)
#prediction_rf$prob

#print(graph_learner2$model)

#test_set[which(prediction_rf$truth!=prediction_rf$response)]


#########################################################
#Determine importance variables
#Use the results from above

learner3 = lrn("classif.ranger", id = "rf",predict_type="prob",importance = "impurity")

learner3$param_set$values$mtry =      388     
learner3$param_set$values$num.trees  =    5096 

learner3$train(task)

#Importance filter
filter = flt("importance", learner = learner3)
filter$calculate(task)
dff=as.data.table(filter)
#write.csv(dff,"C:/Users/ghystad/Documents/Templeton_grant_research/First_paper_material/For_journal_submission/training.rf.importance.filter.csv")




