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
## We apply nested resampling to obtain the generalized prediction performance

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
N=length(fnames2)


M=list()
for(i in 1:N){
	colnames(z[[i]])="mass"
	#remove commas
	z[[i]]=data.frame(do.call("rbind", strsplit(as.character(z[[i]]$mass), ",", fixed = TRUE)))
	z[[i]]=data.frame(lapply(z[[i]],as.numeric))
	colnames(z[[i]])=c("scan",as.character(seq(50,NN,1)))
	z[[i]]=z[[i]] %>% slice(1:MM)       #Selects the first MM rows
  	M[[i]]=z[[i]]
}

#saveRDS(M,"C:/Users/ghystad/Documents/Templeton_grant_research/GC_MS_R_files/Templeton_R_files/M.RDS")
#M=readRDS("C:/Users/ghystad/Documents/Templeton_grant_research/GC_MS_R_files/Templeton_R_files/M.RDS")


################
y=species_n
y=factor(y,labels=c("A","B"))

######################################################################
#Preprocessing
#Detect the significant peaks as local max above four times the signal to noise ratio

#number of m/z values to use
MZ=151
#Create Chromatograms for each sample and each m/z value inside each sample
Sample.list=list()
for (i in 1:N){
  S=list()
	for(j in 1:MZ){
	  S[[j]] = createMassSpectrum(mass=seq(1,MM,1), intensity=unlist(M[[i]][,j+1]),
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
   Processed.Mass_dataframe=as.data.frame(Processed.Mass_dataframe)
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
	Sample.list[[i]] = Processed.peaks
}

mass.scan.list=list()
for(i in 1:N){
  sm=as.numeric(unlist(Sample.list[[i]]))     #scan and mass spectrum for each sample
  mass.scan.list[[i]]=sm
}

#Sample vs mass/scan numbers
data.mass.scan = do.call(rbind, mass.scan.list) #Put the samples into a dataframe

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

L_original=list()  # Scan numbers for each cluster per m/z value selected 
LL_dataframe=list() #The max intensities for the ith m/z value for each bin of scan numbers and each sample

for (i in 1:MZ){
	data.mass.scan_new2 = data.mass.scan_new %>% select(ends_with(MZ_name[i])) #Select columns that ends with m/z=i
	scan_name=sub(basename(colnames(data.mass.scan_new2)),pattern = MZ_name[i], replacement = "" , fixed = TRUE)
	colnames(data.mass.scan_new2) = scan_name  #Name the columns with the scan numbers
	f=function(x){
	y=sum(x)
	}
	N_elements=as.numeric(apply(as.matrix(data.mass.scan_new2),2,f)) #Count the number of  elements for each scan number across the samples
	vec=as.numeric(rep(colnames(data.mass.scan_new2), times=N_elements)) #Repeat the nonzero scan numbers with its frequency
	hc=hclust(dist(vec),method="complete")  # Hiearcial clustering on the scan numbers
	clusters=cutree(hc,h=20)                # Cluster with a distance of 20
	dataframeN=as.data.frame(data.mass.scan) %>% select(ends_with(c(paste0(MZ_name[i]))))
	L.original2=list() #Scan numbers for each cluster for the ith m/z value selected 
	L_dataframe=list() #The max intensity for the jth cluster in the ith m/z value for each sample
		for (j in 1:length(unique(clusters))){
			indd=which(clusters==unique(clusters)[j])
			L.original2[[j]]=as.data.frame(unique(vec[indd]))
			colnames(L.original2[[j]])=MZ_name2[i]
			Lt=vec[indd]
			L_mean=round(mean(Lt))   #Mean number of scan number in each cluster
			dataframeN2=dataframeN[,unique(Lt)]
				va=if(length(unique(Lt))>1){
				apply(dataframeN2,1,max) #Find the max intensity over the scan numbers in one cluster for each sample
				}else {dataframeN2}
			DD=data.frame(matrix(ncol = 1, nrow =  N)) #The max intensity for the jth cluster in the ith m/z value for each sample
			colnames(DD)= paste(outer(as.character(L_mean),MZ_name2[i], paste, sep = ';'))
			DD[,1]=va
			L_dataframe[[j]]=DD
		}
	L_original[[i]]=L.original2
	LL_dataframe[[i]]=do.call(cbind, L_dataframe) #The max intensities for the ith m/z value for each bin of scan numbers and each sample
}

data.mass.scan33=do.call(cbind, LL_dataframe) # Dataframe with the samples as rows and the intertwined m/z values and scan numbers as columns
names2=colnames(data.mass.scan33)


lg=function(x){(length(which(x>0)))/(dim(data.mass.scan33)[2])}
Perc_Nnonzero=apply(data.mass.scan33,1,lg) # Ratio of nonzero feature values for each observation

data.mass.scan33 = data.mass.scan33 %>% mutate(Perc_Nnonzero)

colnames(data.mass.scan33) = make.names(colnames(data.mass.scan33))

#Removing variables with near zero variance, nearZeroVar, and 
#correlated predictors, findCorrelation, are based on the book
#"Applied Predictive Modeling" by Kuhn, M. and Johnson, K., Springer, New York, 2013 

near.zero.variance = nearZeroVar(data.mass.scan33) #Remove features with near zero variance
data.mass.scan33 = data.mass.scan33[, -near.zero.variance] #remove variables with near zero variance
dim(data.mass.scan33) #new dimension 

#Omit correlated predictors from the samples
corr.ma=cor(data.mass.scan33)
Corr.h = findCorrelation(corr.ma, 0.85)
data.mass.scan33 <- data.mass.scan33[, -Corr.h]
dim(data.mass.scan33)  #new dimension 


#Contains all samples with the selected features and corresponding y-values (A or B)
training_transformed=data.frame(data.mass.scan33,y=y)


########################################################
#Machine learning
#The machine learning R-code is based on the mlr3 library and the book https://mlr3book.mlr-org.com/, #Flexible and Robust Machine Learning Using mlr3 in R#
#by Lang, M. et al

#Create a task 
task=as_task_classif(training_transformed, target = "y", positive="B")

#Using stratified random sampling for each resampling
task$col_roles$stratum = "y"

###############################
#Ranger

set.seed(99)
inner_loop = rsmp("cv", folds = 10) # inner loop for nested CV
outer_loop = rsmp("cv", folds = 5)
outer_loop$instantiate(task)

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

at = AutoTuner$new(
  learner = graph_learner,
  resampling = inner_loop,
  measure = msr("classif.auc"),
  terminator = trm("evals", n_evals = 20),
  tuner = tnr("random_search"),
  store_models = TRUE)


rr = resample(task = task, learner = at, resampling = outer_loop, store_models = TRUE)
extract_inner_tuning_results(rr)

#Unbiased prediction performance (testing error)
rr$aggregate()
 

rr$aggregate(msr("classif.auc"))





