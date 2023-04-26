#R script for comparison of methods to relate a latent class exposure to distal outcomes within a counterfactual mediation model
#methods: one-step, bch, modal class, non-inclusive PCD, inclusive PCD and updated PCD (new method)

#additional files needed to run this script:
#mplus data file: sim'i'.dat (generated using Mplus input files: "1a poor entropy sim data.inp", "1b medium entropy sim data.inp", "1c good entropy sim data.inp")
#mplus input files: 
#"2b uncond latent class.inp"
#"2c onestep mediation.inp"
#"2d bch mediation.inp"
#"2e modal mediation.inp"
#"2f npcd mediation.inp"
#"2g inc latent class.inp"
#"2h incpcd mediation.inp"
#"2i upcd mediation.inp"

#using Mplus v8.8 (display order of results in Mplus output files can depend on version)

############### INSTALL AND LOAD PACKAGES ###################

#install the required packages (if not already installed)
list.of.packages <- c("MplusAutomation", "logistf", "mnormt", "gdata")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] 
if(length(new.packages)) install.packages(new.packages)

#call all required packages
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages)

#check versions of required packages
getNamespaceVersion("MplusAutomation") #1.1.0
getNamespaceVersion("logistf") #1.24.1     
getNamespaceVersion("mnormt") #2.0.2
getNamespaceVersion("gdata") #2.18.0.1

################### POOR ENTROPY ############################

############### DEFINE DIRECTORIES ##########################

#Folder where simulation data are saved:
sim.data <- "file location/all sims/data/poor entropy"
#Folder where Mplus input files (uncond, onestep, bch, modal and inc) are saved:
inp.files <- "file location/all sims"
#*Folder where non-inclusive PCD files are saved:
npcd.files <- "file location/all sims"
#*Folder where inclusive PCD files are saved:
incpcd.files <- "file location/all sims"
#*Folder where updated PCD files are saved:
upcd.files <- "file location/all sims"
#*Folder where results are saved:
results.files <- "file location/all sims/results/poor entropy"

#script for simulated datasets with poor entropy (based on ALSPAC)
#############################################################
#X = 4 class latent nominal exposure (conduct trajectories: early onset persistent, adolescent onset, childhood limited and low)
latent.classes <- 4
#M = binary mediator (peer deviance)
#Y = binary outcome (problematic alcohol use)
#U1-U5 = binary latent class indicators (conduct problems from age 4 to 13 years)
indicators <- 5

############### PREPARE MATRICES ##########################

#here we are running 500 simulated datasets
sims <- 500
#simulated data has a sample size of 5000
sample <- 5000
#we wish to compare each latent class with all the others 
class.comp <- 4*3
#we wish to report 3 mediation effects (total, indirect - tie, and direct - pde) and their standard errors
mediation.effects <- 3*2

#create matrices to store the mediation results for each method
#number of rows = number of simulated datasets (here we use 500)
#number of columns is based on number of class comparisons (12), number of mediation effects (total, tie and pde) and SEs (6) and 
#number of extra parameters (e.g. area under the trajectory, class probabilities, or additional checks)

onestep.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes*2+1))

#labels represent the mediation effects (total, tie and pde) for each class comparison, area under the trajectory for each class, and class probabilities
colnames(onestep.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                           "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                           "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                           "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                           "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se", 
                           "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se", 
                           "auc1","auc2","auc3","auc4","largest.th.se","p1","p2","p3","p4")

bch.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(bch.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                       "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                       "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                       "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                       "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                       "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                       "p1","p2","p3","p4")

modal.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(modal.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                         "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                         "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                         "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                         "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                         "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                         "p1","p2","p3","p4")

npcd.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(npcd.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                        "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                        "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                        "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                        "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                        "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                        "p1","p2","p3","p4")

incpcd.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(incpcd.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                          "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                          "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                          "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                          "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                          "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                          "p1","p2","p3","p4")

pcd.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes+6))

#labels represent the mediation effects (total, tie and pde) for each class comparison, 6 flags for potential issues and class probabilities 
colnames(pcd.all) <-  c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                        "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                        "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                        "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                        "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                        "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se", "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp",
                        "p1","p2","p3","p4")

#create matrix to store the entropy of the unconditional model
entropy <- matrix(NA,sims,1)

#create matrix to store the area under the trajectory for each class in the unconditional model (this allows us to know the order of the classes)
auc <- matrix(NA,sims,latent.classes)
colnames(auc) <- c("auc1","auc2","auc3","auc4")

#create matrix to store the area under the trajectory for each class in the inclusive latent class model (this allows us to know the order of the classes)
auc.inc <- matrix(NA,sims,latent.classes)
colnames(auc.inc) <- c("auc1","auc2","auc3","auc4")

############### READ IN THE DATA ##########################

#run loop over 500 simulated datasets 
for(i in 1:sims) {
  
#set working directory using file paths saved at the start
setwd(sim.data)

#read in Mplus .dat file with simulated data (sim1.dat)
data.original <- read.table(file=paste0("sim", i, ".dat"), sep="", header=FALSE, col.names = c("y", #outcome
                                                                                               "u1","u2","u3","u4","u5", #latent class indicators
                                                                                               "m", #mediator
                                                                                               "c" #modal class assignment for the exposure
                                                                                                ))

#replicate the original data and add empty columns to store imputed latent class membership for 4 classes (once generated) 
#and interactions between each latent class and the mediator
all.data <- data.original
all.data$x1 <- NA #empty column to add dummy code for membership in latent class 1
all.data$x2 <- NA #empty column to add dummy code for membership in latent class 2
all.data$x3 <- NA #empty column to add dummy code for membership in latent class 3
all.data$x4 <- NA #empty column to add dummy code for membership in latent class 4
all.data$int1 <- NA  #empty column to add interaction between latent class 1 and mediator
all.data$int2 <- NA  #empty column to add interaction between latent class 2 and mediator
all.data$int3 <- NA  #empty column to add interaction between latent class 3 and mediator
all.data$int4 <- NA  #empty column to add interaction between latent class 4 and mediator

############### UNCONDITIONAL LATENT CLASS MODEL ##########################

#set working directory using file paths saved at the start
setwd(inp.files)

#this is simply renaming the Mplus .dat file and saving in the folder with the Mplus input files
#this step is important when running many simulations, but not necessary otherwise
#sim.dat will be written over each time a new simulated dataset is analysed
prepareMplusData(data.original,"sim.dat")

#run the unconditional latent class model
runModels("2b uncond latent class.inp")

#read in the Mplus output file and name it "model_output" to use later
model_output <- readModels("2b uncond latent class.out")

#save parameters from unconditional latent class model (within-class thresholds for latent class indicators and class intercepts) and their SEs
coef.se.original <- model_output$parameters$`unstandardized`[1:(indicators*latent.classes+(latent.classes-1)),3:4] 
#labels represent five within-class thresholds across 4 classes and 3 class intercepts)
rownames(coef.se.original) <- c("th1.c1","th2.c1","th3.c1","th4.c1","th5.c1", #within-class thresholds for class 1
                                "th1.c2","th2.c2","th3.c2","th4.c2","th5.c2", #within-class thresholds for class 2
                                "th1.c3","th2.c3","th3.c3","th4.c3","th5.c3", #within-class thresholds for class 3
                                "th1.c4","th2.c4","th3.c4","th4.c4","th5.c4", #within-class thresholds for class 4
                                "int.c1","int.c2","int.c3") #latent class intercepts

#save model entropy into the matrix we created at the start of the script
entropy[i,1] <- model_output$summaries$`Entropy`

#create a matrix of class probabilities
p.original <- matrix(NA,latent.classes,1)
#these can be calculated using the 3 class intercepts we saved in "coef.se.original"
p.original[1,1] <- exp(coef.se.original["int.c1","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
p.original[2,1] <- exp(coef.se.original["int.c2","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
p.original[3,1] <- exp(coef.se.original["int.c3","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))  
p.original[4,1] <- 1/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))

#derive area under the trajectory parameters so we know the order of the latent classes in the unconditional model
#these can be calculated using the within-class thresholds we saved in "coef.se.original"
auc[i,"auc1"] <- exp(-1*coef.se.original["th1.c1","est"])/(1+exp(-1*coef.se.original["th1.c1","est"]))+2*exp(-1*coef.se.original["th2.c1","est"])/(1+exp(-1*coef.se.original["th2.c1","est"]))+3*exp(-1*coef.se.original["th3.c1","est"])/(1+exp(-1*coef.se.original["th3.c1","est"]))+4*exp(-1*coef.se.original["th4.c1","est"])/(1+exp(-1*coef.se.original["th4.c1","est"]))+5*exp(-1*coef.se.original["th5.c1","est"])/(1+exp(-1*coef.se.original["th5.c1","est"]))
auc[i,"auc2"] <- exp(-1*coef.se.original["th1.c2","est"])/(1+exp(-1*coef.se.original["th1.c2","est"]))+2*exp(-1*coef.se.original["th2.c2","est"])/(1+exp(-1*coef.se.original["th2.c2","est"]))+3*exp(-1*coef.se.original["th3.c2","est"])/(1+exp(-1*coef.se.original["th3.c2","est"]))+4*exp(-1*coef.se.original["th4.c2","est"])/(1+exp(-1*coef.se.original["th4.c2","est"]))+5*exp(-1*coef.se.original["th5.c2","est"])/(1+exp(-1*coef.se.original["th5.c2","est"]))
auc[i,"auc3"] <- exp(-1*coef.se.original["th1.c3","est"])/(1+exp(-1*coef.se.original["th1.c3","est"]))+2*exp(-1*coef.se.original["th2.c3","est"])/(1+exp(-1*coef.se.original["th2.c3","est"]))+3*exp(-1*coef.se.original["th3.c3","est"])/(1+exp(-1*coef.se.original["th3.c3","est"]))+4*exp(-1*coef.se.original["th4.c3","est"])/(1+exp(-1*coef.se.original["th4.c3","est"]))+5*exp(-1*coef.se.original["th5.c3","est"])/(1+exp(-1*coef.se.original["th5.c3","est"]))
auc[i,"auc4"] <- exp(-1*coef.se.original["th1.c4","est"])/(1+exp(-1*coef.se.original["th1.c4","est"]))+2*exp(-1*coef.se.original["th2.c4","est"])/(1+exp(-1*coef.se.original["th2.c4","est"]))+3*exp(-1*coef.se.original["th3.c4","est"])/(1+exp(-1*coef.se.original["th3.c4","est"]))+4*exp(-1*coef.se.original["th4.c4","est"])/(1+exp(-1*coef.se.original["th4.c4","est"]))+5*exp(-1*coef.se.original["th5.c4","est"])/(1+exp(-1*coef.se.original["th5.c4","est"]))

#import tech3 as a covariance matrix to allow us to use the covariance between parameters for perturbing in updated PCD
#############################################################
#NB the order of parameters in this covariance matrix is different to order in "coef.se.original" (see tech1) so this will need to be addressed before perturbing parameters
#this is due to no missing data in the indicators which affects the numbering of parameters in tech1
#############################################################
#this outputs paramCov as a matrix object in R 
tech3 <- model_output$tech3$paramCov
#this creates a full, symmetrical covariance matrix
upperTriangle(tech3) <- lowerTriangle(tech3, byrow=TRUE)
cov <- tech3

#address any fixed parameters: within class thresholds that have been fixed at 15 or -15 (representing a probability of 0 or 100%) do not have (co)variances
#the code below means that fixed parameters will still get perturbed in updated PCD but only a very small amount
#replace missing (999) in covariance matrix with 0 to represent no covariance for the fixed parameters
cov[cov==999] <- 0
#change the variance for fixed parameters to be very small (0.000000001)
for(j in 1:(indicators*latent.classes+(latent.classes-1))) {
  if (cov[j,j]==0) cov[j,j]<- 0.000000001}

###################### ONESTEP MODEL ##########################

#run the onestep latent class mediation model
runModels("2c onestep mediation.inp")

#read in parameters from one-step latent class model    
est.onestep <- readModels("2c onestep mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects, area under the trajectory for each class, and class probabilities) and their SEs in the matrix created at the start "onestep.all"
onestep.all[i, (match("tot.1v4",colnames(onestep.all))):(match("pde.4v1",colnames(onestep.all)))] <- est.onestep[116:151,"est"] # TOT_1V4 to PDE_4V1
onestep.all[i, (match("tot.1v4.se",colnames(onestep.all))):(match("pde.4v1.se",colnames(onestep.all)))] <- est.onestep[116:151,"se"] # TOT_1V4 to PDE_4V1
onestep.all[i, (match("auc1",colnames(onestep.all))):(match("auc4",colnames(onestep.all)))] <- est.onestep[152:155,"est"] # AUC
onestep.all[i, (match("p1",colnames(onestep.all))):(match("p4",colnames(onestep.all)))] <- est.onestep[88:91,"est"]  # P_X1 to P_X4

#save standard errors for thresholds from onestep model
se <- est.onestep[c(3:7,17:21,31:35,45:49),4] 

#we want a flag so that we know the largest standard error for a threshold in the onestep model
#this flag captures the largest SE
onestep.all[i,"largest.th.se"] <- max(se)

######################## BCH ##########################

#bch can only be used with multiple latent classes if they are combined into one class
#use one 8 class model for X and M and set up mediation model using nom-nom-cat approach (because bch needs to include an auxiliary variable)

#read in bch weights and modal class assignment which were exported out of the unconditional latent class model in "bch.txt"
data.bch <- read.table("bch.txt", sep="", na.strings="*", header=FALSE)
#select the weights for each class and the modal class assignment
data.bch <- data.bch[,c(6:9,14)]
#labels representing one weight for each class and the modal class assignment
colnames(data.bch) <- c("bch1","bch2","bch3","bch4","modal")
#combine the bch weights with the mediator and outcome from the original data
data.bch <- cbind(data.original[,c("y","m")], data.bch)

#create weights for 8 class model by multiplying the bch weight for each class with the observed data on the mediator
data.bch$bch11 <- data.bch$bch1*data.bch$m
data.bch$bch12 <- data.bch$bch1*(1-data.bch$m)
data.bch$bch21 <- data.bch$bch2*data.bch$m
data.bch$bch22 <- data.bch$bch2*(1-data.bch$m)
data.bch$bch31 <- data.bch$bch3*data.bch$m
data.bch$bch32 <- data.bch$bch3*(1-data.bch$m)
data.bch$bch41 <- data.bch$bch4*data.bch$m
data.bch$bch42 <- data.bch$bch4*(1-data.bch$m)

#prepare Mplus .dat file using "data.bch"  
prepareMplusData(data.bch,"bch.dat")
#run mediation model using bch method
runModels("2d bch mediation.inp")

#read in parameters from bch latent class mediation model   
est.bch <- readModels("2d bch mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "bch.all"
bch.all[i,(match("tot.1v4",colnames(bch.all))):(match("pde.4v1",colnames(bch.all)))] <- est.bch[68:103,"est"] # TOT_1V4 to PDE_4V1
bch.all[i,(match("tot.1v4.se",colnames(bch.all))):(match("pde.4v1.se",colnames(bch.all)))] <- est.bch[68:103,"se"] # TOT_1V4 to PDE_4V1
bch.all[i,(match("p1",colnames(bch.all))):(match("p4",colnames(bch.all)))] <- est.bch[40:43,"est"]  # P_X1 to P_X4

######################## MODAL CLASS ASSIGNMENT ##########################

#run mediation model using modal class assignment
runModels("2e modal mediation.inp")

#read in parameters from modal class mediation model  
est.modal <- readModels("2e modal mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "modal.all"
modal.all[i,(match("tot.1v4",colnames(modal.all))):(match("pde.4v1",colnames(modal.all)))] <- est.modal[100:135,"est"] # TOT_1V4 to PDE_4V1
modal.all[i,(match("tot.1v4.se",colnames(bch.all))):(match("pde.4v1.se",colnames(bch.all)))] <- est.modal[100:135,"se"] # TOT_1V4 to PDE_4V1
modal.all[i,(match("p1",colnames(bch.all))):(match("p4",colnames(bch.all)))] <- est.modal[72:75,"est"]  # P_X1 to P_X4

######################## NON-INCLUSIVE PCD ##########################

#set seed for non-inclusive PCD
#this will set the seed for each session in R (not each time run loops)
#to make sure you get the same results each time, run the loop only once within each R session
set.seed(80)

#calculate P(X=x|U) - this will give us the probability of class membership (cprobs) from unconditional latent class model

#for each class, multiply the individual data (responses to 5 binary indicators: U1 to U5) with within-class thresholds
#number of latent class indicators
#class 1
#this creates a matrix with within-class thresholds for class 1 repeated for every individual in the dataset (e.g. repeated 5000 times)
theta <- matrix(rep(coef.se.original[(match("th1.c1",rownames(coef.se.original))):(match("th5.c1",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
#if latent class indicator is present (-1)*threshold is used, if latent class indicator is absent (1)*threshold is used
P1<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
#class 2
#this creates a matrix with within-class thresholds for class 2 repeated for every individual in the dataset (e.g. repeated 5000 times)
theta <- matrix(rep(coef.se.original[(match("th1.c2",rownames(coef.se.original))):(match("th5.c2",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
P2<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
#class 3
#this creates a matrix with within-class thresholds for class 3 repeated for every individual in the dataset (e.g. repeated 5000 times)
theta <- matrix(rep(coef.se.original[(match("th1.c3",rownames(coef.se.original))):(match("th5.c3",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
P3<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
#class 4
#this creates a matrix with within-class thresholds for class 4 repeated for every individual in the dataset (e.g. repeated 5000 times)
theta <- matrix(rep(coef.se.original[(match("th1.c4",rownames(coef.se.original))):(match("th5.c4",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
P4<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))

#multiply the 5 probabilities (for each latent class indicator) in each row of P* then multiply by class probabilities that were saved earlier in "p.original"  
N1<-apply(P1,1,prod)*p.original[1,1]
N2<-apply(P2,1,prod)*p.original[2,1]
N3<-apply(P3,1,prod)*p.original[3,1]
N4<-apply(P4,1,prod)*p.original[4,1]    

#this gives us the probability of class membership for each person in the dataset 
P<-cbind(N1/(N1+N2+N3+N4),N2/(N1+N2+N3+N4),N3/(N1+N2+N3+N4),N4/(N1+N2+N3+N4))

#we could have simply read these in from "bch.txt" (as we did with the bch weights) as they were exported from the unconditional latent class model 
colnames(P) <-  c("cprob1","cprob2","cprob3","cprob4")

#we will now impute class membership 40 times for each person using their probability of class membership
#we will create 40 imputed datasets (chosen to keep Monte Carlo error at less than 10% of standard error for parameters from regression model for Y)
imp.n <- 40
#create a matrix to store imputed class membership for each person (sample*imp.n=5000*40)
imp <- matrix(NA,sample*imp.n,latent.classes+1)
colnames(imp) <-  c("imp","x1","x2","x3","x4")

#first column is simply an indicator for imputation number (range from 1 to 40)
for(h in 0:(imp.n-1)) {
  imp [c(sample*h+1:sample),1]<-h+1}

#for each of 40 imputations we will create a matrix to store class membership, 
#then we will use the probabilities of class membership "P" to randomly assign each individual to a class (X = 1,..k).
#and we will add their imputed class membership into the matrix "imp" we created earlier
for(j in 1:imp.n) {
  x <- matrix(NA,sample,latent.classes)    
  for(k in 1:sample) {x[k,1:latent.classes] <- rmultinom(1,1,P[k,])}
  for(l in 1:imp.n) {
    if (j==l) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x}
} 

#analysis
#can use a separate folder for non-inclusive PCD because of large number of files created (e.g., 40 imputed datasets)
setwd(npcd.files) 

#combine imputed class membership stored in "imp.n" with original data
#prepare 1 mplus .dat file for each imputed dataset to run mediation model (40 .dat files should be created)
for(l in 1:imp.n) {
  imp.subset <- cbind(data.original,subset(imp, imp[,1]==l))
  prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))
} 

#create "imp.txt" file for mplus to call imputed datasets
imp.txt <- matrix(NA,imp.n,1)
for(l in 1:imp.n) {
  imp.txt[l,1] <- paste0("imp_", l, ".dat")
} 
write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#run mediation model using non-inclusive PCD
runModels("2f npcd mediation.inp")

#read in parameters from npcd mediation model 
est.npcd <- readModels("2f npcd mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "npcd.all"
npcd.all[i,(match("tot.1v4",colnames(npcd.all))):(match("pde.4v1",colnames(npcd.all)))] <- est.npcd[100:135,"est"] # TOT_1V4 to PDE_4V1
npcd.all[i,(match("tot.1v4.se",colnames(npcd.all))):(match("pde.4v1.se",colnames(npcd.all)))] <- est.npcd[100:135,"se"] # TOT_1V4 to PDE_4V1
npcd.all[i,(match("p1",colnames(npcd.all))):(match("p4",colnames(npcd.all)))] <- est.npcd[72:75,"est"]  # P_X1 to P_X4

######################## INCLUSIVE PCD ##########################

#set seed for inclusive PCD
set.seed(81)

#set working directory using file paths saved at the start
setwd(inp.files)

#run inclusive latent class model and export out cprobs
runModels("2g inc latent class.inp")

#read in auc parameters from inclusive latent class model
est.inc <- readModels("2g inc latent class.out", what="parameters")$parameters$`unstandardized`
auc.inc[i,(match("auc1",colnames(auc.inc))):(match("auc4",colnames(auc.inc)))] <- est.inc[30:33,"est"]  # AUC1 to AUC4

#read in class probabilities from "inc.txt" (as we did with the bch weights exported from the unconditional latent class model) 
cprobs <- read.table("inc.txt", sep="", na.strings="*", header=FALSE)
#we could also extract modal classes from here (v12) to use the inclusive modal approach
P <- cprobs[,c(8:11)]
colnames(P) <- c("cprob1","cprob2","cprob3","cprob4")

#we will now impute class membership 40 times for each person using their probability of class membership
#create a matrix to store imputed class membership for each person (sample*imp.n=5000*40)
imp <- matrix(NA,sample*imp.n,latent.classes+1)
colnames(imp) <-  c("imp","x1","x2","x3","x4")

#first column is simply an indicator for imputation number (range from 1 to 40)
for(h in 0:(imp.n-1)) {
  imp [c(sample*h+1:sample),1]<-h+1}

#for each of 40 imputations we will create a matrix to store class membership, 
#then we will use the probabilities of class membership "P" to randomly assign each individual to a class (X = 1,..k).
#and we will add their imputed class membership into the matrix "imp" we created earlier
for(j in 1:imp.n) {
  x <- matrix(NA,sample,latent.classes)    
  for(k in 1:sample) {x[k,1:latent.classes] <- rmultinom(1,1,P[k,])}
  for(l in 1:imp.n) {
    if (j==l) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x}
} 

#analysis
#can use a separate folder for inclusive PCD because of large number of files created (e.g., 40 imputed datasets)
setwd(incpcd.files) 

#combine imputed class membership stored in "imp.n" with original data
#prepare 1 mplus .dat file for each imputed dataset to run mediation model (40 .dat files should be created)
for(l in 1:imp.n) {
  imp.subset <- cbind(data.original,subset(imp, imp[,1]==l))
  prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))
} 

#create "imp.txt" file for mplus to call imputed datasets
imp.txt <- matrix(NA,imp.n,1)
for(l in 1:imp.n) {
  imp.txt[l,1] <- paste0("imp_", l, ".dat")
} 
write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#run mediation model using inclusive PCD
runModels("2h incpcd mediation.inp")

#read in parameters from incpcd mediation model 
est.incpcd <- readModels("2h incpcd mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "incpcd.all"
incpcd.all[i,(match("tot.1v4",colnames(incpcd.all))):(match("pde.4v1",colnames(incpcd.all)))] <- est.incpcd[100:135,"est"] # TOT_1V4 to PDE_4V1
incpcd.all[i,(match("tot.1v4.se",colnames(incpcd.all))):(match("pde.4v1.se",colnames(incpcd.all)))] <- est.incpcd[100:135,"se"] # TOT_1V4 to PDE_4V1
incpcd.all[i,(match("p1",colnames(incpcd.all))):(match("p4",colnames(incpcd.all)))] <- est.incpcd[72:75,"est"]  # P_X1 to P_X4

#remove any objects from memory that will not be needed for the rest of the script
rm(model_output, est.onestep, data.bch, est.bch, est.modal, est.npcd, est.incpcd, cprobs, P)

######################## UPDATED PCD ##########################

#set seed for updated PCD
set.seed(82)

#STEP 1: 
#initialise beta by fitting the logistic regression model for Y with beta for classes set to zero 
#initialise alpha by fitting the logistic regression model for M with alpha for classes set to zero
#perturb the estimates from these models with Gaussian noise of mean zero and values from the variance-covariance matrix of these parameter estimates

#firth logistic regression is used to address problems with perfect prediction when bringing in exposure-mediator interactions later in the script
#regress outcome y on mediator m
model.y <- logistf(y~m, data=all.data)
#regression model for mediator m 
model.m <- logistf(m~1, data=all.data)

#perturbing the coefficients once around coefficients in models using variance-covariance matrix 
#initially use zero for exposure coefficients (and exposure-mediator interactions in outcome model)
#once latent class exposure has been imputed in subsequent runs, these will become coefficients from models
beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y)),0,0,0,0,0,0)
beta.m<-c(rmnorm(1,coef(model.m),vcov(model.m)),0,0,0)

#we will create 80 imputed datasets for class membership (chosen to keep Monte Carlo error at less than 10% of standard error for parameters in regression model for Y)
imp.n <- 80
#we will allow 20 iterations between saving out imputed class membership
cycles <- 20
#we will allow a burn in of 100 iterations before starting to save out imputed class membership
burnin <- 100

#create a matrix to save the results from every iteration (all results are recorded to assess convergence later)
#rows=iterations=(cycles*imp.n+burnin)
#columns=iteration number, beta coeficients from regression models, cell sizes and 2 flags for issues=(6+(latent.classes-1)*3+latent.classes*4)=31
results <- matrix(NA,cycles*imp.n+burnin,6+(latent.classes-1)*3+latent.classes*4) 

colnames(results) <- c("iteration", #iteration number
                       "b0.y","b1.y","b2.y","b3.y","b4.y","b5.y","b6.y","b7.y", #coefficients from regression model for the outcome
                       "b0.m","b1.m","b2.m","b3.m", #coefficients from regression model for the mediator
                       "x1.m0.y0","x1.m1.y0","x2.m0.y0","x2.m1.y0","x3.m0.y0","x3.m1.y0","x4.m0.y0","x4.m1.y0","x1.m0.y1","x1.m1.y1","x2.m0.y1","x2.m1.y1","x3.m0.y1","x3.m1.y1","x4.m0.y1","x4.m1.y1", #cell sizes from crosstabs for classes by mediator by outcome
                       "zero.cell", #flag for presence of zero cells in crosstabs
                       "large.th" #flag for a within-class threshold in unconditional latent class model that was out of bounds after perturbing
                        ) 

#the first column is simply an indicator of iteration number (range from 1 to 1700)
results[,1]<-1:(cycles*imp.n+burnin)

#create a matrix to store cell sizes from crosstabs for classes by mediator by outcome 
xmy<-matrix(NA,1,latent.classes*4)
colnames(xmy) <-  c("x1m0y0","x1m1y0","x2m0y0","x2m1y0","x3m0y0","x3m1y0","x4m0y0","x4m1y0","x1m0y1","x1m1y1","x2m0y1","x2m1y1","x3m0y1","x3m1y1","x4m0y1","x4m1y1")
#create a matrix to store presence of zero cells in this crosstabs
zero.cell<-matrix(NA,1,1)
#create a matrix to store presence of a within-class threshold in unconditional latent class model that was out of bounds after perturbing (e.g., not corresponding to 0 to 100% probability) 
large.th<-matrix(NA,1,1)

#create a matrix to store imputed class membership for each person (sample*imp.n=5000*80)
imp <- matrix(NA,sample*imp.n,latent.classes+1)
colnames(imp) <- c("imp","x1","x2","x3","x4") 
#first column is simply an indicator for imputation number (range from 1 to 80)
for(h in 0:(imp.n-1)) {
  imp [c(sample*h+1:sample),1]<-h+1}

#we will now create a loop to repeat steps 2 and 3 below 1700 (imp.n*cycles+burnin) times and save estimates after every 20 iterations (after a 100 iteration burn in)
for(j in 1:(imp.n*cycles+burnin)) {
  results[j,2:9]<-beta.y #save the perturbed coefficients from regression model for the outcome "beta.y"
  results[j,10:13]<-beta.m #save the perturbed coefficients from regression model for the mediator "beta.m"    
  results[j,14:29]<-xmy #save cell sizes from crosstabs for classes by mediator by outcome "xmy"
  results[j,30]<-zero.cell #save a flag for presence of zero cells in this crosstabs "zero.cell"
  results[j,31]<-large.th #save a flag for a within-class threshold in unconditional latent class model that was out of bounds after peturbing "large.th"
  
  #this saves each individual's imputed class membership 80 times (generated later in the script) into the matrix "imp" we created earlier     
  for(l in 1:imp.n) {
    if (j==(cycles*l+burnin)) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x
  } 
  
  #STEP 2a:
  #we will now perturb the parameters (within-class thresholds and class intercepts) that we saved earlier "coef.se.original" from the unconditional latent class model 
  #in order to perturb the parameters using the covariance matrix in tech3 we need to reorder parameters so they match with numbering in tech1
  #because we have no missing data, they are in a different order to what would be expected
  #when using a dataset with missing data on class indicators (as would usually be the case outside of simulated data), this reordering step is not needed
  
  #keep only the parameters (within-class thresholds for latent class indicators and class intercepts) from unconditional latent class model (i.e. drop the SEs)
  coef.mplus.orig<-coef.se.original[,"est"]
  #reorder the parameters so that they are in the same order as is used in the covariances matrix of parameters from tech3
  index<-c(1,5,9,13,17,2,6,10,14,18,3,7,11,15,19,4,8,12,16,20,21,22,23)
  coef.mplus.orig<-coef.mplus.orig[order(index)]
  #perturb these parameters based on their variance-covariance matrix (saved in "cov")
  coef.mplus<-rmnorm(1,coef.mplus.orig,cov)
  #reorder again to preserve original ordering
  index<-c(1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19,5,10,15,20,21,22,23)
  coef<-coef.mplus[order(index)]
  
  #when perturbing within-class thresholds, those with a large standard error can go out of bounds (e.g., corresponding to a probability that is not between 0 and 100%)
  #we will create a flag so we know when this is the case
  large.th[1,1]<-0    
  for(k in 1:(latent.classes*indicators)) {
    if (coef[k]>15||coef[k]<(-15)) large.th[1,1] <- 1} 
  #we will also constrain within-class thresholds to be between -15 and 15 (corresponding to a probability that is between 0 and 100%)
  coef[coef>15] <- 15
  coef[coef<(-15)] <- (-15)
  
  #turning 'coef' into a matrix so can label rows and columns to use in script below (instead of referring to numbers)
  coef<-as.matrix(coef)
  rownames(coef) <- c("th1.c1","th2.c1","th3.c1","th4.c1","th5.c1",
                      "th1.c2","th2.c2","th3.c2","th4.c2","th5.c2",
                      "th1.c3","th2.c3","th3.c3","th4.c3","th5.c3",
                      "th1.c4","th2.c4","th3.c4","th4.c4","th5.c4",
                      "int.c1","int.c2","int.c3")
  colnames(coef) <- c("est")
  
  #now calculate P(X=x|U) - without perturbing this give us the probability of class membership (cprobs) from unconditional latent class model
  
  #for each class, multiply the individual data (responses to 5 binary indicators: U1 to U5) with within-class thresholds
  #class 1
  #this creates a matrix with within-class thresholds for class 1 repeated for every individual in the dataset (e.g. repeated 5000 times)
  theta <- matrix(rep(coef[(match("th1.c1",rownames(coef.se.original))):(match("th5.c1",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
  #if latent class indicator is present (-1)*threshold is used, if latent class indicator is absent (1)*threshold is used
  P1<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  #class 2
  theta <- matrix(rep(coef[(match("th1.c2",rownames(coef.se.original))):(match("th5.c2",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
  P2<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  #class 3
  theta <- matrix(rep(coef[(match("th1.c3",rownames(coef.se.original))):(match("th5.c3",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
  P3<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  #class 4
  theta <- matrix(rep(coef[(match("th1.c4",rownames(coef.se.original))):(match("th5.c4",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
  P4<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  
  #we need to create a new matrix of class probabilities using the perturbed class intercepts
  p <- matrix(NA,latent.classes,1)
  p[1,1] <- exp(coef["int.c1","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
  p[2,1] <- exp(coef["int.c2","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
  p[3,1] <- exp(coef["int.c3","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))  
  p[4,1] <- 1/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
  
  #then multiply "P1" to "P4" by the perturbed class probabilities that we have created above "p" 
  #this gives us the probability of class membership for each person in the dataset
  #these will differ slightly to "cprobs" that can be exported from the unconditional latent class model due to perturbing
  N1<-apply(P1,1,prod)*p[1,1]
  N2<-apply(P2,1,prod)*p[2,1]
  N3<-apply(P3,1,prod)*p[3,1]
  N4<-apply(P4,1,prod)*p[4,1]    
  
  #STEP 2b - Bayes rule:
  
  #combine P(X=x|U) along with coefficients from regression models for outcome and mediator
  #this will calculate probabilities P(X=x|Y,M,U) for each class (x = 1,..k) for each individual
  #we will do this via a number a steps below
  
  #P(Y=1|M,X): calculate probability that the outcome (Y) = 1 given mediator (M) and exposure (X, latent classes) 
  #this uses the coefficients from the regression model for the outcome (intercept (b0), coef for M (b1), coefs for X (b2-b4) and coefs for XM interaction (b5-b7))
  #class 1
  num1 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m+results[j,"b2.y"]+results[j,"b5.y"]*all.data$m)
  Y11 <- num1/(1+num1)  
  #class 2
  num2 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m+results[j,"b3.y"]+results[j,"b6.y"]*all.data$m)
  Y21 <- num2/(1+num2)  
  #class 3
  num3 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m+results[j,"b4.y"]+results[j,"b7.y"]*all.data$m)
  Y31 <- num3/(1+num3)
  #class 4
  num4 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m)
  Y41 <- num4/(1+num4)
  
  #P(Y=0|M,X): calculate probability that the outcome (Y) = 0 given mediator (M) and exposure (X, latent classes) 
  Y10 <- 1/(1+num1)
  Y20 <- 1/(1+num2)
  Y30 <- 1/(1+num3)
  Y40 <- 1/(1+num4)
  
  #P(M=1|X): calculate probability that the mediator (M) = 1 given exposure (X, latent classes) 
  #this uses the coefficients from the regression model for the mediator (intercept (b0) and coefs for X (b1-b3))
  #class 1 
  num1 <- exp(results[j,"b0.m"]+results[j,"b1.m"])
  M11 <- num1/(1+num1)
  #class 2
  num2 <- exp(results[j,"b0.m"]+results[j,"b2.m"])
  M21 <- num2/(1+num2) 
  #class 3
  num3 <- exp(results[j,"b0.m"]+results[j,"b3.m"])
  M31 <- num3/(1+num3)   
  #class 4
  num4 <- exp(results[j,"b0.m"])
  M41 <- num4/(1+num4) 
  
  #P(M=0|X): calculate probability that the mediator (M) = 0 given exposure (X, latent classes)    
  M10 <- 1/(1+num1)
  M20 <- 1/(1+num2)
  M30 <- 1/(1+num3)      
  M40 <- 1/(1+num4)
  
  #P(X=x|Y=1,M=1,U)
  #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
  N111<-N1*Y11*M11
  #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
  N211<-N2*Y21*M21
  #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
  N311<-N3*Y31*M31
  #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
  N411<-N4*Y41*M41
  
  denom <- N111+N211+N311+N411
  Q111<-N111/denom
  Q211<-N211/denom
  Q311<-N311/denom
  Q411<-N411/denom 
  Q11<-cbind(Q111,Q211,Q311,Q411)
  head(Q11)
  
  #P(X=x|Y=0,M=0,U) 
  #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
  N100<-N1*Y10*M10
  #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
  N200<-N2*Y20*M30
  #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
  N300<-N3*Y30*M30
  #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
  N400<-N4*Y40*M40
  
  denom <- N100+N200+N300+N400   
  Q100<-N100/denom
  Q200<-N200/denom
  Q300<-N300/denom
  Q400<-N400/denom    
  Q00<-cbind(Q100,Q200,Q300,Q400)
  head(Q00)
  
  #P(X=x|Y=1,M=0,U)
  #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
  N101<-N1*Y11*M10
  #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
  N201<-N2*Y21*M30
  #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
  N301<-N3*Y31*M30
  #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
  N401<-N4*Y41*M40
  
  denom <- N101+N201+N301+N401     
  Q101<-N101/denom
  Q201<-N201/denom
  Q301<-N301/denom
  Q401<-N401/denom    
  Q01<-cbind(Q101,Q201,Q301,Q401)
  head(Q01)    
  
  #P(X=x|Y=0,M=1,U)
  #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
  N110<-N1*Y10*M11
  #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
  N210<-N2*Y20*M31
  #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
  N310<-N3*Y30*M31
  #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
  N410<-N4*Y40*M41
  
  denom <- N110+N210+N310+N410    
  Q110<-N110/denom
  Q210<-N210/denom
  Q310<-N310/denom
  Q410<-N410/denom   
  Q10<-cbind(Q110,Q210,Q310,Q410)
  head(Q10)   
  
  #derive the probability of class membership for each individual which takes into account the relationship between the classes, mediator and outcome
  #using each individuals observed data on the mediator and outcome
  Q<-Q11
  #probabilities for those with mediator and outcome absent
  Q[all.data$y==0 & all.data$m==0,]<-Q00[all.data$y==0 & all.data$m==0,]
  #probabilities for those with mediator absent and outcome present
  Q[all.data$y==1 & all.data$m==0,]<-Q01[all.data$y==1 & all.data$m==0,]  
  #probabilities for those with mediator present and outcome absent
  Q[all.data$y==0 & all.data$m==1,]<-Q10[all.data$y==0 & all.data$m==1,]
  #these probabilities will be used to impute class membership for each individual 
  head(Q)
  
  colnames(Q) <-  c("cprob1","cprob2","cprob3","cprob4")
  
  #STEP 3:
  
  #create a matrix to store class membership 
  x <- matrix(NA,sample,latent.classes)
  #now we will use the probabilities of class membership "Q" to randomly assign each individual to a class (X = 1,..k).
  #n=1 (number of random vectors to draw); size=1 per person
  for(k in 1:sample) {x[k,1:latent.classes] <- rmultinom(1,1,Q[k,])}
  
  #add imputed class membership to the data
  all.data[,"x1"] <- x[,1]
  all.data[,"x2"] <- x[,2]
  all.data[,"x3"] <- x[,3]
  all.data[,"x4"] <- x[,4]
  #add in exposure-mediator interactions
  all.data[,"int1"] <- all.data[,"x1"]*all.data[,"m"]
  all.data[,"int2"] <- all.data[,"x2"]*all.data[,"m"]
  all.data[,"int3"] <- all.data[,"x3"]*all.data[,"m"]
  all.data[,"int4"] <- all.data[,"x4"]*all.data[,"m"]
  
  ###################################################
  #we will now perform some checks on the cell sizes from crosstabs for classes by mediator by outcome 
  #we will add these cell sizes into "results" to make traceplots to assess convergence
  
  #create a subset of the data for Y=0 
  no.y <- subset(all.data, all.data$y==0)
  #create a subset of the data for Y=1
  yes.y <- subset(all.data, all.data$y==1)
  
  #crosstabs for x1 and m for those with y=0
  x1m.no.y<-table(no.y$x1,no.y$m)
  
  #this is to make sure matrix is 2 by 2 even when there are zero cells
  if (nrow(x1m.no.y)==1) x1m.no.y <- rbind(x1m.no.y,matrix(0,1,2))
  #crosstabs for x2 and m for those with y=0 
  x2m.no.y<-table(no.y$x2,no.y$m)
  if (nrow(x2m.no.y)==1) x2m.no.y <- rbind(x2m.no.y,matrix(0,1,2))
  #crosstabs for x3 and m for those with y=0 
  x3m.no.y<-table(no.y$x3,no.y$m)
  if (nrow(x3m.no.y)==1) x3m.no.y <- rbind(x3m.no.y,matrix(0,1,2))
  #crosstabs for x4 and m for those with y=0 
  x4m.no.y<-table(no.y$x4,no.y$m)
  if (nrow(x4m.no.y)==1) x4m.no.y <- rbind(x4m.no.y,matrix(0,1,2))
  #crosstabs for x1 and m for those with y=1 
  x1m.yes.y<-table(yes.y$x1,yes.y$m)
  if (nrow(x1m.yes.y)==1) x1m.yes.y <- rbind(x1m.yes.y,matrix(0,1,2))    
  x2m.yes.y<-table(yes.y$x2,yes.y$m)
  if (nrow(x2m.yes.y)==1) x2m.yes.y <- rbind(x2m.yes.y,matrix(0,1,2))    
  x3m.yes.y<-table(yes.y$x3,yes.y$m)
  if (nrow(x3m.yes.y)==1) x3m.yes.y <- rbind(x3m.yes.y,matrix(0,1,2))
  x4m.yes.y<-table(yes.y$x4,yes.y$m)
  if (nrow(x4m.yes.y)==1) x4m.yes.y <- rbind(x4m.yes.y,matrix(0,1,2))
  
  #use empty xmy matrix created earlier and fill in with cell sizes
  #cell size for x=1, m=0, y=0
  xmy[1,"x1m0y0"] <- x1m.no.y[2,1]
  #cell size for x=1, m=1, y=0  
  xmy[1,"x1m1y0"] <- x1m.no.y[2,2]
  #cell size for x=2, m=0, y=0  
  xmy[1,"x2m0y0"] <- x2m.no.y[2,1]
  #cell size for x=2, m=1, y=0  
  xmy[1,"x2m1y0"] <- x2m.no.y[2,2]
  #cell size for x=3, m=0, y=0  
  xmy[1,"x3m0y0"] <- x3m.no.y[2,1]
  #cell size for x=3, m=1, y=0  
  xmy[1,"x3m1y0"] <- x3m.no.y[2,2]
  #cell size for x=4, m=0, y=0  
  xmy[1,"x4m0y0"] <- x4m.no.y[2,1]
  #cell size for x=4, m=1, y=0  
  xmy[1,"x4m1y0"] <- x4m.no.y[2,2]
  #cell size for x=1, m=0, y=1  
  xmy[1,"x1m0y1"] <- x1m.yes.y[2,1]
  #cell size for x=1, m=1, y=1  
  xmy[1,"x1m1y1"] <- x1m.yes.y[2,2]
  #cell size for x=2, m=0, y=1  
  xmy[1,"x2m0y1"] <- x2m.yes.y[2,1]
  #cell size for x=2, m=1, y=1  
  xmy[1,"x2m1y1"] <- x2m.yes.y[2,2]
  #cell size for x=3, m=0, y=1  
  xmy[1,"x3m0y1"] <- x3m.yes.y[2,1]
  #cell size for x=3, m=1, y=1  
  xmy[1,"x3m1y1"] <- x3m.yes.y[2,2]
  #cell size for x=4, m=0, y=1  
  xmy[1,"x4m0y1"] <- x4m.yes.y[2,1]
  #cell size for x=4, m=1, y=1  
  xmy[1,"x4m1y1"] <- x4m.yes.y[2,2]
  
  #flag for number of zero cells in crosstabs for classes by mediator by outcome
  zero.cell[1,1] <- length(which(xmy == 0))
  ###################################################
  
  #fit the firth logistic regression model for P(Y|X,M) to obtain updated parameter estimates now classes have been imputed
  if (zero.cell<2) model.y <- logistf(y~m+x1+x2+x3+int1+int2+int3, data=all.data)
  #if there is more than 1 zero cell, even firth logistic regression does not converge, therefore it is necessary to remove XM interactions from regression model
  if (zero.cell>1) model.y <- logistf(y~m+x1+x2+x3, data=all.data)
  
  #fit the logistic regression model for P(M|X) to obtain updated parameter estimates
  model.m <- logistf(m~x1+x2+x3, data=all.data)
  
  #perturbing the beta coefficients once around coefficients in models using variance-covariance matrix 
  if (zero.cell>1) beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y)),0,0,0)
  if (zero.cell<2) beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y))) 
  beta.m<-c(rmnorm(1,coef(model.m),vcov(model.m)))
  
} #end of iteration loop

################################################
#assessing convergence

#create plot of parameters and cell sizes across each iteration
results.df <- as.data.frame(results)

#remove the first iteration when no data for X (latent classes) 
results.df <- subset(results.df, iteration>1)

#to check autocorrelation for each beta after burn in of 100 iterations
b2.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b2.y"]
b3.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b3.y"]
b4.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b4.y"]
b5.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b5.y"]
b6.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b6.y"]
b7.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b7.y"]
b1.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b1.m"]
b2.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b2.m"]
b3.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b3.m"]

#can use a separate folder for updated PCD because of large number of files created (e.g., 80 imputed datasets)
setwd(upcd.files)

#create a pdf showing traceplots, histograms and autocorrelation plots for each beta, and traceplots and histograms for cell sizes across all iterations

pdf(file=paste0("trace_plots", i, ".pdf")) 
plot(results.df$iteration, results.df$b2.y,
     xlab = "iteration",
     ylab = "b2.y (y on class 1 vs class 4)")
lines(results.df$iteration, results.df$b2.y)
hist (results.df$b2.y, xlab = "b2.y (y on class 1 vs class 4)")
acf(b2.y, lag.max=50)
plot(results.df$iteration, results.df$b3.y,
     xlab = "iteration",
     ylab = "b3.y (y on class 2 vs class 4)")
lines(results.df$iteration, results.df$b3.y)
hist (results.df$b3.y, xlab = "b3.y (y on class 2 vs class 4)")
acf(b3.y, lag.max=50)
plot(results.df$iteration, results.df$b4.y,
     xlab = "iteration",
     ylab = "b4.y (y on class 3 vs class 4)")
lines(results.df$iteration, results.df$b4.y)
hist (results.df$b4.y, xlab = "b4.y (y on class 3 vs class 4)")
acf(b4.y, lag.max=50)
plot(results.df$iteration, results.df$b5.y,
     xlab = "iteration",
     ylab = "b5.y (y on class 1 x m interaction)")
lines(results.df$iteration, results.df$b5.y)
hist (results.df$b5.y, xlab = "b5.y (y on class 1 x m interaction)")
acf(b5.y, lag.max=50)
plot(results.df$iteration, results.df$b6.y,
     xlab = "iteration",
     ylab = "b6.y (y on class 2 x m interaction)")
lines(results.df$iteration, results.df$b6.y)
hist (results.df$b6.y, xlab = "b6.y (y on class 2 x m interaction)")
acf(b6.y, lag.max=50)
plot(results.df$iteration, results.df$b7.y,
     xlab = "iteration",
     ylab = "b7.y (y on class 3 x m interaction)")
lines(results.df$iteration, results.df$b7.y)
hist (results.df$b7.y, xlab = "b7.y (y on class 3 x m interaction)")
acf(b7.y, lag.max=50)
plot(results.df$iteration, results.df$b1.m,
     xlab = "iteration",
     ylab = "b1.m (m on class 1 vs class 4)")
lines(results.df$iteration, results.df$b1.m)
hist (results.df$b1.m, xlab = "b1.m (m on class 1 vs class 4)")
acf(b1.m, lag.max=50)
plot(results.df$iteration, results.df$b2.m,
     xlab = "iteration",
     ylab = "b2.m (m on class 2 vs class 4)")
lines(results.df$iteration, results.df$b2.m)
hist (results.df$b2.m, xlab = "b2.m (m on class 2 vs class 4)")
acf(b2.m, lag.max=50)
plot(results.df$iteration, results.df$b3.m,
     xlab = "iteration",
     ylab = "b3.m (m on class 3 vs class 4")
lines(results.df$iteration, results.df$b3.m)
hist (results.df$b3.m, xlab = "b3.m (m on class 3 vs class 4)")
acf(b3.m, lag.max=50)
plot(results.df$iteration, results.df$x1.m0.y0,
     xlab = "iteration",
     ylab = "cell size x=1 m=0 y=0")
lines(results.df$iteration, results.df$x1.m0.y0)
hist (results.df$x1.m0.y0, xlab = "cell size x=1 m=0 y=0")
plot(results.df$iteration, results.df$x1.m1.y0,
     xlab = "iteration",
     ylab = "cell size x=1 m=1 y=0")
lines(results.df$iteration, results.df$x1.m1.y0)
hist (results.df$x1.m1.y0, xlab = "cell size x=1 m=1 y=0")
plot(results.df$iteration, results.df$x2.m0.y0,
     xlab = "iteration",
     ylab = "cell size x=2 m=0 y=0")
lines(results.df$iteration, results.df$x2.m0.y0)
hist (results.df$x2.m0.y0, xlab = "cell size x=2 m=0 y=0")
plot(results.df$iteration, results.df$x2.m1.y0,
     xlab = "iteration",
     ylab = "cell size x=2 m=1 y=0")
lines(results.df$iteration, results.df$x2.m1.y0)
hist (results.df$x2.m1.y0, xlab = "cell size x=2 m=1 y=0")
plot(results.df$iteration, results.df$x3.m0.y0,
     xlab = "iteration",
     ylab = "cell size x=3 m=0 y=0")
lines(results.df$iteration, results.df$x3.m0.y0)
hist (results.df$x3.m0.y0, xlab = "cell size x=3 m=0 y=0")
plot(results.df$iteration, results.df$x3.m1.y0,
     xlab = "iteration",
     ylab = "cell size x=3 m=1 y=0")
lines(results.df$iteration, results.df$x3.m1.y0)
hist (results.df$x3.m1.y0, xlab = "cell size x=3 m=1 y=0")
plot(results.df$iteration, results.df$x4.m0.y0,
     xlab = "iteration",
     ylab = "cell size x=4 m=0 y=0")
lines(results.df$iteration, results.df$x4.m0.y0)
hist (results.df$x4.m0.y0, xlab = "cell size x=4 m=0 y=0")
plot(results.df$iteration, results.df$x4.m1.y0,
     xlab = "iteration",
     ylab = "cell size x=4 m=1 y=0")
lines(results.df$iteration, results.df$x4.m1.y0)
hist (results.df$x4.m1.y0, xlab = "cell size x=4 m=1 y=0")
plot(results.df$iteration, results.df$x1.m0.y1,
     xlab = "iteration",
     ylab = "cell size x=1 m=0 y=1")
lines(results.df$iteration, results.df$x1.m0.y1)
hist (results.df$x1.m0.y1, xlab = "cell size x=1 m=0 y=1")
plot(results.df$iteration, results.df$x1.m1.y1,
     xlab = "iteration",
     ylab = "cell size x=1 m=1 y=1")
lines(results.df$iteration, results.df$x1.m1.y1)
hist (results.df$x1.m1.y1, xlab = "cell size x=1 m=1 y=1")
plot(results.df$iteration, results.df$x2.m0.y1,
     xlab = "iteration",
     ylab = "cell size x=2 m=0 y=1")
lines(results.df$iteration, results.df$x2.m0.y1)
hist (results.df$x2.m0.y1, xlab = "cell size x=2 m=0 y=1")
plot(results.df$iteration, results.df$x2.m1.y1,
     xlab = "iteration",
     ylab = "cell size x=2 m=1 y=1")
lines(results.df$iteration, results.df$x2.m1.y1)
hist (results.df$x2.m1.y1, xlab = "cell size x=2 m=1 y=1")
plot(results.df$iteration, results.df$x3.m0.y1,
     xlab = "iteration",
     ylab = "cell size x=3 m=0 y=1")
lines(results.df$iteration, results.df$x3.m0.y1)
hist (results.df$x3.m0.y1, xlab = "cell size x=3 m=0 y=1")
plot(results.df$iteration, results.df$x3.m1.y1,
     xlab = "iteration",
     ylab = "cell size x=3 m=1 y=1")
lines(results.df$iteration, results.df$x3.m1.y1)
hist (results.df$x3.m1.y1, xlab = "cell size x=3 m=1 y=1")
plot(results.df$iteration, results.df$x4.m0.y1,
     xlab = "iteration",
     ylab = "cell size x=4 m=0 y=1")
lines(results.df$iteration, results.df$x4.m0.y1)
hist (results.df$x4.m0.y1, xlab = "cell size x=4 m=0 y=1")
plot(results.df$iteration, results.df$x4.m1.y1,
     xlab = "iteration",
     ylab = "cell size x=4 m=1 y=1")
lines(results.df$iteration, results.df$x4.m1.y1)
hist (results.df$x4.m1.y1, xlab = "cell size x=4 m=1 y=1")
dev.off()

#analysis

#combine imputed class membership stored in "imp.n" with original data
#prepare 1 mplus .dat file for each imputed dataset to run mediation model (80 .dat files should be created)
for(l in 1:imp.n) {
  imp.subset <- cbind(data.original,subset(imp, imp[,1]==l))
  prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))}

#create "imp.txt" file for mplus to call imputed datasets
imp.txt <- matrix(NA,imp.n,1)
for(l in 1:imp.n) {
  imp.txt[l,1] <- paste0("imp_", l, ".dat")}
write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#run mediation model using updated PCD
runModels("2i pcd mediation.inp")

#read in parameters from updated PCD mediation model 
est.pcd <- readModels("2i pcd mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "upcd.all"
pcd.all[i,(match("tot.1v4",colnames(pcd.all))):(match("pde.4v1",colnames(pcd.all)))] <- est.pcd[100:135,"est"] # TOT_1V4 to PDE_4V1
pcd.all[i,(match("tot.1v4.se",colnames(pcd.all))):(match("pde.4v1.se",colnames(pcd.all)))] <- est.pcd[100:135,"se"] # TOT_1V4 to PDE_4V1
pcd.all[i,(match("p1",colnames(pcd.all))):(match("p4",colnames(pcd.all)))] <- est.pcd[72:75,"est"]  # P_X1 to P_X4

#flags for potential issues: "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp"

#flag if fixed threshold in unconditional model which results in missing value in tech3 covariance matrix of parameters
pcd.all[i,"fixed.th"] <- 0 #set flag to 0 to start
for(k in 1:(indicators*latent.classes+(latent.classes-1))*(indicators*latent.classes+(latent.classes-1))) { #size of tech3 matrix is 23x23 (529) 
  if (tech3[k]==999) pcd.all[i,"fixed.th"] <- 1}

#when perturbing within-class thresholds, those with a large standard error can go out of bounds (e.g., corresponding to a probability that is not between 0 and 100%)
#we created a flag for this in "results" so here we will move this into "pcd.all" and record whether this was the case in any iteration
pcd.all[i,"large.th"] <- 0 #set flag to 0 to start
results[1,"large.th"]<-0 
for(k in 1:(cycles*imp.n+burnin)) {
  if (results[k,"large.th"]==1) pcd.all[i,"large.th"] <- 1}

#we also want a flag so that we know the largest standard error for a threshold in the unconditional model (and the threshold this SE corresponds to)
#this flag captures the largest SE
pcd.all[i,"largest.th.se"] <- max(coef.se.original[,"se"])
#this flag records the threshold that corresponds to largest SE
pcd.all[i,"largest.th"] <- coef.se.original[which.max(coef.se.original[,"se"]),"est"]

#we created a flag for number of zero cells in the crosstabs for classes by mediator by outcome in "results"
#we will move this flag to "pcd.all" and record the largest number of zero cells in any iteration
results[1,"zero.cell"]<-0 #change row 1 in "results" to 0 as latent classes not imputed in first iteration
pcd.all[i,"zero.cell"] <- max(results[,"zero.cell"])

#we will also flag if there was a zero cell (in the crosstabs for classes by mediator by outcome) in any of the iterations that we saved as one of the 80 imputed datasets
zero.cell.imp <- matrix(NA,1,imp.n)
for(l in 1:imp.n) {
  zero.cell.imp[1,l] <- results[cycles*l+burnin,"zero.cell"]}

pcd.all[i,"zero.cell.imp"] <- 0 #set flag to 0 to start
for(k in 1:imp.n) {
  if (zero.cell.imp[1,k]>0) pcd.all[1,"zero.cell.imp"] <- 1}

cat(i,"\n")} #end of simulations loop

#############################################################

#saving results for comparisons of interest (e.g. with low class as the reference group) and exporting out mediation effects for each method

#onestep
onestep <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes+1)
colnames(onestep) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                       "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4","largest.th.se")

#move over class probabilities
onestep[,"p1"]<-onestep.all[,"p1"]
onestep[,"p2"]<-onestep.all[,"p2"]
onestep[,"p3"]<-onestep.all[,"p3"]
onestep[,"p4"]<-onestep.all[,"p4"]
onestep[,"largest.th.se"]<-onestep.all[,"largest.th.se"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory parameters (auc1-auc4))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc3"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc3"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc2"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc2"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#bch
bch <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(bch) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                   "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
bch[,"p1"]<-bch.all[,"p1"]
bch[,"p2"]<-bch.all[,"p2"]
bch[,"p3"]<-bch.all[,"p3"]
bch[,"p4"]<-bch.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}


#modal
modal <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(modal) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                     "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
modal[,"p1"]<-modal.all[,"p1"]
modal[,"p2"]<-modal.all[,"p2"]
modal[,"p3"]<-modal.all[,"p3"]
modal[,"p4"]<-modal.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#npcd
npcd <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(npcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                    "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
npcd[,"p1"]<-npcd.all[,"p1"]
npcd[,"p2"]<-npcd.all[,"p2"]
npcd[,"p3"]<-npcd.all[,"p3"]
npcd[,"p4"]<-npcd.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#incpcd
incpcd <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(incpcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                      "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
incpcd[,"p1"]<-incpcd.all[,"p1"]
incpcd[,"p2"]<-incpcd.all[,"p2"]
incpcd[,"p3"]<-incpcd.all[,"p3"]
incpcd[,"p4"]<-incpcd.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low) 
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc3"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc3"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc2"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc2"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#pcd
pcd <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes+6)
colnames(pcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                   "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4",
                   "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp")

#move over class probabilities
pcd[,"p1"]<-pcd.all[,"p1"]
pcd[,"p2"]<-pcd.all[,"p2"]
pcd[,"p3"]<-pcd.all[,"p3"]
pcd[,"p4"]<-pcd.all[,"p4"]
pcd[,"fixed.th"]<-pcd.all[,"fixed.th"]
pcd[,"large.th"]<-pcd.all[,"large.th"]
pcd[,"largest.th"]<-pcd.all[,"largest.th"]
pcd[,"largest.th.se"]<-pcd.all[,"largest.th.se"]
pcd[,"zero.cell"]<-pcd.all[,"zero.cell"]
pcd[,"zero.cell.imp"]<-pcd.all[,"zero.cell.imp"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#put results in a table and export to a results folder
setwd(results.files)

write.table(onestep, file="table.onestep.txt", sep = ",")
write.table(bch, file="table.bch.txt", sep = ",")
write.table(modal, file="table.modal.txt", sep = ",")
write.table(npcd, file="table.npcd.txt", sep = ",")
write.table(incpcd, file="table.incpcd.txt", sep = ",")
write.table(pcd, file="table.pcd.txt", sep = ",")

write.table(entropy, file="table.entropy.txt", sep = ",")
write.table(auc, file="table.auc.txt", sep = ",")

################################################end of script#####################################################################

############### MEDIUM ENTROPY ##########################

############### DEFINE DIRECTORIES ##########################

#Folder where simulation data are saved:
sim.data <- "file location/all sims/data/medium entropy"
#Folder where Mplus input files (uncond, onestep, bch, modal and inc) are saved:
inp.files <- "file location/all sims"
#*Folder where non-inclusive PCD files are saved:
npcd.files <- "file location/all sims"
#*Folder where inclusive PCD files are saved:
incpcd.files <- "file location/all sims"
#*Folder where updated PCD files are saved:
upcd.files <- "file location/all sims"
#*Folder where results are saved:
results.files <- "file location/all sims/results/medium entropy"

#script for simulated datasets with medium entropy (based on ALSPAC)
#############################################################
#X = 4 class latent nominal exposure (conduct trajectories: early onset persistent, adolescent onset, childhood limited and low)
latent.classes <- 4
#M = binary mediator (peer deviance)
#Y = binary outcome (problematic alcohol use)
#U1-U5 = binary latent class indicators (conduct problems from age 4 to 13 years)
indicators <- 5

############### PREPARE MATRICES ##########################

#here we are running 500 simulated datasets
sims <- 500
#simulated data has a sample size of 5000
sample <- 5000
#we wish to compare each latent class with all the others 
class.comp <- 4*3
#we wish to report 3 mediation effects (total, indirect - tie, and direct - pde) and their standard errors
mediation.effects <- 3*2

#create matrices to store the mediation results for each method
#number of rows = number of simulated datasets (here we use 500)
#number of columns is based on number of class comparisons (12), number of mediation effects (total, tie and pde) and SEs (6) and 
#number of extra parameters (e.g. area under the trajectory, class probabilities, or additional checks)

onestep.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes*2+1))

#labels represent the mediation effects (total, tie and pde) for each class comparison, area under the trajectory for each class, and class probabilities
colnames(onestep.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                           "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                           "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                           "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                           "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se", 
                           "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se", 
                           "auc1","auc2","auc3","auc4","largest.th.se","p1","p2","p3","p4")

bch.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(bch.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                       "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                       "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                       "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                       "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                       "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                       "p1","p2","p3","p4")

modal.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(modal.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                         "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                         "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                         "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                         "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                         "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                         "p1","p2","p3","p4")

npcd.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(npcd.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                        "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                        "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                        "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                        "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                        "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                        "p1","p2","p3","p4")

incpcd.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(incpcd.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                          "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                          "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                          "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                          "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                          "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                          "p1","p2","p3","p4")

pcd.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes+6))

#labels represent the mediation effects (total, tie and pde) for each class comparison, 6 flags for potential issues and class probabilities 
colnames(pcd.all) <-  c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                        "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                        "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                        "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                        "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                        "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se", "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp",
                        "p1","p2","p3","p4")

#create matrix to store the entropy of the unconditional model
entropy <- matrix(NA,sims,1)

#create matrix to store the area under the trajectory for each class in the unconditional model (this allows us to know the order of the classes)
auc <- matrix(NA,sims,latent.classes)
colnames(auc) <- c("auc1","auc2","auc3","auc4")

#create matrix to store the area under the trajectory for each class in the inclusive latent class model (this allows us to know the order of the classes)
auc.inc <- matrix(NA,sims,latent.classes)
colnames(auc.inc) <- c("auc1","auc2","auc3","auc4")

############### READ IN THE DATA ##########################

#run loop over 500 simulated datasets 
for(i in 1:sims) {
  
  #set working directory using file paths saved at the start
  setwd(sim.data)
  
  #read in Mplus .dat file with simulated data (sim1.dat)
  data.original <- read.table(file=paste0("sim", i, ".dat"), sep="", header=FALSE, col.names = c("y", #outcome
                                                                                                 "u1","u2","u3","u4","u5", #latent class indicators
                                                                                                 "m", #mediator
                                                                                                 "c" #modal class assignment for the exposure
  ))
  
  #replicate the original data and add empty columns to store imputed latent class membership for 4 classes (once generated) 
  #and interactions between each latent class and the mediator
  all.data <- data.original
  all.data$x1 <- NA #empty column to add dummy code for membership in latent class 1
  all.data$x2 <- NA #empty column to add dummy code for membership in latent class 2
  all.data$x3 <- NA #empty column to add dummy code for membership in latent class 3
  all.data$x4 <- NA #empty column to add dummy code for membership in latent class 4
  all.data$int1 <- NA  #empty column to add interaction between latent class 1 and mediator
  all.data$int2 <- NA  #empty column to add interaction between latent class 2 and mediator
  all.data$int3 <- NA  #empty column to add interaction between latent class 3 and mediator
  all.data$int4 <- NA  #empty column to add interaction between latent class 4 and mediator
  
  ############### UNCONDITIONAL LATENT CLASS MODEL ##########################
  
  #set working directory using file paths saved at the start
  setwd(inp.files)
  
  #this is simply renaming the Mplus .dat file and saving in the folder with the Mplus input files
  #this step is important when running many simulations, but not necessary otherwise
  #sim.dat will be written over each time a new simulated dataset is analysed
  prepareMplusData(data.original,"sim.dat")
  
  #run the unconditional latent class model
  runModels("2b uncond latent class.inp")
  
  #read in the Mplus output file and name it "model_output" to use later
  model_output <- readModels("2b uncond latent class.out")
  
  #save parameters from unconditional latent class model (within-class thresholds for latent class indicators and class intercepts) and their SEs
  coef.se.original <- model_output$parameters$`unstandardized`[1:(indicators*latent.classes+(latent.classes-1)),3:4] 
  #labels represent five within-class thresholds across 4 classes and 3 class intercepts)
  rownames(coef.se.original) <- c("th1.c1","th2.c1","th3.c1","th4.c1","th5.c1", #within-class thresholds for class 1
                                  "th1.c2","th2.c2","th3.c2","th4.c2","th5.c2", #within-class thresholds for class 2
                                  "th1.c3","th2.c3","th3.c3","th4.c3","th5.c3", #within-class thresholds for class 3
                                  "th1.c4","th2.c4","th3.c4","th4.c4","th5.c4", #within-class thresholds for class 4
                                  "int.c1","int.c2","int.c3") #latent class intercepts
  
  #save model entropy into the matrix we created at the start of the script
  entropy[i,1] <- model_output$summaries$`Entropy`
  
  #create a matrix of class probabilities
  p.original <- matrix(NA,latent.classes,1)
  #these can be calculated using the 3 class intercepts we saved in "coef.se.original"
  p.original[1,1] <- exp(coef.se.original["int.c1","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
  p.original[2,1] <- exp(coef.se.original["int.c2","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
  p.original[3,1] <- exp(coef.se.original["int.c3","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))  
  p.original[4,1] <- 1/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
  
  #derive area under the trajectory parameters so we know the order of the latent classes in the unconditional model
  #these can be calculated using the within-class thresholds we saved in "coef.se.original"
  auc[i,"auc1"] <- exp(-1*coef.se.original["th1.c1","est"])/(1+exp(-1*coef.se.original["th1.c1","est"]))+2*exp(-1*coef.se.original["th2.c1","est"])/(1+exp(-1*coef.se.original["th2.c1","est"]))+3*exp(-1*coef.se.original["th3.c1","est"])/(1+exp(-1*coef.se.original["th3.c1","est"]))+4*exp(-1*coef.se.original["th4.c1","est"])/(1+exp(-1*coef.se.original["th4.c1","est"]))+5*exp(-1*coef.se.original["th5.c1","est"])/(1+exp(-1*coef.se.original["th5.c1","est"]))
  auc[i,"auc2"] <- exp(-1*coef.se.original["th1.c2","est"])/(1+exp(-1*coef.se.original["th1.c2","est"]))+2*exp(-1*coef.se.original["th2.c2","est"])/(1+exp(-1*coef.se.original["th2.c2","est"]))+3*exp(-1*coef.se.original["th3.c2","est"])/(1+exp(-1*coef.se.original["th3.c2","est"]))+4*exp(-1*coef.se.original["th4.c2","est"])/(1+exp(-1*coef.se.original["th4.c2","est"]))+5*exp(-1*coef.se.original["th5.c2","est"])/(1+exp(-1*coef.se.original["th5.c2","est"]))
  auc[i,"auc3"] <- exp(-1*coef.se.original["th1.c3","est"])/(1+exp(-1*coef.se.original["th1.c3","est"]))+2*exp(-1*coef.se.original["th2.c3","est"])/(1+exp(-1*coef.se.original["th2.c3","est"]))+3*exp(-1*coef.se.original["th3.c3","est"])/(1+exp(-1*coef.se.original["th3.c3","est"]))+4*exp(-1*coef.se.original["th4.c3","est"])/(1+exp(-1*coef.se.original["th4.c3","est"]))+5*exp(-1*coef.se.original["th5.c3","est"])/(1+exp(-1*coef.se.original["th5.c3","est"]))
  auc[i,"auc4"] <- exp(-1*coef.se.original["th1.c4","est"])/(1+exp(-1*coef.se.original["th1.c4","est"]))+2*exp(-1*coef.se.original["th2.c4","est"])/(1+exp(-1*coef.se.original["th2.c4","est"]))+3*exp(-1*coef.se.original["th3.c4","est"])/(1+exp(-1*coef.se.original["th3.c4","est"]))+4*exp(-1*coef.se.original["th4.c4","est"])/(1+exp(-1*coef.se.original["th4.c4","est"]))+5*exp(-1*coef.se.original["th5.c4","est"])/(1+exp(-1*coef.se.original["th5.c4","est"]))
  
  #import tech3 as a covariance matrix to allow us to use the covariance between parameters for perturbing in updated PCD
  #############################################################
  #NB the order of parameters in this covariance matrix is different to order in "coef.se.original" (see tech1) so this will need to be addressed before perturbing parameters
  #this is due to no missing data in the indicators which affects the numbering of parameters in tech1
  #############################################################
  #this outputs paramCov as a matrix object in R 
  tech3 <- model_output$tech3$paramCov
  #this creates a full, symmetrical covariance matrix
  upperTriangle(tech3) <- lowerTriangle(tech3, byrow=TRUE)
  cov <- tech3
  
  #address any fixed parameters: within class thresholds that have been fixed at 15 or -15 (representing a probability of 0 or 100%) do not have (co)variances
  #the code below means that fixed parameters will still get perturbed in updated PCD but only a very small amount
  #replace missing (999) in covariance matrix with 0 to represent no covariance for the fixed parameters
  cov[cov==999] <- 0
  #change the variance for fixed parameters to be very small (0.000000001)
  for(j in 1:(indicators*latent.classes+(latent.classes-1))) {
    if (cov[j,j]==0) cov[j,j]<- 0.000000001}
  
  ###################### ONESTEP MODEL ##########################
  
  #run the onestep latent class mediation model
  runModels("2c onestep mediation.inp")
  
  #read in parameters from one-step latent class model    
  est.onestep <- readModels("2c onestep mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects, area under the trajectory for each class, and class probabilities) and their SEs in the matrix created at the start "onestep.all"
  onestep.all[i, (match("tot.1v4",colnames(onestep.all))):(match("pde.4v1",colnames(onestep.all)))] <- est.onestep[116:151,"est"] # TOT_1V4 to PDE_4V1
  onestep.all[i, (match("tot.1v4.se",colnames(onestep.all))):(match("pde.4v1.se",colnames(onestep.all)))] <- est.onestep[116:151,"se"] # TOT_1V4 to PDE_4V1
  onestep.all[i, (match("auc1",colnames(onestep.all))):(match("auc4",colnames(onestep.all)))] <- est.onestep[152:155,"est"] # AUC
  onestep.all[i, (match("p1",colnames(onestep.all))):(match("p4",colnames(onestep.all)))] <- est.onestep[88:91,"est"]  # P_X1 to P_X4
  
  #save standard errors for thresholds from onestep model
  se <- est.onestep[c(3:7,17:21,31:35,45:49),4] 
  
  #we want a flag so that we know the largest standard error for a threshold in the onestep model
  #this flag captures the largest SE
  onestep.all[i,"largest.th.se"] <- max(se)
  
  ######################## BCH ##########################
  
  #bch can only be used with multiple latent classes if they are combined into one class
  #use one 8 class model for X and M and set up mediation model using nom-nom-cat approach (because bch needs to include an auxiliary variable)
  
  #read in bch weights and modal class assignment which were exported out of the unconditional latent class model in "bch.txt"
  data.bch <- read.table("bch.txt", sep="", na.strings="*", header=FALSE)
  #select the weights for each class and the modal class assignment
  data.bch <- data.bch[,c(6:9,14)]
  #labels representing one weight for each class and the modal class assignment
  colnames(data.bch) <- c("bch1","bch2","bch3","bch4","modal")
  #combine the bch weights with the mediator and outcome from the original data
  data.bch <- cbind(data.original[,c("y","m")], data.bch)
  
  #create weights for 8 class model by multiplying the bch weight for each class with the observed data on the mediator
  data.bch$bch11 <- data.bch$bch1*data.bch$m
  data.bch$bch12 <- data.bch$bch1*(1-data.bch$m)
  data.bch$bch21 <- data.bch$bch2*data.bch$m
  data.bch$bch22 <- data.bch$bch2*(1-data.bch$m)
  data.bch$bch31 <- data.bch$bch3*data.bch$m
  data.bch$bch32 <- data.bch$bch3*(1-data.bch$m)
  data.bch$bch41 <- data.bch$bch4*data.bch$m
  data.bch$bch42 <- data.bch$bch4*(1-data.bch$m)
  
  #prepare Mplus .dat file using "data.bch"  
  prepareMplusData(data.bch,"bch.dat")
  #run mediation model using bch method
  runModels("2d bch mediation.inp")
  
  #read in parameters from bch latent class mediation model   
  est.bch <- readModels("2d bch mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "bch.all"
  bch.all[i,(match("tot.1v4",colnames(bch.all))):(match("pde.4v1",colnames(bch.all)))] <- est.bch[68:103,"est"] # TOT_1V4 to PDE_4V1
  bch.all[i,(match("tot.1v4.se",colnames(bch.all))):(match("pde.4v1.se",colnames(bch.all)))] <- est.bch[68:103,"se"] # TOT_1V4 to PDE_4V1
  bch.all[i,(match("p1",colnames(bch.all))):(match("p4",colnames(bch.all)))] <- est.bch[40:43,"est"]  # P_X1 to P_X4
  
  ######################## MODAL CLASS ASSIGNMENT ##########################
  
  #run mediation model using modal class assignment
  runModels("2e modal mediation.inp")
  
  #read in parameters from modal class mediation model  
  est.modal <- readModels("2e modal mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "modal.all"
  modal.all[i,(match("tot.1v4",colnames(modal.all))):(match("pde.4v1",colnames(modal.all)))] <- est.modal[100:135,"est"] # TOT_1V4 to PDE_4V1
  modal.all[i,(match("tot.1v4.se",colnames(bch.all))):(match("pde.4v1.se",colnames(bch.all)))] <- est.modal[100:135,"se"] # TOT_1V4 to PDE_4V1
  modal.all[i,(match("p1",colnames(bch.all))):(match("p4",colnames(bch.all)))] <- est.modal[72:75,"est"]  # P_X1 to P_X4
  
  ######################## NON-INCLUSIVE PCD ##########################
  
  #set seed for non-inclusive PCD
  #this will set the seed for each session in R (not each time run loops)
  #to make sure you get the same results each time, run the loop only once within each R session
  set.seed(80)
  
  #calculate P(X=x|U) - this will give us the probability of class membership (cprobs) from unconditional latent class model
  
  #for each class, multiply the individual data (responses to 5 binary indicators: U1 to U5) with within-class thresholds
  #number of latent class indicators
  #class 1
  #this creates a matrix with within-class thresholds for class 1 repeated for every individual in the dataset (e.g. repeated 5000 times)
  theta <- matrix(rep(coef.se.original[(match("th1.c1",rownames(coef.se.original))):(match("th5.c1",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
  #if latent class indicator is present (-1)*threshold is used, if latent class indicator is absent (1)*threshold is used
  P1<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  #class 2
  #this creates a matrix with within-class thresholds for class 2 repeated for every individual in the dataset (e.g. repeated 5000 times)
  theta <- matrix(rep(coef.se.original[(match("th1.c2",rownames(coef.se.original))):(match("th5.c2",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
  P2<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  #class 3
  #this creates a matrix with within-class thresholds for class 3 repeated for every individual in the dataset (e.g. repeated 5000 times)
  theta <- matrix(rep(coef.se.original[(match("th1.c3",rownames(coef.se.original))):(match("th5.c3",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
  P3<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  #class 4
  #this creates a matrix with within-class thresholds for class 4 repeated for every individual in the dataset (e.g. repeated 5000 times)
  theta <- matrix(rep(coef.se.original[(match("th1.c4",rownames(coef.se.original))):(match("th5.c4",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
  P4<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  
  #multiply the 5 probabilities (for each latent class indicator) in each row of P* then multiply by class probabilities that were saved earlier in "p.original"  
  N1<-apply(P1,1,prod)*p.original[1,1]
  N2<-apply(P2,1,prod)*p.original[2,1]
  N3<-apply(P3,1,prod)*p.original[3,1]
  N4<-apply(P4,1,prod)*p.original[4,1]    
  
  #this gives us the probability of class membership for each person in the dataset 
  P<-cbind(N1/(N1+N2+N3+N4),N2/(N1+N2+N3+N4),N3/(N1+N2+N3+N4),N4/(N1+N2+N3+N4))
  
  #we could have simply read these in from "bch.txt" (as we did with the bch weights) as they were exported from the unconditional latent class model 
  colnames(P) <-  c("cprob1","cprob2","cprob3","cprob4")
  
  #we will now impute class membership 40 times for each person using their probability of class membership
  #we will create 40 imputed datasets (chosen to keep Monte Carlo error at less than 10% of standard error for parameters from regression model for Y)
  imp.n <- 40
  #create a matrix to store imputed class membership for each person (sample*imp.n=5000*40)
  imp <- matrix(NA,sample*imp.n,latent.classes+1)
  colnames(imp) <-  c("imp","x1","x2","x3","x4")
  
  #first column is simply an indicator for imputation number (range from 1 to 40)
  for(h in 0:(imp.n-1)) {
    imp [c(sample*h+1:sample),1]<-h+1}
  
  #for each of 40 imputations we will create a matrix to store class membership, 
  #then we will use the probabilities of class membership "P" to randomly assign each individual to a class (X = 1,..k).
  #and we will add their imputed class membership into the matrix "imp" we created earlier
  for(j in 1:imp.n) {
    x <- matrix(NA,sample,latent.classes)    
    for(k in 1:sample) {x[k,1:latent.classes] <- rmultinom(1,1,P[k,])}
    for(l in 1:imp.n) {
      if (j==l) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x}
  } 
  
  #analysis
  #can use a separate folder for non-inclusive PCD because of large number of files created (e.g., 40 imputed datasets)
  setwd(npcd.files) 
  
  #combine imputed class membership stored in "imp.n" with original data
  #prepare 1 mplus .dat file for each imputed dataset to run mediation model (40 .dat files should be created)
  for(l in 1:imp.n) {
    imp.subset <- cbind(data.original,subset(imp, imp[,1]==l))
    prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))
  } 
  
  #create "imp.txt" file for mplus to call imputed datasets
  imp.txt <- matrix(NA,imp.n,1)
  for(l in 1:imp.n) {
    imp.txt[l,1] <- paste0("imp_", l, ".dat")
  } 
  write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #run mediation model using non-inclusive PCD
  runModels("2f npcd mediation.inp")
  
  #read in parameters from npcd mediation model 
  est.npcd <- readModels("2f npcd mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "npcd.all"
  npcd.all[i,(match("tot.1v4",colnames(npcd.all))):(match("pde.4v1",colnames(npcd.all)))] <- est.npcd[100:135,"est"] # TOT_1V4 to PDE_4V1
  npcd.all[i,(match("tot.1v4.se",colnames(npcd.all))):(match("pde.4v1.se",colnames(npcd.all)))] <- est.npcd[100:135,"se"] # TOT_1V4 to PDE_4V1
  npcd.all[i,(match("p1",colnames(npcd.all))):(match("p4",colnames(npcd.all)))] <- est.npcd[72:75,"est"]  # P_X1 to P_X4
  
  ######################## INCLUSIVE PCD ##########################
  
  #set seed for inclusive PCD
  set.seed(81)
  
  #set working directory using file paths saved at the start
  setwd(inp.files)
  
  #run inclusive latent class model and export out cprobs
  runModels("2g inc latent class.inp")
  
  #read in auc parameters from inclusive latent class model
  est.inc <- readModels("2g inc latent class.out", what="parameters")$parameters$`unstandardized`
  auc.inc[i,(match("auc1",colnames(auc.inc))):(match("auc4",colnames(auc.inc)))] <- est.inc[30:33,"est"]  # AUC1 to AUC4
  
  #read in class probabilities from "inc.txt" (as we did with the bch weights exported from the unconditional latent class model) 
  cprobs <- read.table("inc.txt", sep="", na.strings="*", header=FALSE)
  #we could also extract modal classes from here (v12) to use the inclusive modal approach
  P <- cprobs[,c(8:11)]
  colnames(P) <- c("cprob1","cprob2","cprob3","cprob4")
  
  #we will now impute class membership 40 times for each person using their probability of class membership
  #create a matrix to store imputed class membership for each person (sample*imp.n=5000*40)
  imp <- matrix(NA,sample*imp.n,latent.classes+1)
  colnames(imp) <-  c("imp","x1","x2","x3","x4")
  
  #first column is simply an indicator for imputation number (range from 1 to 40)
  for(h in 0:(imp.n-1)) {
    imp [c(sample*h+1:sample),1]<-h+1}
  
  #for each of 40 imputations we will create a matrix to store class membership, 
  #then we will use the probabilities of class membership "P" to randomly assign each individual to a class (X = 1,..k).
  #and we will add their imputed class membership into the matrix "imp" we created earlier
  for(j in 1:imp.n) {
    x <- matrix(NA,sample,latent.classes)    
    for(k in 1:sample) {x[k,1:latent.classes] <- rmultinom(1,1,P[k,])}
    for(l in 1:imp.n) {
      if (j==l) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x}
  } 
  
  #analysis
  #can use a separate folder for inclusive PCD because of large number of files created (e.g., 40 imputed datasets)
  setwd(incpcd.files) 
  
  #combine imputed class membership stored in "imp.n" with original data
  #prepare 1 mplus .dat file for each imputed dataset to run mediation model (40 .dat files should be created)
  for(l in 1:imp.n) {
    imp.subset <- cbind(data.original,subset(imp, imp[,1]==l))
    prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))
  } 
  
  #create "imp.txt" file for mplus to call imputed datasets
  imp.txt <- matrix(NA,imp.n,1)
  for(l in 1:imp.n) {
    imp.txt[l,1] <- paste0("imp_", l, ".dat")
  } 
  write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #run mediation model using inclusive PCD
  runModels("2h incpcd mediation.inp")
  
  #read in parameters from incpcd mediation model 
  est.incpcd <- readModels("2h incpcd mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "incpcd.all"
  incpcd.all[i,(match("tot.1v4",colnames(incpcd.all))):(match("pde.4v1",colnames(incpcd.all)))] <- est.incpcd[100:135,"est"] # TOT_1V4 to PDE_4V1
  incpcd.all[i,(match("tot.1v4.se",colnames(incpcd.all))):(match("pde.4v1.se",colnames(incpcd.all)))] <- est.incpcd[100:135,"se"] # TOT_1V4 to PDE_4V1
  incpcd.all[i,(match("p1",colnames(incpcd.all))):(match("p4",colnames(incpcd.all)))] <- est.incpcd[72:75,"est"]  # P_X1 to P_X4
  
  #remove any objects from memory that will not be needed for the rest of the script
  rm(model_output, est.onestep, data.bch, est.bch, est.modal, est.npcd, est.incpcd, cprobs, P)
  
  ######################## UPDATED PCD ##########################
  
  #set seed for updated PCD
  set.seed(82)
  
  #STEP 1: 
  #initialise beta by fitting the logistic regression model for Y with beta for classes set to zero 
  #initialise alpha by fitting the logistic regression model for M with alpha for classes set to zero
  #perturb the estimates from these models with Gaussian noise of mean zero and values from the variance-covariance matrix of these parameter estimates
  
  #firth logistic regression is used to address problems with perfect prediction when bringing in exposure-mediator interactions later in the script
  #regress outcome y on mediator m
  model.y <- logistf(y~m, data=all.data)
  #regression model for mediator m 
  model.m <- logistf(m~1, data=all.data)
  
  #perturbing the coefficients once around coefficients in models using variance-covariance matrix 
  #initially use zero for exposure coefficients (and exposure-mediator interactions in outcome model)
  #once latent class exposure has been imputed in subsequent runs, these will become coefficients from models
  beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y)),0,0,0,0,0,0)
  beta.m<-c(rmnorm(1,coef(model.m),vcov(model.m)),0,0,0)
  
  #we will create 80 imputed datasets for class membership (chosen to keep Monte Carlo error at less than 10% of standard error for parameters in regression model for Y)
  imp.n <- 80
  #we will allow 20 iterations between saving out imputed class membership
  cycles <- 20
  #we will allow a burn in of 100 iterations before starting to save out imputed class membership
  burnin <- 100
  
  #create a matrix to save the results from every iteration (all results are recorded to assess convergence later)
  #rows=iterations=(cycles*imp.n+burnin)
  #columns=iteration number, beta coeficients from regression models, cell sizes and 2 flags for issues=(6+(latent.classes-1)*3+latent.classes*4)=31
  results <- matrix(NA,cycles*imp.n+burnin,6+(latent.classes-1)*3+latent.classes*4) 
  
  colnames(results) <- c("iteration", #iteration number
                         "b0.y","b1.y","b2.y","b3.y","b4.y","b5.y","b6.y","b7.y", #coefficients from regression model for the outcome
                         "b0.m","b1.m","b2.m","b3.m", #coefficients from regression model for the mediator
                         "x1.m0.y0","x1.m1.y0","x2.m0.y0","x2.m1.y0","x3.m0.y0","x3.m1.y0","x4.m0.y0","x4.m1.y0","x1.m0.y1","x1.m1.y1","x2.m0.y1","x2.m1.y1","x3.m0.y1","x3.m1.y1","x4.m0.y1","x4.m1.y1", #cell sizes from crosstabs for classes by mediator by outcome
                         "zero.cell", #flag for presence of zero cells in crosstabs
                         "large.th" #flag for a within-class threshold in unconditional latent class model that was out of bounds after perturbing
  ) 
  
  #the first column is simply an indicator of iteration number (range from 1 to 1700)
  results[,1]<-1:(cycles*imp.n+burnin)
  
  #create a matrix to store cell sizes from crosstabs for classes by mediator by outcome 
  xmy<-matrix(NA,1,latent.classes*4)
  colnames(xmy) <-  c("x1m0y0","x1m1y0","x2m0y0","x2m1y0","x3m0y0","x3m1y0","x4m0y0","x4m1y0","x1m0y1","x1m1y1","x2m0y1","x2m1y1","x3m0y1","x3m1y1","x4m0y1","x4m1y1")
  #create a matrix to store presence of zero cells in this crosstabs
  zero.cell<-matrix(NA,1,1)
  #create a matrix to store presence of a within-class threshold in unconditional latent class model that was out of bounds after perturbing (e.g., not corresponding to 0 to 100% probability) 
  large.th<-matrix(NA,1,1)
  
  #create a matrix to store imputed class membership for each person (sample*imp.n=5000*80)
  imp <- matrix(NA,sample*imp.n,latent.classes+1)
  colnames(imp) <- c("imp","x1","x2","x3","x4") 
  #first column is simply an indicator for imputation number (range from 1 to 80)
  for(h in 0:(imp.n-1)) {
    imp [c(sample*h+1:sample),1]<-h+1}
  
  #we will now create a loop to repeat steps 2 and 3 below 1700 (imp.n*cycles+burnin) times and save estimates after every 20 iterations (after a 100 iteration burn in)
  for(j in 1:(imp.n*cycles+burnin)) {
    results[j,2:9]<-beta.y #save the perturbed coefficients from regression model for the outcome "beta.y"
    results[j,10:13]<-beta.m #save the perturbed coefficients from regression model for the mediator "beta.m"    
    results[j,14:29]<-xmy #save cell sizes from crosstabs for classes by mediator by outcome "xmy"
    results[j,30]<-zero.cell #save a flag for presence of zero cells in this crosstabs "zero.cell"
    results[j,31]<-large.th #save a flag for a within-class threshold in unconditional latent class model that was out of bounds after peturbing "large.th"
    
    #this saves each individual's imputed class membership 80 times (generated later in the script) into the matrix "imp" we created earlier     
    for(l in 1:imp.n) {
      if (j==(cycles*l+burnin)) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x
    } 
    
    #STEP 2a:
    #we will now perturb the parameters (within-class thresholds and class intercepts) that we saved earlier "coef.se.original" from the unconditional latent class model 
    #in order to perturb the parameters using the covariance matrix in tech3 we need to reorder parameters so they match with numbering in tech1
    #because we have no missing data, they are in a different order to what would be expected
    #when using a dataset with missing data on class indicators (as would usually be the case outside of simulated data), this reordering step is not needed
    
    #keep only the parameters (within-class thresholds for latent class indicators and class intercepts) from unconditional latent class model (i.e. drop the SEs)
    coef.mplus.orig<-coef.se.original[,"est"]
    #reorder the parameters so that they are in the same order as is used in the covariances matrix of parameters from tech3
    index<-c(1,5,9,13,17,2,6,10,14,18,3,7,11,15,19,4,8,12,16,20,21,22,23)
    coef.mplus.orig<-coef.mplus.orig[order(index)]
    #perturb these parameters based on their variance-covariance matrix (saved in "cov")
    coef.mplus<-rmnorm(1,coef.mplus.orig,cov)
    #reorder again to preserve original ordering
    index<-c(1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19,5,10,15,20,21,22,23)
    coef<-coef.mplus[order(index)]
    
    #when perturbing within-class thresholds, those with a large standard error can go out of bounds (e.g., corresponding to a probability that is not between 0 and 100%)
    #we will create a flag so we know when this is the case
    large.th[1,1]<-0    
    for(k in 1:(latent.classes*indicators)) {
      if (coef[k]>15||coef[k]<(-15)) large.th[1,1] <- 1} 
    #we will also constrain within-class thresholds to be between -15 and 15 (corresponding to a probability that is between 0 and 100%)
    coef[coef>15] <- 15
    coef[coef<(-15)] <- (-15)
    
    #turning 'coef' into a matrix so can label rows and columns to use in script below (instead of referring to numbers)
    coef<-as.matrix(coef)
    rownames(coef) <- c("th1.c1","th2.c1","th3.c1","th4.c1","th5.c1",
                        "th1.c2","th2.c2","th3.c2","th4.c2","th5.c2",
                        "th1.c3","th2.c3","th3.c3","th4.c3","th5.c3",
                        "th1.c4","th2.c4","th3.c4","th4.c4","th5.c4",
                        "int.c1","int.c2","int.c3")
    colnames(coef) <- c("est")
    
    #now calculate P(X=x|U) - without perturbing this give us the probability of class membership (cprobs) from unconditional latent class model
    
    #for each class, multiply the individual data (responses to 5 binary indicators: U1 to U5) with within-class thresholds
    #class 1
    #this creates a matrix with within-class thresholds for class 1 repeated for every individual in the dataset (e.g. repeated 5000 times)
    theta <- matrix(rep(coef[(match("th1.c1",rownames(coef.se.original))):(match("th5.c1",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
    #if latent class indicator is present (-1)*threshold is used, if latent class indicator is absent (1)*threshold is used
    P1<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
    #class 2
    theta <- matrix(rep(coef[(match("th1.c2",rownames(coef.se.original))):(match("th5.c2",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
    P2<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
    #class 3
    theta <- matrix(rep(coef[(match("th1.c3",rownames(coef.se.original))):(match("th5.c3",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
    P3<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
    #class 4
    theta <- matrix(rep(coef[(match("th1.c4",rownames(coef.se.original))):(match("th5.c4",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
    P4<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
    
    #we need to create a new matrix of class probabilities using the perturbed class intercepts
    p <- matrix(NA,latent.classes,1)
    p[1,1] <- exp(coef["int.c1","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
    p[2,1] <- exp(coef["int.c2","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
    p[3,1] <- exp(coef["int.c3","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))  
    p[4,1] <- 1/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
    
    #then multiply "P1" to "P4" by the perturbed class probabilities that we have created above "p" 
    #this gives us the probability of class membership for each person in the dataset
    #these will differ slightly to "cprobs" that can be exported from the unconditional latent class model due to perturbing
    N1<-apply(P1,1,prod)*p[1,1]
    N2<-apply(P2,1,prod)*p[2,1]
    N3<-apply(P3,1,prod)*p[3,1]
    N4<-apply(P4,1,prod)*p[4,1]    
    
    #STEP 2b - Bayes rule:
    
    #combine P(X=x|U) along with coefficients from regression models for outcome and mediator
    #this will calculate probabilities P(X=x|Y,M,U) for each class (x = 1,..k) for each individual
    #we will do this via a number a steps below
    
    #P(Y=1|M,X): calculate probability that the outcome (Y) = 1 given mediator (M) and exposure (X, latent classes) 
    #this uses the coefficients from the regression model for the outcome (intercept (b0), coef for M (b1), coefs for X (b2-b4) and coefs for XM interaction (b5-b7))
    #class 1
    num1 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m+results[j,"b2.y"]+results[j,"b5.y"]*all.data$m)
    Y11 <- num1/(1+num1)  
    #class 2
    num2 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m+results[j,"b3.y"]+results[j,"b6.y"]*all.data$m)
    Y21 <- num2/(1+num2)  
    #class 3
    num3 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m+results[j,"b4.y"]+results[j,"b7.y"]*all.data$m)
    Y31 <- num3/(1+num3)
    #class 4
    num4 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m)
    Y41 <- num4/(1+num4)
    
    #P(Y=0|M,X): calculate probability that the outcome (Y) = 0 given mediator (M) and exposure (X, latent classes) 
    Y10 <- 1/(1+num1)
    Y20 <- 1/(1+num2)
    Y30 <- 1/(1+num3)
    Y40 <- 1/(1+num4)
    
    #P(M=1|X): calculate probability that the mediator (M) = 1 given exposure (X, latent classes) 
    #this uses the coefficients from the regression model for the mediator (intercept (b0) and coefs for X (b1-b3))
    #class 1 
    num1 <- exp(results[j,"b0.m"]+results[j,"b1.m"])
    M11 <- num1/(1+num1)
    #class 2
    num2 <- exp(results[j,"b0.m"]+results[j,"b2.m"])
    M21 <- num2/(1+num2) 
    #class 3
    num3 <- exp(results[j,"b0.m"]+results[j,"b3.m"])
    M31 <- num3/(1+num3)   
    #class 4
    num4 <- exp(results[j,"b0.m"])
    M41 <- num4/(1+num4) 
    
    #P(M=0|X): calculate probability that the mediator (M) = 0 given exposure (X, latent classes)    
    M10 <- 1/(1+num1)
    M20 <- 1/(1+num2)
    M30 <- 1/(1+num3)      
    M40 <- 1/(1+num4)
    
    #P(X=x|Y=1,M=1,U)
    #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
    N111<-N1*Y11*M11
    #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
    N211<-N2*Y21*M21
    #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
    N311<-N3*Y31*M31
    #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
    N411<-N4*Y41*M41
    
    denom <- N111+N211+N311+N411
    Q111<-N111/denom
    Q211<-N211/denom
    Q311<-N311/denom
    Q411<-N411/denom 
    Q11<-cbind(Q111,Q211,Q311,Q411)
    head(Q11)
    
    #P(X=x|Y=0,M=0,U) 
    #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
    N100<-N1*Y10*M10
    #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
    N200<-N2*Y20*M30
    #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
    N300<-N3*Y30*M30
    #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
    N400<-N4*Y40*M40
    
    denom <- N100+N200+N300+N400   
    Q100<-N100/denom
    Q200<-N200/denom
    Q300<-N300/denom
    Q400<-N400/denom    
    Q00<-cbind(Q100,Q200,Q300,Q400)
    head(Q00)
    
    #P(X=x|Y=1,M=0,U)
    #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
    N101<-N1*Y11*M10
    #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
    N201<-N2*Y21*M30
    #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
    N301<-N3*Y31*M30
    #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
    N401<-N4*Y41*M40
    
    denom <- N101+N201+N301+N401     
    Q101<-N101/denom
    Q201<-N201/denom
    Q301<-N301/denom
    Q401<-N401/denom    
    Q01<-cbind(Q101,Q201,Q301,Q401)
    head(Q01)    
    
    #P(X=x|Y=0,M=1,U)
    #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
    N110<-N1*Y10*M11
    #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
    N210<-N2*Y20*M31
    #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
    N310<-N3*Y30*M31
    #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
    N410<-N4*Y40*M41
    
    denom <- N110+N210+N310+N410    
    Q110<-N110/denom
    Q210<-N210/denom
    Q310<-N310/denom
    Q410<-N410/denom   
    Q10<-cbind(Q110,Q210,Q310,Q410)
    head(Q10)   
    
    #derive the probability of class membership for each individual which takes into account the relationship between the classes, mediator and outcome
    #using each individuals observed data on the mediator and outcome
    Q<-Q11
    #probabilities for those with mediator and outcome absent
    Q[all.data$y==0 & all.data$m==0,]<-Q00[all.data$y==0 & all.data$m==0,]
    #probabilities for those with mediator absent and outcome present
    Q[all.data$y==1 & all.data$m==0,]<-Q01[all.data$y==1 & all.data$m==0,]  
    #probabilities for those with mediator present and outcome absent
    Q[all.data$y==0 & all.data$m==1,]<-Q10[all.data$y==0 & all.data$m==1,]
    #these probabilities will be used to impute class membership for each individual 
    head(Q)
    
    colnames(Q) <-  c("cprob1","cprob2","cprob3","cprob4")
    
    #STEP 3:
    
    #create a matrix to store class membership 
    x <- matrix(NA,sample,latent.classes)
    #now we will use the probabilities of class membership "Q" to randomly assign each individual to a class (X = 1,..k).
    #n=1 (number of random vectors to draw); size=1 per person
    for(k in 1:sample) {x[k,1:latent.classes] <- rmultinom(1,1,Q[k,])}
    
    #add imputed class membership to the data
    all.data[,"x1"] <- x[,1]
    all.data[,"x2"] <- x[,2]
    all.data[,"x3"] <- x[,3]
    all.data[,"x4"] <- x[,4]
    #add in exposure-mediator interactions
    all.data[,"int1"] <- all.data[,"x1"]*all.data[,"m"]
    all.data[,"int2"] <- all.data[,"x2"]*all.data[,"m"]
    all.data[,"int3"] <- all.data[,"x3"]*all.data[,"m"]
    all.data[,"int4"] <- all.data[,"x4"]*all.data[,"m"]
    
    ###################################################
    #we will now perform some checks on the cell sizes from crosstabs for classes by mediator by outcome 
    #we will add these cell sizes into "results" to make traceplots to assess convergence
    
    #create a subset of the data for Y=0 
    no.y <- subset(all.data, all.data$y==0)
    #create a subset of the data for Y=1
    yes.y <- subset(all.data, all.data$y==1)
    
    #crosstabs for x1 and m for those with y=0
    x1m.no.y<-table(no.y$x1,no.y$m)
    
    #this is to make sure matrix is 2 by 2 even when there are zero cells
    if (nrow(x1m.no.y)==1) x1m.no.y <- rbind(x1m.no.y,matrix(0,1,2))
    #crosstabs for x2 and m for those with y=0 
    x2m.no.y<-table(no.y$x2,no.y$m)
    if (nrow(x2m.no.y)==1) x2m.no.y <- rbind(x2m.no.y,matrix(0,1,2))
    #crosstabs for x3 and m for those with y=0 
    x3m.no.y<-table(no.y$x3,no.y$m)
    if (nrow(x3m.no.y)==1) x3m.no.y <- rbind(x3m.no.y,matrix(0,1,2))
    #crosstabs for x4 and m for those with y=0 
    x4m.no.y<-table(no.y$x4,no.y$m)
    if (nrow(x4m.no.y)==1) x4m.no.y <- rbind(x4m.no.y,matrix(0,1,2))
    #crosstabs for x1 and m for those with y=1 
    x1m.yes.y<-table(yes.y$x1,yes.y$m)
    if (nrow(x1m.yes.y)==1) x1m.yes.y <- rbind(x1m.yes.y,matrix(0,1,2))    
    x2m.yes.y<-table(yes.y$x2,yes.y$m)
    if (nrow(x2m.yes.y)==1) x2m.yes.y <- rbind(x2m.yes.y,matrix(0,1,2))    
    x3m.yes.y<-table(yes.y$x3,yes.y$m)
    if (nrow(x3m.yes.y)==1) x3m.yes.y <- rbind(x3m.yes.y,matrix(0,1,2))
    x4m.yes.y<-table(yes.y$x4,yes.y$m)
    if (nrow(x4m.yes.y)==1) x4m.yes.y <- rbind(x4m.yes.y,matrix(0,1,2))
    
    #use empty xmy matrix created earlier and fill in with cell sizes
    #cell size for x=1, m=0, y=0
    xmy[1,"x1m0y0"] <- x1m.no.y[2,1]
    #cell size for x=1, m=1, y=0  
    xmy[1,"x1m1y0"] <- x1m.no.y[2,2]
    #cell size for x=2, m=0, y=0  
    xmy[1,"x2m0y0"] <- x2m.no.y[2,1]
    #cell size for x=2, m=1, y=0  
    xmy[1,"x2m1y0"] <- x2m.no.y[2,2]
    #cell size for x=3, m=0, y=0  
    xmy[1,"x3m0y0"] <- x3m.no.y[2,1]
    #cell size for x=3, m=1, y=0  
    xmy[1,"x3m1y0"] <- x3m.no.y[2,2]
    #cell size for x=4, m=0, y=0  
    xmy[1,"x4m0y0"] <- x4m.no.y[2,1]
    #cell size for x=4, m=1, y=0  
    xmy[1,"x4m1y0"] <- x4m.no.y[2,2]
    #cell size for x=1, m=0, y=1  
    xmy[1,"x1m0y1"] <- x1m.yes.y[2,1]
    #cell size for x=1, m=1, y=1  
    xmy[1,"x1m1y1"] <- x1m.yes.y[2,2]
    #cell size for x=2, m=0, y=1  
    xmy[1,"x2m0y1"] <- x2m.yes.y[2,1]
    #cell size for x=2, m=1, y=1  
    xmy[1,"x2m1y1"] <- x2m.yes.y[2,2]
    #cell size for x=3, m=0, y=1  
    xmy[1,"x3m0y1"] <- x3m.yes.y[2,1]
    #cell size for x=3, m=1, y=1  
    xmy[1,"x3m1y1"] <- x3m.yes.y[2,2]
    #cell size for x=4, m=0, y=1  
    xmy[1,"x4m0y1"] <- x4m.yes.y[2,1]
    #cell size for x=4, m=1, y=1  
    xmy[1,"x4m1y1"] <- x4m.yes.y[2,2]
    
    #flag for number of zero cells in crosstabs for classes by mediator by outcome
    zero.cell[1,1] <- length(which(xmy == 0))
    ###################################################
    
    #fit the firth logistic regression model for P(Y|X,M) to obtain updated parameter estimates now classes have been imputed
    if (zero.cell<2) model.y <- logistf(y~m+x1+x2+x3+int1+int2+int3, data=all.data)
    #if there is more than 1 zero cell, even firth logistic regression does not converge, therefore it is necessary to remove XM interactions from regression model
    if (zero.cell>1) model.y <- logistf(y~m+x1+x2+x3, data=all.data)
    
    #fit the logistic regression model for P(M|X) to obtain updated parameter estimates
    model.m <- logistf(m~x1+x2+x3, data=all.data)
    
    #perturbing the beta coefficients once around coefficients in models using variance-covariance matrix 
    if (zero.cell>1) beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y)),0,0,0)
    if (zero.cell<2) beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y))) 
    beta.m<-c(rmnorm(1,coef(model.m),vcov(model.m)))
    
  } #end of iteration loop
  
  ################################################
  #assessing convergence
  
  #create plot of parameters and cell sizes across each iteration
  results.df <- as.data.frame(results)
  
  #remove the first iteration when no data for X (latent classes) 
  results.df <- subset(results.df, iteration>1)
  
  #to check autocorrelation for each beta after burn in of 100 iterations
  b2.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b2.y"]
  b3.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b3.y"]
  b4.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b4.y"]
  b5.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b5.y"]
  b6.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b6.y"]
  b7.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b7.y"]
  b1.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b1.m"]
  b2.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b2.m"]
  b3.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b3.m"]
  
  #can use a separate folder for updated PCD because of large number of files created (e.g., 80 imputed datasets)
  setwd(upcd.files)
  
  #create a pdf showing traceplots, histograms and autocorrelation plots for each beta, and traceplots and histograms for cell sizes across all iterations
  
  pdf(file=paste0("trace_plots", i, ".pdf")) 
  plot(results.df$iteration, results.df$b2.y,
       xlab = "iteration",
       ylab = "b2.y (y on class 1 vs class 4)")
  lines(results.df$iteration, results.df$b2.y)
  hist (results.df$b2.y, xlab = "b2.y (y on class 1 vs class 4)")
  acf(b2.y, lag.max=50)
  plot(results.df$iteration, results.df$b3.y,
       xlab = "iteration",
       ylab = "b3.y (y on class 2 vs class 4)")
  lines(results.df$iteration, results.df$b3.y)
  hist (results.df$b3.y, xlab = "b3.y (y on class 2 vs class 4)")
  acf(b3.y, lag.max=50)
  plot(results.df$iteration, results.df$b4.y,
       xlab = "iteration",
       ylab = "b4.y (y on class 3 vs class 4)")
  lines(results.df$iteration, results.df$b4.y)
  hist (results.df$b4.y, xlab = "b4.y (y on class 3 vs class 4)")
  acf(b4.y, lag.max=50)
  plot(results.df$iteration, results.df$b5.y,
       xlab = "iteration",
       ylab = "b5.y (y on class 1 x m interaction)")
  lines(results.df$iteration, results.df$b5.y)
  hist (results.df$b5.y, xlab = "b5.y (y on class 1 x m interaction)")
  acf(b5.y, lag.max=50)
  plot(results.df$iteration, results.df$b6.y,
       xlab = "iteration",
       ylab = "b6.y (y on class 2 x m interaction)")
  lines(results.df$iteration, results.df$b6.y)
  hist (results.df$b6.y, xlab = "b6.y (y on class 2 x m interaction)")
  acf(b6.y, lag.max=50)
  plot(results.df$iteration, results.df$b7.y,
       xlab = "iteration",
       ylab = "b7.y (y on class 3 x m interaction)")
  lines(results.df$iteration, results.df$b7.y)
  hist (results.df$b7.y, xlab = "b7.y (y on class 3 x m interaction)")
  acf(b7.y, lag.max=50)
  plot(results.df$iteration, results.df$b1.m,
       xlab = "iteration",
       ylab = "b1.m (m on class 1 vs class 4)")
  lines(results.df$iteration, results.df$b1.m)
  hist (results.df$b1.m, xlab = "b1.m (m on class 1 vs class 4)")
  acf(b1.m, lag.max=50)
  plot(results.df$iteration, results.df$b2.m,
       xlab = "iteration",
       ylab = "b2.m (m on class 2 vs class 4)")
  lines(results.df$iteration, results.df$b2.m)
  hist (results.df$b2.m, xlab = "b2.m (m on class 2 vs class 4)")
  acf(b2.m, lag.max=50)
  plot(results.df$iteration, results.df$b3.m,
       xlab = "iteration",
       ylab = "b3.m (m on class 3 vs class 4")
  lines(results.df$iteration, results.df$b3.m)
  hist (results.df$b3.m, xlab = "b3.m (m on class 3 vs class 4)")
  acf(b3.m, lag.max=50)
  plot(results.df$iteration, results.df$x1.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=1 m=0 y=0")
  lines(results.df$iteration, results.df$x1.m0.y0)
  hist (results.df$x1.m0.y0, xlab = "cell size x=1 m=0 y=0")
  plot(results.df$iteration, results.df$x1.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=1 m=1 y=0")
  lines(results.df$iteration, results.df$x1.m1.y0)
  hist (results.df$x1.m1.y0, xlab = "cell size x=1 m=1 y=0")
  plot(results.df$iteration, results.df$x2.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=2 m=0 y=0")
  lines(results.df$iteration, results.df$x2.m0.y0)
  hist (results.df$x2.m0.y0, xlab = "cell size x=2 m=0 y=0")
  plot(results.df$iteration, results.df$x2.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=2 m=1 y=0")
  lines(results.df$iteration, results.df$x2.m1.y0)
  hist (results.df$x2.m1.y0, xlab = "cell size x=2 m=1 y=0")
  plot(results.df$iteration, results.df$x3.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=3 m=0 y=0")
  lines(results.df$iteration, results.df$x3.m0.y0)
  hist (results.df$x3.m0.y0, xlab = "cell size x=3 m=0 y=0")
  plot(results.df$iteration, results.df$x3.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=3 m=1 y=0")
  lines(results.df$iteration, results.df$x3.m1.y0)
  hist (results.df$x3.m1.y0, xlab = "cell size x=3 m=1 y=0")
  plot(results.df$iteration, results.df$x4.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=4 m=0 y=0")
  lines(results.df$iteration, results.df$x4.m0.y0)
  hist (results.df$x4.m0.y0, xlab = "cell size x=4 m=0 y=0")
  plot(results.df$iteration, results.df$x4.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=4 m=1 y=0")
  lines(results.df$iteration, results.df$x4.m1.y0)
  hist (results.df$x4.m1.y0, xlab = "cell size x=4 m=1 y=0")
  plot(results.df$iteration, results.df$x1.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=1 m=0 y=1")
  lines(results.df$iteration, results.df$x1.m0.y1)
  hist (results.df$x1.m0.y1, xlab = "cell size x=1 m=0 y=1")
  plot(results.df$iteration, results.df$x1.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=1 m=1 y=1")
  lines(results.df$iteration, results.df$x1.m1.y1)
  hist (results.df$x1.m1.y1, xlab = "cell size x=1 m=1 y=1")
  plot(results.df$iteration, results.df$x2.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=2 m=0 y=1")
  lines(results.df$iteration, results.df$x2.m0.y1)
  hist (results.df$x2.m0.y1, xlab = "cell size x=2 m=0 y=1")
  plot(results.df$iteration, results.df$x2.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=2 m=1 y=1")
  lines(results.df$iteration, results.df$x2.m1.y1)
  hist (results.df$x2.m1.y1, xlab = "cell size x=2 m=1 y=1")
  plot(results.df$iteration, results.df$x3.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=3 m=0 y=1")
  lines(results.df$iteration, results.df$x3.m0.y1)
  hist (results.df$x3.m0.y1, xlab = "cell size x=3 m=0 y=1")
  plot(results.df$iteration, results.df$x3.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=3 m=1 y=1")
  lines(results.df$iteration, results.df$x3.m1.y1)
  hist (results.df$x3.m1.y1, xlab = "cell size x=3 m=1 y=1")
  plot(results.df$iteration, results.df$x4.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=4 m=0 y=1")
  lines(results.df$iteration, results.df$x4.m0.y1)
  hist (results.df$x4.m0.y1, xlab = "cell size x=4 m=0 y=1")
  plot(results.df$iteration, results.df$x4.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=4 m=1 y=1")
  lines(results.df$iteration, results.df$x4.m1.y1)
  hist (results.df$x4.m1.y1, xlab = "cell size x=4 m=1 y=1")
  dev.off()
  
  #analysis
  
  #combine imputed class membership stored in "imp.n" with original data
  #prepare 1 mplus .dat file for each imputed dataset to run mediation model (80 .dat files should be created)
  for(l in 1:imp.n) {
    imp.subset <- cbind(data.original,subset(imp, imp[,1]==l))
    prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))}
  
  #create "imp.txt" file for mplus to call imputed datasets
  imp.txt <- matrix(NA,imp.n,1)
  for(l in 1:imp.n) {
    imp.txt[l,1] <- paste0("imp_", l, ".dat")}
  write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #run mediation model using updated PCD
  runModels("2i pcd mediation.inp")
  
  #read in parameters from updated PCD mediation model 
  est.pcd <- readModels("2i pcd mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "upcd.all"
  pcd.all[i,(match("tot.1v4",colnames(pcd.all))):(match("pde.4v1",colnames(pcd.all)))] <- est.pcd[100:135,"est"] # TOT_1V4 to PDE_4V1
  pcd.all[i,(match("tot.1v4.se",colnames(pcd.all))):(match("pde.4v1.se",colnames(pcd.all)))] <- est.pcd[100:135,"se"] # TOT_1V4 to PDE_4V1
  pcd.all[i,(match("p1",colnames(pcd.all))):(match("p4",colnames(pcd.all)))] <- est.pcd[72:75,"est"]  # P_X1 to P_X4
  
  #flags for potential issues: "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp"
  
  #flag if fixed threshold in unconditional model which results in missing value in tech3 covariance matrix of parameters
  pcd.all[i,"fixed.th"] <- 0 #set flag to 0 to start
  for(k in 1:(indicators*latent.classes+(latent.classes-1))*(indicators*latent.classes+(latent.classes-1))) { #size of tech3 matrix is 23x23 (529) 
    if (tech3[k]==999) pcd.all[i,"fixed.th"] <- 1}
  
  #when perturbing within-class thresholds, those with a large standard error can go out of bounds (e.g., corresponding to a probability that is not between 0 and 100%)
  #we created a flag for this in "results" so here we will move this into "pcd.all" and record whether this was the case in any iteration
  pcd.all[i,"large.th"] <- 0 #set flag to 0 to start
  results[1,"large.th"]<-0 
  for(k in 1:(cycles*imp.n+burnin)) {
    if (results[k,"large.th"]==1) pcd.all[i,"large.th"] <- 1}
  
  #we also want a flag so that we know the largest standard error for a threshold in the unconditional model (and the threshold this SE corresponds to)
  #this flag captures the largest SE
  pcd.all[i,"largest.th.se"] <- max(coef.se.original[,"se"])
  #this flag records the threshold that corresponds to largest SE
  pcd.all[i,"largest.th"] <- coef.se.original[which.max(coef.se.original[,"se"]),"est"]
  
  #we created a flag for number of zero cells in the crosstabs for classes by mediator by outcome in "results"
  #we will move this flag to "pcd.all" and record the largest number of zero cells in any iteration
  results[1,"zero.cell"]<-0 #change row 1 in "results" to 0 as latent classes not imputed in first iteration
  pcd.all[i,"zero.cell"] <- max(results[,"zero.cell"])
  
  #we will also flag if there was a zero cell (in the crosstabs for classes by mediator by outcome) in any of the iterations that we saved as one of the 80 imputed datasets
  zero.cell.imp <- matrix(NA,1,imp.n)
  for(l in 1:imp.n) {
    zero.cell.imp[1,l] <- results[cycles*l+burnin,"zero.cell"]}
  
  pcd.all[i,"zero.cell.imp"] <- 0 #set flag to 0 to start
  for(k in 1:imp.n) {
    if (zero.cell.imp[1,k]>0) pcd.all[1,"zero.cell.imp"] <- 1}
  
  cat(i,"\n")} #end of simulations loop

#############################################################

#saving results for comparisons of interest (e.g. with low class as the reference group) and exporting out mediation effects for each method

#onestep
onestep <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes+1)
colnames(onestep) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                       "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4","largest.th.se")

#move over class probabilities
onestep[,"p1"]<-onestep.all[,"p1"]
onestep[,"p2"]<-onestep.all[,"p2"]
onestep[,"p3"]<-onestep.all[,"p3"]
onestep[,"p4"]<-onestep.all[,"p4"]
onestep[,"largest.th.se"]<-onestep.all[,"largest.th.se"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory parameters (auc1-auc4))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc3"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc3"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc2"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc2"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#bch
bch <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(bch) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                   "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
bch[,"p1"]<-bch.all[,"p1"]
bch[,"p2"]<-bch.all[,"p2"]
bch[,"p3"]<-bch.all[,"p3"]
bch[,"p4"]<-bch.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}


#modal
modal <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(modal) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                     "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
modal[,"p1"]<-modal.all[,"p1"]
modal[,"p2"]<-modal.all[,"p2"]
modal[,"p3"]<-modal.all[,"p3"]
modal[,"p4"]<-modal.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#npcd
npcd <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(npcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                    "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
npcd[,"p1"]<-npcd.all[,"p1"]
npcd[,"p2"]<-npcd.all[,"p2"]
npcd[,"p3"]<-npcd.all[,"p3"]
npcd[,"p4"]<-npcd.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#incpcd
incpcd <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(incpcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                      "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
incpcd[,"p1"]<-incpcd.all[,"p1"]
incpcd[,"p2"]<-incpcd.all[,"p2"]
incpcd[,"p3"]<-incpcd.all[,"p3"]
incpcd[,"p4"]<-incpcd.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low) 
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc3"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc3"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc2"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc2"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#pcd
pcd <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes+6)
colnames(pcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                   "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4",
                   "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp")

#move over class probabilities
pcd[,"p1"]<-pcd.all[,"p1"]
pcd[,"p2"]<-pcd.all[,"p2"]
pcd[,"p3"]<-pcd.all[,"p3"]
pcd[,"p4"]<-pcd.all[,"p4"]
pcd[,"fixed.th"]<-pcd.all[,"fixed.th"]
pcd[,"large.th"]<-pcd.all[,"large.th"]
pcd[,"largest.th"]<-pcd.all[,"largest.th"]
pcd[,"largest.th.se"]<-pcd.all[,"largest.th.se"]
pcd[,"zero.cell"]<-pcd.all[,"zero.cell"]
pcd[,"zero.cell.imp"]<-pcd.all[,"zero.cell.imp"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#put results in a table and export to a results folder
setwd(results.files)

write.table(onestep, file="table.onestep.txt", sep = ",")
write.table(bch, file="table.bch.txt", sep = ",")
write.table(modal, file="table.modal.txt", sep = ",")
write.table(npcd, file="table.npcd.txt", sep = ",")
write.table(incpcd, file="table.incpcd.txt", sep = ",")
write.table(pcd, file="table.pcd.txt", sep = ",")

write.table(entropy, file="table.entropy.txt", sep = ",")
write.table(auc, file="table.auc.txt", sep = ",")

################################################end of script#####################################################################

############### GOOD ENTROPY ##########################

############### DEFINE DIRECTORIES ##########################

#Folder where simulation data are saved:
sim.data <- "file location/all sims/data/good entropy"
#Folder where Mplus input files (uncond, onestep, bch, modal and inc) are saved:
inp.files <- "file location/all sims"
#*Folder where non-inclusive PCD files are saved:
npcd.files <- "file location/all sims"
#*Folder where inclusive PCD files are saved:
incpcd.files <- "file location/all sims"
#*Folder where updated PCD files are saved:
upcd.files <- "file location/all sims"
#*Folder where results are saved:
results.files <- "file location/all sims/results/good entropy"

#script for simulated datasets with good entropy (based on ALSPAC)
#############################################################
#X = 4 class latent nominal exposure (conduct trajectories: early onset persistent, adolescent onset, childhood limited and low)
latent.classes <- 4
#M = binary mediator (peer deviance)
#Y = binary outcome (problematic alcohol use)
#U1-U5 = binary latent class indicators (conduct problems from age 4 to 13 years)
indicators <- 5

############### PREPARE MATRICES ##########################

#here we are running 500 simulated datasets
sims <- 500
#simulated data has a sample size of 5000
sample <- 5000
#we wish to compare each latent class with all the others 
class.comp <- 4*3
#we wish to report 3 mediation effects (total, indirect - tie, and direct - pde) and their standard errors
mediation.effects <- 3*2

#create matrices to store the mediation results for each method
#number of rows = number of simulated datasets (here we use 500)
#number of columns is based on number of class comparisons (12), number of mediation effects (total, tie and pde) and SEs (6) and 
#number of extra parameters (e.g. area under the trajectory, class probabilities, or additional checks)

onestep.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes*2+1))

#labels represent the mediation effects (total, tie and pde) for each class comparison, area under the trajectory for each class, and class probabilities
colnames(onestep.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                           "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                           "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                           "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                           "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se", 
                           "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se", 
                           "auc1","auc2","auc3","auc4","largest.th.se","p1","p2","p3","p4")

bch.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(bch.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                       "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                       "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                       "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                       "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                       "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                       "p1","p2","p3","p4")

modal.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(modal.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                         "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                         "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                         "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                         "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                         "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                         "p1","p2","p3","p4")

npcd.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(npcd.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                        "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                        "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                        "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                        "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                        "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                        "p1","p2","p3","p4")

incpcd.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(incpcd.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                          "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                          "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                          "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                          "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                          "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                          "p1","p2","p3","p4")

pcd.all <- matrix(NA,sims,(class.comp*mediation.effects+latent.classes+6))

#labels represent the mediation effects (total, tie and pde) for each class comparison, 6 flags for potential issues and class probabilities 
colnames(pcd.all) <-  c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                        "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                        "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                        "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                        "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                        "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se", "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp",
                        "p1","p2","p3","p4")

#create matrix to store the entropy of the unconditional model
entropy <- matrix(NA,sims,1)

#create matrix to store the area under the trajectory for each class in the unconditional model (this allows us to know the order of the classes)
auc <- matrix(NA,sims,latent.classes)
colnames(auc) <- c("auc1","auc2","auc3","auc4")

#create matrix to store the area under the trajectory for each class in the inclusive latent class model (this allows us to know the order of the classes)
auc.inc <- matrix(NA,sims,latent.classes)
colnames(auc.inc) <- c("auc1","auc2","auc3","auc4")

############### READ IN THE DATA ##########################

#run loop over 500 simulated datasets 
for(i in 1:sims) {
  
  #set working directory using file paths saved at the start
  setwd(sim.data)
  
  #read in Mplus .dat file with simulated data (sim1.dat)
  data.original <- read.table(file=paste0("sim", i, ".dat"), sep="", header=FALSE, col.names = c("y", #outcome
                                                                                                 "u1","u2","u3","u4","u5", #latent class indicators
                                                                                                 "m", #mediator
                                                                                                 "c" #modal class assignment for the exposure
  ))
  
  #replicate the original data and add empty columns to store imputed latent class membership for 4 classes (once generated) 
  #and interactions between each latent class and the mediator
  all.data <- data.original
  all.data$x1 <- NA #empty column to add dummy code for membership in latent class 1
  all.data$x2 <- NA #empty column to add dummy code for membership in latent class 2
  all.data$x3 <- NA #empty column to add dummy code for membership in latent class 3
  all.data$x4 <- NA #empty column to add dummy code for membership in latent class 4
  all.data$int1 <- NA  #empty column to add interaction between latent class 1 and mediator
  all.data$int2 <- NA  #empty column to add interaction between latent class 2 and mediator
  all.data$int3 <- NA  #empty column to add interaction between latent class 3 and mediator
  all.data$int4 <- NA  #empty column to add interaction between latent class 4 and mediator
  
  ############### UNCONDITIONAL LATENT CLASS MODEL ##########################
  
  #set working directory using file paths saved at the start
  setwd(inp.files)
  
  #this is simply renaming the Mplus .dat file and saving in the folder with the Mplus input files
  #this step is important when running many simulations, but not necessary otherwise
  #sim.dat will be written over each time a new simulated dataset is analysed
  prepareMplusData(data.original,"sim.dat")
  
  #run the unconditional latent class model
  runModels("2b uncond latent class.inp")
  
  #read in the Mplus output file and name it "model_output" to use later
  model_output <- readModels("2b uncond latent class.out")
  
  #save parameters from unconditional latent class model (within-class thresholds for latent class indicators and class intercepts) and their SEs
  coef.se.original <- model_output$parameters$`unstandardized`[1:(indicators*latent.classes+(latent.classes-1)),3:4] 
  #labels represent five within-class thresholds across 4 classes and 3 class intercepts)
  rownames(coef.se.original) <- c("th1.c1","th2.c1","th3.c1","th4.c1","th5.c1", #within-class thresholds for class 1
                                  "th1.c2","th2.c2","th3.c2","th4.c2","th5.c2", #within-class thresholds for class 2
                                  "th1.c3","th2.c3","th3.c3","th4.c3","th5.c3", #within-class thresholds for class 3
                                  "th1.c4","th2.c4","th3.c4","th4.c4","th5.c4", #within-class thresholds for class 4
                                  "int.c1","int.c2","int.c3") #latent class intercepts
  
  #save model entropy into the matrix we created at the start of the script
  entropy[i,1] <- model_output$summaries$`Entropy`
  
  #create a matrix of class probabilities
  p.original <- matrix(NA,latent.classes,1)
  #these can be calculated using the 3 class intercepts we saved in "coef.se.original"
  p.original[1,1] <- exp(coef.se.original["int.c1","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
  p.original[2,1] <- exp(coef.se.original["int.c2","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
  p.original[3,1] <- exp(coef.se.original["int.c3","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))  
  p.original[4,1] <- 1/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
  
  #derive area under the trajectory parameters so we know the order of the latent classes in the unconditional model
  #these can be calculated using the within-class thresholds we saved in "coef.se.original"
  auc[i,"auc1"] <- exp(-1*coef.se.original["th1.c1","est"])/(1+exp(-1*coef.se.original["th1.c1","est"]))+2*exp(-1*coef.se.original["th2.c1","est"])/(1+exp(-1*coef.se.original["th2.c1","est"]))+3*exp(-1*coef.se.original["th3.c1","est"])/(1+exp(-1*coef.se.original["th3.c1","est"]))+4*exp(-1*coef.se.original["th4.c1","est"])/(1+exp(-1*coef.se.original["th4.c1","est"]))+5*exp(-1*coef.se.original["th5.c1","est"])/(1+exp(-1*coef.se.original["th5.c1","est"]))
  auc[i,"auc2"] <- exp(-1*coef.se.original["th1.c2","est"])/(1+exp(-1*coef.se.original["th1.c2","est"]))+2*exp(-1*coef.se.original["th2.c2","est"])/(1+exp(-1*coef.se.original["th2.c2","est"]))+3*exp(-1*coef.se.original["th3.c2","est"])/(1+exp(-1*coef.se.original["th3.c2","est"]))+4*exp(-1*coef.se.original["th4.c2","est"])/(1+exp(-1*coef.se.original["th4.c2","est"]))+5*exp(-1*coef.se.original["th5.c2","est"])/(1+exp(-1*coef.se.original["th5.c2","est"]))
  auc[i,"auc3"] <- exp(-1*coef.se.original["th1.c3","est"])/(1+exp(-1*coef.se.original["th1.c3","est"]))+2*exp(-1*coef.se.original["th2.c3","est"])/(1+exp(-1*coef.se.original["th2.c3","est"]))+3*exp(-1*coef.se.original["th3.c3","est"])/(1+exp(-1*coef.se.original["th3.c3","est"]))+4*exp(-1*coef.se.original["th4.c3","est"])/(1+exp(-1*coef.se.original["th4.c3","est"]))+5*exp(-1*coef.se.original["th5.c3","est"])/(1+exp(-1*coef.se.original["th5.c3","est"]))
  auc[i,"auc4"] <- exp(-1*coef.se.original["th1.c4","est"])/(1+exp(-1*coef.se.original["th1.c4","est"]))+2*exp(-1*coef.se.original["th2.c4","est"])/(1+exp(-1*coef.se.original["th2.c4","est"]))+3*exp(-1*coef.se.original["th3.c4","est"])/(1+exp(-1*coef.se.original["th3.c4","est"]))+4*exp(-1*coef.se.original["th4.c4","est"])/(1+exp(-1*coef.se.original["th4.c4","est"]))+5*exp(-1*coef.se.original["th5.c4","est"])/(1+exp(-1*coef.se.original["th5.c4","est"]))
  
  #import tech3 as a covariance matrix to allow us to use the covariance between parameters for perturbing in updated PCD
  #############################################################
  #NB the order of parameters in this covariance matrix is different to order in "coef.se.original" (see tech1) so this will need to be addressed before perturbing parameters
  #this is due to no missing data in the indicators which affects the numbering of parameters in tech1
  #############################################################
  #this outputs paramCov as a matrix object in R 
  tech3 <- model_output$tech3$paramCov
  #this creates a full, symmetrical covariance matrix
  upperTriangle(tech3) <- lowerTriangle(tech3, byrow=TRUE)
  cov <- tech3
  
  #address any fixed parameters: within class thresholds that have been fixed at 15 or -15 (representing a probability of 0 or 100%) do not have (co)variances
  #the code below means that fixed parameters will still get perturbed in updated PCD but only a very small amount
  #replace missing (999) in covariance matrix with 0 to represent no covariance for the fixed parameters
  cov[cov==999] <- 0
  #change the variance for fixed parameters to be very small (0.000000001)
  for(j in 1:(indicators*latent.classes+(latent.classes-1))) {
    if (cov[j,j]==0) cov[j,j]<- 0.000000001}
  
  ###################### ONESTEP MODEL ##########################
  
  #run the onestep latent class mediation model
  runModels("2c onestep mediation.inp")
  
  #read in parameters from one-step latent class model    
  est.onestep <- readModels("2c onestep mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects, area under the trajectory for each class, and class probabilities) and their SEs in the matrix created at the start "onestep.all"
  onestep.all[i, (match("tot.1v4",colnames(onestep.all))):(match("pde.4v1",colnames(onestep.all)))] <- est.onestep[116:151,"est"] # TOT_1V4 to PDE_4V1
  onestep.all[i, (match("tot.1v4.se",colnames(onestep.all))):(match("pde.4v1.se",colnames(onestep.all)))] <- est.onestep[116:151,"se"] # TOT_1V4 to PDE_4V1
  onestep.all[i, (match("auc1",colnames(onestep.all))):(match("auc4",colnames(onestep.all)))] <- est.onestep[152:155,"est"] # AUC
  onestep.all[i, (match("p1",colnames(onestep.all))):(match("p4",colnames(onestep.all)))] <- est.onestep[88:91,"est"]  # P_X1 to P_X4
  
  #save standard errors for thresholds from onestep model
  se <- est.onestep[c(3:7,17:21,31:35,45:49),4] 
  
  #we want a flag so that we know the largest standard error for a threshold in the onestep model
  #this flag captures the largest SE
  onestep.all[i,"largest.th.se"] <- max(se)
  
  ######################## BCH ##########################
  
  #bch can only be used with multiple latent classes if they are combined into one class
  #use one 8 class model for X and M and set up mediation model using nom-nom-cat approach (because bch needs to include an auxiliary variable)
  
  #read in bch weights and modal class assignment which were exported out of the unconditional latent class model in "bch.txt"
  data.bch <- read.table("bch.txt", sep="", na.strings="*", header=FALSE)
  #select the weights for each class and the modal class assignment
  data.bch <- data.bch[,c(6:9,14)]
  #labels representing one weight for each class and the modal class assignment
  colnames(data.bch) <- c("bch1","bch2","bch3","bch4","modal")
  #combine the bch weights with the mediator and outcome from the original data
  data.bch <- cbind(data.original[,c("y","m")], data.bch)
  
  #create weights for 8 class model by multiplying the bch weight for each class with the observed data on the mediator
  data.bch$bch11 <- data.bch$bch1*data.bch$m
  data.bch$bch12 <- data.bch$bch1*(1-data.bch$m)
  data.bch$bch21 <- data.bch$bch2*data.bch$m
  data.bch$bch22 <- data.bch$bch2*(1-data.bch$m)
  data.bch$bch31 <- data.bch$bch3*data.bch$m
  data.bch$bch32 <- data.bch$bch3*(1-data.bch$m)
  data.bch$bch41 <- data.bch$bch4*data.bch$m
  data.bch$bch42 <- data.bch$bch4*(1-data.bch$m)
  
  #prepare Mplus .dat file using "data.bch"  
  prepareMplusData(data.bch,"bch.dat")
  #run mediation model using bch method
  runModels("2d bch mediation.inp")
  
  #read in parameters from bch latent class mediation model   
  est.bch <- readModels("2d bch mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "bch.all"
  bch.all[i,(match("tot.1v4",colnames(bch.all))):(match("pde.4v1",colnames(bch.all)))] <- est.bch[68:103,"est"] # TOT_1V4 to PDE_4V1
  bch.all[i,(match("tot.1v4.se",colnames(bch.all))):(match("pde.4v1.se",colnames(bch.all)))] <- est.bch[68:103,"se"] # TOT_1V4 to PDE_4V1
  bch.all[i,(match("p1",colnames(bch.all))):(match("p4",colnames(bch.all)))] <- est.bch[40:43,"est"]  # P_X1 to P_X4
  
  ######################## MODAL CLASS ASSIGNMENT ##########################
  
  #run mediation model using modal class assignment
  runModels("2e modal mediation.inp")
  
  #read in parameters from modal class mediation model  
  est.modal <- readModels("2e modal mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "modal.all"
  modal.all[i,(match("tot.1v4",colnames(modal.all))):(match("pde.4v1",colnames(modal.all)))] <- est.modal[100:135,"est"] # TOT_1V4 to PDE_4V1
  modal.all[i,(match("tot.1v4.se",colnames(bch.all))):(match("pde.4v1.se",colnames(bch.all)))] <- est.modal[100:135,"se"] # TOT_1V4 to PDE_4V1
  modal.all[i,(match("p1",colnames(bch.all))):(match("p4",colnames(bch.all)))] <- est.modal[72:75,"est"]  # P_X1 to P_X4
  
  ######################## NON-INCLUSIVE PCD ##########################
  
  #set seed for non-inclusive PCD
  #this will set the seed for each session in R (not each time run loops)
  #to make sure you get the same results each time, run the loop only once within each R session
  set.seed(80)
  
  #calculate P(X=x|U) - this will give us the probability of class membership (cprobs) from unconditional latent class model
  
  #for each class, multiply the individual data (responses to 5 binary indicators: U1 to U5) with within-class thresholds
  #number of latent class indicators
  #class 1
  #this creates a matrix with within-class thresholds for class 1 repeated for every individual in the dataset (e.g. repeated 5000 times)
  theta <- matrix(rep(coef.se.original[(match("th1.c1",rownames(coef.se.original))):(match("th5.c1",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
  #if latent class indicator is present (-1)*threshold is used, if latent class indicator is absent (1)*threshold is used
  P1<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  #class 2
  #this creates a matrix with within-class thresholds for class 2 repeated for every individual in the dataset (e.g. repeated 5000 times)
  theta <- matrix(rep(coef.se.original[(match("th1.c2",rownames(coef.se.original))):(match("th5.c2",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
  P2<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  #class 3
  #this creates a matrix with within-class thresholds for class 3 repeated for every individual in the dataset (e.g. repeated 5000 times)
  theta <- matrix(rep(coef.se.original[(match("th1.c3",rownames(coef.se.original))):(match("th5.c3",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
  P3<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  #class 4
  #this creates a matrix with within-class thresholds for class 4 repeated for every individual in the dataset (e.g. repeated 5000 times)
  theta <- matrix(rep(coef.se.original[(match("th1.c4",rownames(coef.se.original))):(match("th5.c4",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
  P4<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
  
  #multiply the 5 probabilities (for each latent class indicator) in each row of P* then multiply by class probabilities that were saved earlier in "p.original"  
  N1<-apply(P1,1,prod)*p.original[1,1]
  N2<-apply(P2,1,prod)*p.original[2,1]
  N3<-apply(P3,1,prod)*p.original[3,1]
  N4<-apply(P4,1,prod)*p.original[4,1]    
  
  #this gives us the probability of class membership for each person in the dataset 
  P<-cbind(N1/(N1+N2+N3+N4),N2/(N1+N2+N3+N4),N3/(N1+N2+N3+N4),N4/(N1+N2+N3+N4))
  
  #we could have simply read these in from "bch.txt" (as we did with the bch weights) as they were exported from the unconditional latent class model 
  colnames(P) <-  c("cprob1","cprob2","cprob3","cprob4")
  
  #we will now impute class membership 40 times for each person using their probability of class membership
  #we will create 40 imputed datasets (chosen to keep Monte Carlo error at less than 10% of standard error for parameters from regression model for Y)
  imp.n <- 40
  #create a matrix to store imputed class membership for each person (sample*imp.n=5000*40)
  imp <- matrix(NA,sample*imp.n,latent.classes+1)
  colnames(imp) <-  c("imp","x1","x2","x3","x4")
  
  #first column is simply an indicator for imputation number (range from 1 to 40)
  for(h in 0:(imp.n-1)) {
    imp [c(sample*h+1:sample),1]<-h+1}
  
  #for each of 40 imputations we will create a matrix to store class membership, 
  #then we will use the probabilities of class membership "P" to randomly assign each individual to a class (X = 1,..k).
  #and we will add their imputed class membership into the matrix "imp" we created earlier
  for(j in 1:imp.n) {
    x <- matrix(NA,sample,latent.classes)    
    for(k in 1:sample) {x[k,1:latent.classes] <- rmultinom(1,1,P[k,])}
    for(l in 1:imp.n) {
      if (j==l) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x}
  } 
  
  #analysis
  #can use a separate folder for non-inclusive PCD because of large number of files created (e.g., 40 imputed datasets)
  setwd(npcd.files) 
  
  #combine imputed class membership stored in "imp.n" with original data
  #prepare 1 mplus .dat file for each imputed dataset to run mediation model (40 .dat files should be created)
  for(l in 1:imp.n) {
    imp.subset <- cbind(data.original,subset(imp, imp[,1]==l))
    prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))
  } 
  
  #create "imp.txt" file for mplus to call imputed datasets
  imp.txt <- matrix(NA,imp.n,1)
  for(l in 1:imp.n) {
    imp.txt[l,1] <- paste0("imp_", l, ".dat")
  } 
  write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #run mediation model using non-inclusive PCD
  runModels("2f npcd mediation.inp")
  
  #read in parameters from npcd mediation model 
  est.npcd <- readModels("2f npcd mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "npcd.all"
  npcd.all[i,(match("tot.1v4",colnames(npcd.all))):(match("pde.4v1",colnames(npcd.all)))] <- est.npcd[100:135,"est"] # TOT_1V4 to PDE_4V1
  npcd.all[i,(match("tot.1v4.se",colnames(npcd.all))):(match("pde.4v1.se",colnames(npcd.all)))] <- est.npcd[100:135,"se"] # TOT_1V4 to PDE_4V1
  npcd.all[i,(match("p1",colnames(npcd.all))):(match("p4",colnames(npcd.all)))] <- est.npcd[72:75,"est"]  # P_X1 to P_X4
  
  ######################## INCLUSIVE PCD ##########################
  
  #set seed for inclusive PCD
  set.seed(81)
  
  #set working directory using file paths saved at the start
  setwd(inp.files)
  
  #run inclusive latent class model and export out cprobs
  runModels("2g inc latent class.inp")
  
  #read in auc parameters from inclusive latent class model
  est.inc <- readModels("2g inc latent class.out", what="parameters")$parameters$`unstandardized`
  auc.inc[i,(match("auc1",colnames(auc.inc))):(match("auc4",colnames(auc.inc)))] <- est.inc[30:33,"est"]  # AUC1 to AUC4
  
  #read in class probabilities from "inc.txt" (as we did with the bch weights exported from the unconditional latent class model) 
  cprobs <- read.table("inc.txt", sep="", na.strings="*", header=FALSE)
  #we could also extract modal classes from here (v12) to use the inclusive modal approach
  P <- cprobs[,c(8:11)]
  colnames(P) <- c("cprob1","cprob2","cprob3","cprob4")
  
  #we will now impute class membership 40 times for each person using their probability of class membership
  #create a matrix to store imputed class membership for each person (sample*imp.n=5000*40)
  imp <- matrix(NA,sample*imp.n,latent.classes+1)
  colnames(imp) <-  c("imp","x1","x2","x3","x4")
  
  #first column is simply an indicator for imputation number (range from 1 to 40)
  for(h in 0:(imp.n-1)) {
    imp [c(sample*h+1:sample),1]<-h+1}
  
  #for each of 40 imputations we will create a matrix to store class membership, 
  #then we will use the probabilities of class membership "P" to randomly assign each individual to a class (X = 1,..k).
  #and we will add their imputed class membership into the matrix "imp" we created earlier
  for(j in 1:imp.n) {
    x <- matrix(NA,sample,latent.classes)    
    for(k in 1:sample) {x[k,1:latent.classes] <- rmultinom(1,1,P[k,])}
    for(l in 1:imp.n) {
      if (j==l) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x}
  } 
  
  #analysis
  #can use a separate folder for inclusive PCD because of large number of files created (e.g., 40 imputed datasets)
  setwd(incpcd.files) 
  
  #combine imputed class membership stored in "imp.n" with original data
  #prepare 1 mplus .dat file for each imputed dataset to run mediation model (40 .dat files should be created)
  for(l in 1:imp.n) {
    imp.subset <- cbind(data.original,subset(imp, imp[,1]==l))
    prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))
  } 
  
  #create "imp.txt" file for mplus to call imputed datasets
  imp.txt <- matrix(NA,imp.n,1)
  for(l in 1:imp.n) {
    imp.txt[l,1] <- paste0("imp_", l, ".dat")
  } 
  write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #run mediation model using inclusive PCD
  runModels("2h incpcd mediation.inp")
  
  #read in parameters from incpcd mediation model 
  est.incpcd <- readModels("2h incpcd mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "incpcd.all"
  incpcd.all[i,(match("tot.1v4",colnames(incpcd.all))):(match("pde.4v1",colnames(incpcd.all)))] <- est.incpcd[100:135,"est"] # TOT_1V4 to PDE_4V1
  incpcd.all[i,(match("tot.1v4.se",colnames(incpcd.all))):(match("pde.4v1.se",colnames(incpcd.all)))] <- est.incpcd[100:135,"se"] # TOT_1V4 to PDE_4V1
  incpcd.all[i,(match("p1",colnames(incpcd.all))):(match("p4",colnames(incpcd.all)))] <- est.incpcd[72:75,"est"]  # P_X1 to P_X4
  
  #remove any objects from memory that will not be needed for the rest of the script
  rm(model_output, est.onestep, data.bch, est.bch, est.modal, est.npcd, est.incpcd, cprobs, P)
  
  ######################## UPDATED PCD ##########################
  
  #set seed for updated PCD
  set.seed(82)
  
  #STEP 1: 
  #initialise beta by fitting the logistic regression model for Y with beta for classes set to zero 
  #initialise alpha by fitting the logistic regression model for M with alpha for classes set to zero
  #perturb the estimates from these models with Gaussian noise of mean zero and values from the variance-covariance matrix of these parameter estimates
  
  #firth logistic regression is used to address problems with perfect prediction when bringing in exposure-mediator interactions later in the script
  #regress outcome y on mediator m
  model.y <- logistf(y~m, data=all.data)
  #regression model for mediator m 
  model.m <- logistf(m~1, data=all.data)
  
  #perturbing the coefficients once around coefficients in models using variance-covariance matrix 
  #initially use zero for exposure coefficients (and exposure-mediator interactions in outcome model)
  #once latent class exposure has been imputed in subsequent runs, these will become coefficients from models
  beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y)),0,0,0,0,0,0)
  beta.m<-c(rmnorm(1,coef(model.m),vcov(model.m)),0,0,0)
  
  #we will create 80 imputed datasets for class membership (chosen to keep Monte Carlo error at less than 10% of standard error for parameters in regression model for Y)
  imp.n <- 80
  #we will allow 20 iterations between saving out imputed class membership
  cycles <- 20
  #we will allow a burn in of 100 iterations before starting to save out imputed class membership
  burnin <- 100
  
  #create a matrix to save the results from every iteration (all results are recorded to assess convergence later)
  #rows=iterations=(cycles*imp.n+burnin)
  #columns=iteration number, beta coeficients from regression models, cell sizes and 2 flags for issues=(6+(latent.classes-1)*3+latent.classes*4)=31
  results <- matrix(NA,cycles*imp.n+burnin,6+(latent.classes-1)*3+latent.classes*4) 
  
  colnames(results) <- c("iteration", #iteration number
                         "b0.y","b1.y","b2.y","b3.y","b4.y","b5.y","b6.y","b7.y", #coefficients from regression model for the outcome
                         "b0.m","b1.m","b2.m","b3.m", #coefficients from regression model for the mediator
                         "x1.m0.y0","x1.m1.y0","x2.m0.y0","x2.m1.y0","x3.m0.y0","x3.m1.y0","x4.m0.y0","x4.m1.y0","x1.m0.y1","x1.m1.y1","x2.m0.y1","x2.m1.y1","x3.m0.y1","x3.m1.y1","x4.m0.y1","x4.m1.y1", #cell sizes from crosstabs for classes by mediator by outcome
                         "zero.cell", #flag for presence of zero cells in crosstabs
                         "large.th" #flag for a within-class threshold in unconditional latent class model that was out of bounds after perturbing
  ) 
  
  #the first column is simply an indicator of iteration number (range from 1 to 1700)
  results[,1]<-1:(cycles*imp.n+burnin)
  
  #create a matrix to store cell sizes from crosstabs for classes by mediator by outcome 
  xmy<-matrix(NA,1,latent.classes*4)
  colnames(xmy) <-  c("x1m0y0","x1m1y0","x2m0y0","x2m1y0","x3m0y0","x3m1y0","x4m0y0","x4m1y0","x1m0y1","x1m1y1","x2m0y1","x2m1y1","x3m0y1","x3m1y1","x4m0y1","x4m1y1")
  #create a matrix to store presence of zero cells in this crosstabs
  zero.cell<-matrix(NA,1,1)
  #create a matrix to store presence of a within-class threshold in unconditional latent class model that was out of bounds after perturbing (e.g., not corresponding to 0 to 100% probability) 
  large.th<-matrix(NA,1,1)
  
  #create a matrix to store imputed class membership for each person (sample*imp.n=5000*80)
  imp <- matrix(NA,sample*imp.n,latent.classes+1)
  colnames(imp) <- c("imp","x1","x2","x3","x4") 
  #first column is simply an indicator for imputation number (range from 1 to 80)
  for(h in 0:(imp.n-1)) {
    imp [c(sample*h+1:sample),1]<-h+1}
  
  #we will now create a loop to repeat steps 2 and 3 below 1700 (imp.n*cycles+burnin) times and save estimates after every 20 iterations (after a 100 iteration burn in)
  for(j in 1:(imp.n*cycles+burnin)) {
    results[j,2:9]<-beta.y #save the perturbed coefficients from regression model for the outcome "beta.y"
    results[j,10:13]<-beta.m #save the perturbed coefficients from regression model for the mediator "beta.m"    
    results[j,14:29]<-xmy #save cell sizes from crosstabs for classes by mediator by outcome "xmy"
    results[j,30]<-zero.cell #save a flag for presence of zero cells in this crosstabs "zero.cell"
    results[j,31]<-large.th #save a flag for a within-class threshold in unconditional latent class model that was out of bounds after peturbing "large.th"
    
    #this saves each individual's imputed class membership 80 times (generated later in the script) into the matrix "imp" we created earlier     
    for(l in 1:imp.n) {
      if (j==(cycles*l+burnin)) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x
    } 
    
    #STEP 2a:
    #we will now perturb the parameters (within-class thresholds and class intercepts) that we saved earlier "coef.se.original" from the unconditional latent class model 
    #in order to perturb the parameters using the covariance matrix in tech3 we need to reorder parameters so they match with numbering in tech1
    #because we have no missing data, they are in a different order to what would be expected
    #when using a dataset with missing data on class indicators (as would usually be the case outside of simulated data), this reordering step is not needed
    
    #keep only the parameters (within-class thresholds for latent class indicators and class intercepts) from unconditional latent class model (i.e. drop the SEs)
    coef.mplus.orig<-coef.se.original[,"est"]
    #reorder the parameters so that they are in the same order as is used in the covariances matrix of parameters from tech3
    index<-c(1,5,9,13,17,2,6,10,14,18,3,7,11,15,19,4,8,12,16,20,21,22,23)
    coef.mplus.orig<-coef.mplus.orig[order(index)]
    #perturb these parameters based on their variance-covariance matrix (saved in "cov")
    coef.mplus<-rmnorm(1,coef.mplus.orig,cov)
    #reorder again to preserve original ordering
    index<-c(1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19,5,10,15,20,21,22,23)
    coef<-coef.mplus[order(index)]
    
    #when perturbing within-class thresholds, those with a large standard error can go out of bounds (e.g., corresponding to a probability that is not between 0 and 100%)
    #we will create a flag so we know when this is the case
    large.th[1,1]<-0    
    for(k in 1:(latent.classes*indicators)) {
      if (coef[k]>15||coef[k]<(-15)) large.th[1,1] <- 1} 
    #we will also constrain within-class thresholds to be between -15 and 15 (corresponding to a probability that is between 0 and 100%)
    coef[coef>15] <- 15
    coef[coef<(-15)] <- (-15)
    
    #turning 'coef' into a matrix so can label rows and columns to use in script below (instead of referring to numbers)
    coef<-as.matrix(coef)
    rownames(coef) <- c("th1.c1","th2.c1","th3.c1","th4.c1","th5.c1",
                        "th1.c2","th2.c2","th3.c2","th4.c2","th5.c2",
                        "th1.c3","th2.c3","th3.c3","th4.c3","th5.c3",
                        "th1.c4","th2.c4","th3.c4","th4.c4","th5.c4",
                        "int.c1","int.c2","int.c3")
    colnames(coef) <- c("est")
    
    #now calculate P(X=x|U) - without perturbing this give us the probability of class membership (cprobs) from unconditional latent class model
    
    #for each class, multiply the individual data (responses to 5 binary indicators: U1 to U5) with within-class thresholds
    #class 1
    #this creates a matrix with within-class thresholds for class 1 repeated for every individual in the dataset (e.g. repeated 5000 times)
    theta <- matrix(rep(coef[(match("th1.c1",rownames(coef.se.original))):(match("th5.c1",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
    #if latent class indicator is present (-1)*threshold is used, if latent class indicator is absent (1)*threshold is used
    P1<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
    #class 2
    theta <- matrix(rep(coef[(match("th1.c2",rownames(coef.se.original))):(match("th5.c2",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
    P2<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
    #class 3
    theta <- matrix(rep(coef[(match("th1.c3",rownames(coef.se.original))):(match("th5.c3",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
    P3<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
    #class 4
    theta <- matrix(rep(coef[(match("th1.c4",rownames(coef.se.original))):(match("th5.c4",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
    P4<-exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("u1",colnames(all.data))):(match("u5",colnames(all.data)))]*theta))
    
    #we need to create a new matrix of class probabilities using the perturbed class intercepts
    p <- matrix(NA,latent.classes,1)
    p[1,1] <- exp(coef["int.c1","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
    p[2,1] <- exp(coef["int.c2","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
    p[3,1] <- exp(coef["int.c3","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))  
    p[4,1] <- 1/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
    
    #then multiply "P1" to "P4" by the perturbed class probabilities that we have created above "p" 
    #this gives us the probability of class membership for each person in the dataset
    #these will differ slightly to "cprobs" that can be exported from the unconditional latent class model due to perturbing
    N1<-apply(P1,1,prod)*p[1,1]
    N2<-apply(P2,1,prod)*p[2,1]
    N3<-apply(P3,1,prod)*p[3,1]
    N4<-apply(P4,1,prod)*p[4,1]    
    
    #STEP 2b - Bayes rule:
    
    #combine P(X=x|U) along with coefficients from regression models for outcome and mediator
    #this will calculate probabilities P(X=x|Y,M,U) for each class (x = 1,..k) for each individual
    #we will do this via a number a steps below
    
    #P(Y=1|M,X): calculate probability that the outcome (Y) = 1 given mediator (M) and exposure (X, latent classes) 
    #this uses the coefficients from the regression model for the outcome (intercept (b0), coef for M (b1), coefs for X (b2-b4) and coefs for XM interaction (b5-b7))
    #class 1
    num1 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m+results[j,"b2.y"]+results[j,"b5.y"]*all.data$m)
    Y11 <- num1/(1+num1)  
    #class 2
    num2 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m+results[j,"b3.y"]+results[j,"b6.y"]*all.data$m)
    Y21 <- num2/(1+num2)  
    #class 3
    num3 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m+results[j,"b4.y"]+results[j,"b7.y"]*all.data$m)
    Y31 <- num3/(1+num3)
    #class 4
    num4 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$m)
    Y41 <- num4/(1+num4)
    
    #P(Y=0|M,X): calculate probability that the outcome (Y) = 0 given mediator (M) and exposure (X, latent classes) 
    Y10 <- 1/(1+num1)
    Y20 <- 1/(1+num2)
    Y30 <- 1/(1+num3)
    Y40 <- 1/(1+num4)
    
    #P(M=1|X): calculate probability that the mediator (M) = 1 given exposure (X, latent classes) 
    #this uses the coefficients from the regression model for the mediator (intercept (b0) and coefs for X (b1-b3))
    #class 1 
    num1 <- exp(results[j,"b0.m"]+results[j,"b1.m"])
    M11 <- num1/(1+num1)
    #class 2
    num2 <- exp(results[j,"b0.m"]+results[j,"b2.m"])
    M21 <- num2/(1+num2) 
    #class 3
    num3 <- exp(results[j,"b0.m"]+results[j,"b3.m"])
    M31 <- num3/(1+num3)   
    #class 4
    num4 <- exp(results[j,"b0.m"])
    M41 <- num4/(1+num4) 
    
    #P(M=0|X): calculate probability that the mediator (M) = 0 given exposure (X, latent classes)    
    M10 <- 1/(1+num1)
    M20 <- 1/(1+num2)
    M30 <- 1/(1+num3)      
    M40 <- 1/(1+num4)
    
    #P(X=x|Y=1,M=1,U)
    #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
    N111<-N1*Y11*M11
    #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
    N211<-N2*Y21*M21
    #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
    N311<-N3*Y31*M31
    #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 1, mediator (M) = 1, and latent class indicators (U)
    N411<-N4*Y41*M41
    
    denom <- N111+N211+N311+N411
    Q111<-N111/denom
    Q211<-N211/denom
    Q311<-N311/denom
    Q411<-N411/denom 
    Q11<-cbind(Q111,Q211,Q311,Q411)
    head(Q11)
    
    #P(X=x|Y=0,M=0,U) 
    #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
    N100<-N1*Y10*M10
    #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
    N200<-N2*Y20*M30
    #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
    N300<-N3*Y30*M30
    #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 0, mediator (M) = 0, and latent class indicators (U)
    N400<-N4*Y40*M40
    
    denom <- N100+N200+N300+N400   
    Q100<-N100/denom
    Q200<-N200/denom
    Q300<-N300/denom
    Q400<-N400/denom    
    Q00<-cbind(Q100,Q200,Q300,Q400)
    head(Q00)
    
    #P(X=x|Y=1,M=0,U)
    #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
    N101<-N1*Y11*M10
    #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
    N201<-N2*Y21*M30
    #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
    N301<-N3*Y31*M30
    #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 1, mediator (M) = 0, and latent class indicators (U)
    N401<-N4*Y41*M40
    
    denom <- N101+N201+N301+N401     
    Q101<-N101/denom
    Q201<-N201/denom
    Q301<-N301/denom
    Q401<-N401/denom    
    Q01<-cbind(Q101,Q201,Q301,Q401)
    head(Q01)    
    
    #P(X=x|Y=0,M=1,U)
    #calcuate probability that exposure (X, latent classes) = 1 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
    N110<-N1*Y10*M11
    #calcuate probability that exposure (X, latent classes) = 2 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
    N210<-N2*Y20*M31
    #calcuate probability that exposure (X, latent classes) = 3 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
    N310<-N3*Y30*M31
    #calcuate probability that exposure (X, latent classes) = 4 given outcome (Y) = 0, mediator (M) = 1, and latent class indicators (U)
    N410<-N4*Y40*M41
    
    denom <- N110+N210+N310+N410    
    Q110<-N110/denom
    Q210<-N210/denom
    Q310<-N310/denom
    Q410<-N410/denom   
    Q10<-cbind(Q110,Q210,Q310,Q410)
    head(Q10)   
    
    #derive the probability of class membership for each individual which takes into account the relationship between the classes, mediator and outcome
    #using each individuals observed data on the mediator and outcome
    Q<-Q11
    #probabilities for those with mediator and outcome absent
    Q[all.data$y==0 & all.data$m==0,]<-Q00[all.data$y==0 & all.data$m==0,]
    #probabilities for those with mediator absent and outcome present
    Q[all.data$y==1 & all.data$m==0,]<-Q01[all.data$y==1 & all.data$m==0,]  
    #probabilities for those with mediator present and outcome absent
    Q[all.data$y==0 & all.data$m==1,]<-Q10[all.data$y==0 & all.data$m==1,]
    #these probabilities will be used to impute class membership for each individual 
    head(Q)
    
    colnames(Q) <-  c("cprob1","cprob2","cprob3","cprob4")
    
    #STEP 3:
    
    #create a matrix to store class membership 
    x <- matrix(NA,sample,latent.classes)
    #now we will use the probabilities of class membership "Q" to randomly assign each individual to a class (X = 1,..k).
    #n=1 (number of random vectors to draw); size=1 per person
    for(k in 1:sample) {x[k,1:latent.classes] <- rmultinom(1,1,Q[k,])}
    
    #add imputed class membership to the data
    all.data[,"x1"] <- x[,1]
    all.data[,"x2"] <- x[,2]
    all.data[,"x3"] <- x[,3]
    all.data[,"x4"] <- x[,4]
    #add in exposure-mediator interactions
    all.data[,"int1"] <- all.data[,"x1"]*all.data[,"m"]
    all.data[,"int2"] <- all.data[,"x2"]*all.data[,"m"]
    all.data[,"int3"] <- all.data[,"x3"]*all.data[,"m"]
    all.data[,"int4"] <- all.data[,"x4"]*all.data[,"m"]
    
    ###################################################
    #we will now perform some checks on the cell sizes from crosstabs for classes by mediator by outcome 
    #we will add these cell sizes into "results" to make traceplots to assess convergence
    
    #create a subset of the data for Y=0 
    no.y <- subset(all.data, all.data$y==0)
    #create a subset of the data for Y=1
    yes.y <- subset(all.data, all.data$y==1)
    
    #crosstabs for x1 and m for those with y=0
    x1m.no.y<-table(no.y$x1,no.y$m)
    
    #this is to make sure matrix is 2 by 2 even when there are zero cells
    if (nrow(x1m.no.y)==1) x1m.no.y <- rbind(x1m.no.y,matrix(0,1,2))
    #crosstabs for x2 and m for those with y=0 
    x2m.no.y<-table(no.y$x2,no.y$m)
    if (nrow(x2m.no.y)==1) x2m.no.y <- rbind(x2m.no.y,matrix(0,1,2))
    #crosstabs for x3 and m for those with y=0 
    x3m.no.y<-table(no.y$x3,no.y$m)
    if (nrow(x3m.no.y)==1) x3m.no.y <- rbind(x3m.no.y,matrix(0,1,2))
    #crosstabs for x4 and m for those with y=0 
    x4m.no.y<-table(no.y$x4,no.y$m)
    if (nrow(x4m.no.y)==1) x4m.no.y <- rbind(x4m.no.y,matrix(0,1,2))
    #crosstabs for x1 and m for those with y=1 
    x1m.yes.y<-table(yes.y$x1,yes.y$m)
    if (nrow(x1m.yes.y)==1) x1m.yes.y <- rbind(x1m.yes.y,matrix(0,1,2))    
    x2m.yes.y<-table(yes.y$x2,yes.y$m)
    if (nrow(x2m.yes.y)==1) x2m.yes.y <- rbind(x2m.yes.y,matrix(0,1,2))    
    x3m.yes.y<-table(yes.y$x3,yes.y$m)
    if (nrow(x3m.yes.y)==1) x3m.yes.y <- rbind(x3m.yes.y,matrix(0,1,2))
    x4m.yes.y<-table(yes.y$x4,yes.y$m)
    if (nrow(x4m.yes.y)==1) x4m.yes.y <- rbind(x4m.yes.y,matrix(0,1,2))
    
    #use empty xmy matrix created earlier and fill in with cell sizes
    #cell size for x=1, m=0, y=0
    xmy[1,"x1m0y0"] <- x1m.no.y[2,1]
    #cell size for x=1, m=1, y=0  
    xmy[1,"x1m1y0"] <- x1m.no.y[2,2]
    #cell size for x=2, m=0, y=0  
    xmy[1,"x2m0y0"] <- x2m.no.y[2,1]
    #cell size for x=2, m=1, y=0  
    xmy[1,"x2m1y0"] <- x2m.no.y[2,2]
    #cell size for x=3, m=0, y=0  
    xmy[1,"x3m0y0"] <- x3m.no.y[2,1]
    #cell size for x=3, m=1, y=0  
    xmy[1,"x3m1y0"] <- x3m.no.y[2,2]
    #cell size for x=4, m=0, y=0  
    xmy[1,"x4m0y0"] <- x4m.no.y[2,1]
    #cell size for x=4, m=1, y=0  
    xmy[1,"x4m1y0"] <- x4m.no.y[2,2]
    #cell size for x=1, m=0, y=1  
    xmy[1,"x1m0y1"] <- x1m.yes.y[2,1]
    #cell size for x=1, m=1, y=1  
    xmy[1,"x1m1y1"] <- x1m.yes.y[2,2]
    #cell size for x=2, m=0, y=1  
    xmy[1,"x2m0y1"] <- x2m.yes.y[2,1]
    #cell size for x=2, m=1, y=1  
    xmy[1,"x2m1y1"] <- x2m.yes.y[2,2]
    #cell size for x=3, m=0, y=1  
    xmy[1,"x3m0y1"] <- x3m.yes.y[2,1]
    #cell size for x=3, m=1, y=1  
    xmy[1,"x3m1y1"] <- x3m.yes.y[2,2]
    #cell size for x=4, m=0, y=1  
    xmy[1,"x4m0y1"] <- x4m.yes.y[2,1]
    #cell size for x=4, m=1, y=1  
    xmy[1,"x4m1y1"] <- x4m.yes.y[2,2]
    
    #flag for number of zero cells in crosstabs for classes by mediator by outcome
    zero.cell[1,1] <- length(which(xmy == 0))
    ###################################################
    
    #fit the firth logistic regression model for P(Y|X,M) to obtain updated parameter estimates now classes have been imputed
    if (zero.cell<2) model.y <- logistf(y~m+x1+x2+x3+int1+int2+int3, data=all.data)
    #if there is more than 1 zero cell, even firth logistic regression does not converge, therefore it is necessary to remove XM interactions from regression model
    if (zero.cell>1) model.y <- logistf(y~m+x1+x2+x3, data=all.data)
    
    #fit the logistic regression model for P(M|X) to obtain updated parameter estimates
    model.m <- logistf(m~x1+x2+x3, data=all.data)
    
    #perturbing the beta coefficients once around coefficients in models using variance-covariance matrix 
    if (zero.cell>1) beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y)),0,0,0)
    if (zero.cell<2) beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y))) 
    beta.m<-c(rmnorm(1,coef(model.m),vcov(model.m)))
    
  } #end of iteration loop
  
  ################################################
  #assessing convergence
  
  #create plot of parameters and cell sizes across each iteration
  results.df <- as.data.frame(results)
  
  #remove the first iteration when no data for X (latent classes) 
  results.df <- subset(results.df, iteration>1)
  
  #to check autocorrelation for each beta after burn in of 100 iterations
  b2.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b2.y"]
  b3.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b3.y"]
  b4.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b4.y"]
  b5.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b5.y"]
  b6.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b6.y"]
  b7.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b7.y"]
  b1.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b1.m"]
  b2.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b2.m"]
  b3.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b3.m"]
  
  #can use a separate folder for updated PCD because of large number of files created (e.g., 80 imputed datasets)
  setwd(upcd.files)
  
  #create a pdf showing traceplots, histograms and autocorrelation plots for each beta, and traceplots and histograms for cell sizes across all iterations
  
  pdf(file=paste0("trace_plots", i, ".pdf")) 
  plot(results.df$iteration, results.df$b2.y,
       xlab = "iteration",
       ylab = "b2.y (y on class 1 vs class 4)")
  lines(results.df$iteration, results.df$b2.y)
  hist (results.df$b2.y, xlab = "b2.y (y on class 1 vs class 4)")
  acf(b2.y, lag.max=50)
  plot(results.df$iteration, results.df$b3.y,
       xlab = "iteration",
       ylab = "b3.y (y on class 2 vs class 4)")
  lines(results.df$iteration, results.df$b3.y)
  hist (results.df$b3.y, xlab = "b3.y (y on class 2 vs class 4)")
  acf(b3.y, lag.max=50)
  plot(results.df$iteration, results.df$b4.y,
       xlab = "iteration",
       ylab = "b4.y (y on class 3 vs class 4)")
  lines(results.df$iteration, results.df$b4.y)
  hist (results.df$b4.y, xlab = "b4.y (y on class 3 vs class 4)")
  acf(b4.y, lag.max=50)
  plot(results.df$iteration, results.df$b5.y,
       xlab = "iteration",
       ylab = "b5.y (y on class 1 x m interaction)")
  lines(results.df$iteration, results.df$b5.y)
  hist (results.df$b5.y, xlab = "b5.y (y on class 1 x m interaction)")
  acf(b5.y, lag.max=50)
  plot(results.df$iteration, results.df$b6.y,
       xlab = "iteration",
       ylab = "b6.y (y on class 2 x m interaction)")
  lines(results.df$iteration, results.df$b6.y)
  hist (results.df$b6.y, xlab = "b6.y (y on class 2 x m interaction)")
  acf(b6.y, lag.max=50)
  plot(results.df$iteration, results.df$b7.y,
       xlab = "iteration",
       ylab = "b7.y (y on class 3 x m interaction)")
  lines(results.df$iteration, results.df$b7.y)
  hist (results.df$b7.y, xlab = "b7.y (y on class 3 x m interaction)")
  acf(b7.y, lag.max=50)
  plot(results.df$iteration, results.df$b1.m,
       xlab = "iteration",
       ylab = "b1.m (m on class 1 vs class 4)")
  lines(results.df$iteration, results.df$b1.m)
  hist (results.df$b1.m, xlab = "b1.m (m on class 1 vs class 4)")
  acf(b1.m, lag.max=50)
  plot(results.df$iteration, results.df$b2.m,
       xlab = "iteration",
       ylab = "b2.m (m on class 2 vs class 4)")
  lines(results.df$iteration, results.df$b2.m)
  hist (results.df$b2.m, xlab = "b2.m (m on class 2 vs class 4)")
  acf(b2.m, lag.max=50)
  plot(results.df$iteration, results.df$b3.m,
       xlab = "iteration",
       ylab = "b3.m (m on class 3 vs class 4")
  lines(results.df$iteration, results.df$b3.m)
  hist (results.df$b3.m, xlab = "b3.m (m on class 3 vs class 4)")
  acf(b3.m, lag.max=50)
  plot(results.df$iteration, results.df$x1.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=1 m=0 y=0")
  lines(results.df$iteration, results.df$x1.m0.y0)
  hist (results.df$x1.m0.y0, xlab = "cell size x=1 m=0 y=0")
  plot(results.df$iteration, results.df$x1.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=1 m=1 y=0")
  lines(results.df$iteration, results.df$x1.m1.y0)
  hist (results.df$x1.m1.y0, xlab = "cell size x=1 m=1 y=0")
  plot(results.df$iteration, results.df$x2.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=2 m=0 y=0")
  lines(results.df$iteration, results.df$x2.m0.y0)
  hist (results.df$x2.m0.y0, xlab = "cell size x=2 m=0 y=0")
  plot(results.df$iteration, results.df$x2.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=2 m=1 y=0")
  lines(results.df$iteration, results.df$x2.m1.y0)
  hist (results.df$x2.m1.y0, xlab = "cell size x=2 m=1 y=0")
  plot(results.df$iteration, results.df$x3.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=3 m=0 y=0")
  lines(results.df$iteration, results.df$x3.m0.y0)
  hist (results.df$x3.m0.y0, xlab = "cell size x=3 m=0 y=0")
  plot(results.df$iteration, results.df$x3.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=3 m=1 y=0")
  lines(results.df$iteration, results.df$x3.m1.y0)
  hist (results.df$x3.m1.y0, xlab = "cell size x=3 m=1 y=0")
  plot(results.df$iteration, results.df$x4.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=4 m=0 y=0")
  lines(results.df$iteration, results.df$x4.m0.y0)
  hist (results.df$x4.m0.y0, xlab = "cell size x=4 m=0 y=0")
  plot(results.df$iteration, results.df$x4.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=4 m=1 y=0")
  lines(results.df$iteration, results.df$x4.m1.y0)
  hist (results.df$x4.m1.y0, xlab = "cell size x=4 m=1 y=0")
  plot(results.df$iteration, results.df$x1.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=1 m=0 y=1")
  lines(results.df$iteration, results.df$x1.m0.y1)
  hist (results.df$x1.m0.y1, xlab = "cell size x=1 m=0 y=1")
  plot(results.df$iteration, results.df$x1.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=1 m=1 y=1")
  lines(results.df$iteration, results.df$x1.m1.y1)
  hist (results.df$x1.m1.y1, xlab = "cell size x=1 m=1 y=1")
  plot(results.df$iteration, results.df$x2.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=2 m=0 y=1")
  lines(results.df$iteration, results.df$x2.m0.y1)
  hist (results.df$x2.m0.y1, xlab = "cell size x=2 m=0 y=1")
  plot(results.df$iteration, results.df$x2.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=2 m=1 y=1")
  lines(results.df$iteration, results.df$x2.m1.y1)
  hist (results.df$x2.m1.y1, xlab = "cell size x=2 m=1 y=1")
  plot(results.df$iteration, results.df$x3.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=3 m=0 y=1")
  lines(results.df$iteration, results.df$x3.m0.y1)
  hist (results.df$x3.m0.y1, xlab = "cell size x=3 m=0 y=1")
  plot(results.df$iteration, results.df$x3.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=3 m=1 y=1")
  lines(results.df$iteration, results.df$x3.m1.y1)
  hist (results.df$x3.m1.y1, xlab = "cell size x=3 m=1 y=1")
  plot(results.df$iteration, results.df$x4.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=4 m=0 y=1")
  lines(results.df$iteration, results.df$x4.m0.y1)
  hist (results.df$x4.m0.y1, xlab = "cell size x=4 m=0 y=1")
  plot(results.df$iteration, results.df$x4.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=4 m=1 y=1")
  lines(results.df$iteration, results.df$x4.m1.y1)
  hist (results.df$x4.m1.y1, xlab = "cell size x=4 m=1 y=1")
  dev.off()
  
  #analysis
  
  #combine imputed class membership stored in "imp.n" with original data
  #prepare 1 mplus .dat file for each imputed dataset to run mediation model (80 .dat files should be created)
  for(l in 1:imp.n) {
    imp.subset <- cbind(data.original,subset(imp, imp[,1]==l))
    prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))}
  
  #create "imp.txt" file for mplus to call imputed datasets
  imp.txt <- matrix(NA,imp.n,1)
  for(l in 1:imp.n) {
    imp.txt[l,1] <- paste0("imp_", l, ".dat")}
  write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #run mediation model using updated PCD
  runModels("2i pcd mediation.inp")
  
  #read in parameters from updated PCD mediation model 
  est.pcd <- readModels("2i pcd mediation.out", what="parameters")$parameters$`unstandardized`
  #store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "upcd.all"
  pcd.all[i,(match("tot.1v4",colnames(pcd.all))):(match("pde.4v1",colnames(pcd.all)))] <- est.pcd[100:135,"est"] # TOT_1V4 to PDE_4V1
  pcd.all[i,(match("tot.1v4.se",colnames(pcd.all))):(match("pde.4v1.se",colnames(pcd.all)))] <- est.pcd[100:135,"se"] # TOT_1V4 to PDE_4V1
  pcd.all[i,(match("p1",colnames(pcd.all))):(match("p4",colnames(pcd.all)))] <- est.pcd[72:75,"est"]  # P_X1 to P_X4
  
  #flags for potential issues: "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp"
  
  #flag if fixed threshold in unconditional model which results in missing value in tech3 covariance matrix of parameters
  pcd.all[i,"fixed.th"] <- 0 #set flag to 0 to start
  for(k in 1:(indicators*latent.classes+(latent.classes-1))*(indicators*latent.classes+(latent.classes-1))) { #size of tech3 matrix is 23x23 (529) 
    if (tech3[k]==999) pcd.all[i,"fixed.th"] <- 1}
  
  #when perturbing within-class thresholds, those with a large standard error can go out of bounds (e.g., corresponding to a probability that is not between 0 and 100%)
  #we created a flag for this in "results" so here we will move this into "pcd.all" and record whether this was the case in any iteration
  pcd.all[i,"large.th"] <- 0 #set flag to 0 to start
  results[1,"large.th"]<-0 
  for(k in 1:(cycles*imp.n+burnin)) {
    if (results[k,"large.th"]==1) pcd.all[i,"large.th"] <- 1}
  
  #we also want a flag so that we know the largest standard error for a threshold in the unconditional model (and the threshold this SE corresponds to)
  #this flag captures the largest SE
  pcd.all[i,"largest.th.se"] <- max(coef.se.original[,"se"])
  #this flag records the threshold that corresponds to largest SE
  pcd.all[i,"largest.th"] <- coef.se.original[which.max(coef.se.original[,"se"]),"est"]
  
  #we created a flag for number of zero cells in the crosstabs for classes by mediator by outcome in "results"
  #we will move this flag to "pcd.all" and record the largest number of zero cells in any iteration
  results[1,"zero.cell"]<-0 #change row 1 in "results" to 0 as latent classes not imputed in first iteration
  pcd.all[i,"zero.cell"] <- max(results[,"zero.cell"])
  
  #we will also flag if there was a zero cell (in the crosstabs for classes by mediator by outcome) in any of the iterations that we saved as one of the 80 imputed datasets
  zero.cell.imp <- matrix(NA,1,imp.n)
  for(l in 1:imp.n) {
    zero.cell.imp[1,l] <- results[cycles*l+burnin,"zero.cell"]}
  
  pcd.all[i,"zero.cell.imp"] <- 0 #set flag to 0 to start
  for(k in 1:imp.n) {
    if (zero.cell.imp[1,k]>0) pcd.all[1,"zero.cell.imp"] <- 1}
  
  cat(i,"\n")} #end of simulations loop

#############################################################

#saving results for comparisons of interest (e.g. with low class as the reference group) and exporting out mediation effects for each method

#onestep
onestep <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes+1)
colnames(onestep) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                       "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4","largest.th.se")

#move over class probabilities
onestep[,"p1"]<-onestep.all[,"p1"]
onestep[,"p2"]<-onestep.all[,"p2"]
onestep[,"p3"]<-onestep.all[,"p3"]
onestep[,"p4"]<-onestep.all[,"p4"]
onestep[,"largest.th.se"]<-onestep.all[,"largest.th.se"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory parameters (auc1-auc4))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc1"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc1"]>onestep.all[i,"auc4"])
    onestep[i,1:18] <- onestep.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc1"]>onestep.all[i,"auc3"])
    onestep[i,1:18] <- onestep.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc1"]>onestep.all[i,"auc2"])
    onestep[i,1:18] <- onestep.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"] & onestep.all[i,"auc3"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc4"] & onestep.all[i,"auc4"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc2"]>onestep.all[i,"auc1"] & onestep.all[i,"auc2"]>onestep.all[i,"auc3"] & onestep.all[i,"auc3"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc3"]>onestep.all[i,"auc4"]
     & onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc2"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(onestep.all[i,"auc4"]>onestep.all[i,"auc1"] & onestep.all[i,"auc4"]>onestep.all[i,"auc2"] & onestep.all[i,"auc4"]>onestep.all[i,"auc3"]
     & onestep.all[i,"auc3"]>onestep.all[i,"auc1"] & onestep.all[i,"auc3"]>onestep.all[i,"auc2"] & onestep.all[i,"auc2"]>onestep.all[i,"auc1"])
    onestep[i,1:18] <- onestep.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#bch
bch <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(bch) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                   "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
bch[,"p1"]<-bch.all[,"p1"]
bch[,"p2"]<-bch.all[,"p2"]
bch[,"p3"]<-bch.all[,"p3"]
bch[,"p4"]<-bch.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    bch[i,1:18] <- bch.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    bch[i,1:18] <- bch.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    bch[i,1:18] <- bch.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    bch[i,1:18] <- bch.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}


#modal
modal <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(modal) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                     "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
modal[,"p1"]<-modal.all[,"p1"]
modal[,"p2"]<-modal.all[,"p2"]
modal[,"p3"]<-modal.all[,"p3"]
modal[,"p4"]<-modal.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    modal[i,1:18] <- modal.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    modal[i,1:18] <- modal.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    modal[i,1:18] <- modal.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    modal[i,1:18] <- modal.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#npcd
npcd <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(npcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                    "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
npcd[,"p1"]<-npcd.all[,"p1"]
npcd[,"p2"]<-npcd.all[,"p2"]
npcd[,"p3"]<-npcd.all[,"p3"]
npcd[,"p4"]<-npcd.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    npcd[i,1:18] <- npcd.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    npcd[i,1:18] <- npcd.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    npcd[i,1:18] <- npcd.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    npcd[i,1:18] <- npcd.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#incpcd
incpcd <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes)
colnames(incpcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                      "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#move over class probabilities
incpcd[,"p1"]<-incpcd.all[,"p1"]
incpcd[,"p2"]<-incpcd.all[,"p2"]
incpcd[,"p3"]<-incpcd.all[,"p3"]
incpcd[,"p4"]<-incpcd.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low) 
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc1"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc1"]>auc.inc[i,"auc4"])
    incpcd[i,1:18] <- incpcd.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc1"]>auc.inc[i,"auc3"])
    incpcd[i,1:18] <- incpcd.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc1"]>auc.inc[i,"auc2"])
    incpcd[i,1:18] <- incpcd.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"] & auc.inc[i,"auc3"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc4"] & auc.inc[i,"auc4"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc2"]>auc.inc[i,"auc1"] & auc.inc[i,"auc2"]>auc.inc[i,"auc3"] & auc.inc[i,"auc3"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc3"]>auc.inc[i,"auc4"]
     & auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc2"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc.inc[i,"auc4"]>auc.inc[i,"auc1"] & auc.inc[i,"auc4"]>auc.inc[i,"auc2"] & auc.inc[i,"auc4"]>auc.inc[i,"auc3"]
     & auc.inc[i,"auc3"]>auc.inc[i,"auc1"] & auc.inc[i,"auc3"]>auc.inc[i,"auc2"] & auc.inc[i,"auc2"]>auc.inc[i,"auc1"])
    incpcd[i,1:18] <- incpcd.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#pcd
pcd <- matrix(NA,sims,(latent.classes-1)*mediation.effects+latent.classes+6)
colnames(pcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                   "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4",
                   "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp")

#move over class probabilities
pcd[,"p1"]<-pcd.all[,"p1"]
pcd[,"p2"]<-pcd.all[,"p2"]
pcd[,"p3"]<-pcd.all[,"p3"]
pcd[,"p4"]<-pcd.all[,"p4"]
pcd[,"fixed.th"]<-pcd.all[,"fixed.th"]
pcd[,"large.th"]<-pcd.all[,"large.th"]
pcd[,"largest.th"]<-pcd.all[,"largest.th"]
pcd[,"largest.th.se"]<-pcd.all[,"largest.th.se"]
pcd[,"zero.cell"]<-pcd.all[,"zero.cell"]
pcd[,"zero.cell.imp"]<-pcd.all[,"zero.cell.imp"]

#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
for(i in 1:sims) {
  #1234 (eop,ao,cl,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc3"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc2"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc1"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc1"]>auc[i,"auc4"])
    pcd[i,1:18] <- pcd.all[i,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc1"]>auc[i,"auc3"])
    pcd[i,1:18] <- pcd.all[i,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc1"]>auc[i,"auc2"])
    pcd[i,1:18] <- pcd.all[i,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc2"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc4"] & auc[i,"auc4"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc2"]>auc[i,"auc1"] & auc[i,"auc2"]>auc[i,"auc3"] & auc[i,"auc3"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc3"]>auc[i,"auc4"]
     & auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[i,"auc4"]>auc[i,"auc1"] & auc[i,"auc4"]>auc[i,"auc2"] & auc[i,"auc4"]>auc[i,"auc3"]
     & auc[i,"auc3"]>auc[i,"auc1"] & auc[i,"auc3"]>auc[i,"auc2"] & auc[i,"auc2"]>auc[i,"auc1"])
    pcd[i,1:18] <- pcd.all[i,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]}

#put results in a table and export to a results folder
setwd(results.files)

write.table(onestep, file="table.onestep.txt", sep = ",")
write.table(bch, file="table.bch.txt", sep = ",")
write.table(modal, file="table.modal.txt", sep = ",")
write.table(npcd, file="table.npcd.txt", sep = ",")
write.table(incpcd, file="table.incpcd.txt", sep = ",")
write.table(pcd, file="table.pcd.txt", sep = ",")

write.table(entropy, file="table.entropy.txt", sep = ",")
write.table(auc, file="table.auc.txt", sep = ",")

################################################end of script#####################################################################

