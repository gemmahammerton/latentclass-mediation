#R script for comparison of methods to relate a latent class exposure to distal outcomes within a counterfactual mediation model
#methods: one-step, bch, modal class, non-inclusive PCD, inclusive PCD and updated PCD (new method)

#additional files needed to run this script:
#mplus data file: applied.dat (not publicly available)
#mplus input files: 
#"4b uncond latent class.inp"
#"4c onestep mediation.inp"
#"4d bch mediation.inp"
#"4e modal mediation.inp"
#"4f npcd mediation.inp"
#"4g inc latent class.inp"
#"4h incpcd mediation.inp"
#"4i pcd mediation.inp"

#using Mplus v8.8 (display order of results in Mplus output files can depend on version)

############### INSTALL AND LOAD PACKAGES ###################

#install the required packages (if not already installed)
list.of.packages <- c("MplusAutomation", "logistf", "mnormt", "gdata", "nnet", "brglm2", "detectseparation", "ggplot2", "forcats")
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
getNamespaceVersion("nnet") #7.3-16
getNamespaceVersion("brglm2") #0.8.2
getNamespaceVersion("detectseparation") #0.2
getNamespaceVersion("ggplot2") 
getNamespaceVersion("forcats") 

############### DEFINE DIRECTORIES ##########################

#Folder where applied data and mplus input files (uncond, onestep, bch, modal and inc) are saved:
sim.data <- "file location/applied example"
#*Folder where non-inclusive PCD files are saved:
npcd.files <- "file location/applied example"
#*Folder where inclusive PCD files are saved:
incpcd.files <- "file location/applied example"
#*Folder where updated PCD files are saved:
upcd.files <- "file location/applied example"
#*Folder where results are saved:
results.files <- "file location/applied example"

#script for applied example examining whether association between conduct trajectories and internalising problems is mediated by drug use (using ALSPAC data)
#############################################################
#X = 4 class latent nominal exposure (conduct trajectories: early onset persistent, adolescent onset, childhood limited and low)
latent.classes <- 4
#M = binary mediator (illicit drug use)
#Y = binary outcome (internalising problems)
#U1-U6 = binary latent class indicators (conduct problems from age 4 to 13 years)
indicators <- 6
#confounders (sex and cumulative sociodemographic risk score)
confounders <- 2

############### PREPARE MATRICES ##########################

#we wish to compare each latent class with all the others 
class.comp <- 4*3
#we wish to report 3 mediation effects (total, indirect - tie, and direct - pde) and their standard errors
mediation.effects <- 3*2

#create matrices to store the mediation results for each method
#number of rows = 1
#number of columns is based on number of class comparisons (12), number of mediation effects (total, tie and pde) and SEs (6) and 
#number of extra parameters (e.g. area under the trajectory, class probabilities, or additional checks)

onestep.all <- matrix(NA,1,(class.comp*mediation.effects+latent.classes*2))

#labels represent the mediation effects (total, tie and pde) for each class comparison, area under the trajectory for each class, and class probabilities
colnames(onestep.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                           "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                           "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                           "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                           "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se", 
                           "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se", 
                           "auc1","auc2","auc3","auc4","p1","p2","p3","p4")

onestep <- matrix(NA,1,(latent.classes-1)*mediation.effects+latent.classes)
colnames(onestep) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                       "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

bch.all <- matrix(NA,1,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(bch.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                       "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                       "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                       "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                       "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                       "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                       "p1","p2","p3","p4")

bch <- matrix(NA,1,(latent.classes-1)*mediation.effects+latent.classes)
colnames(bch) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                       "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

modal.all <- matrix(NA,1,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(modal.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                         "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                         "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                         "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                         "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                         "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                         "p1","p2","p3","p4")

modal <- matrix(NA,1,(latent.classes-1)*mediation.effects+latent.classes)
colnames(modal) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                       "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

npcd.all <- matrix(NA,1,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(npcd.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                        "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                        "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                        "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                        "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                        "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                        "p1","p2","p3","p4")

npcd <- matrix(NA,1,(latent.classes-1)*mediation.effects+latent.classes)
colnames(npcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                    "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

incpcd.all <- matrix(NA,1,(class.comp*mediation.effects+latent.classes))

#labels represent the mediation effects (total, tie and pde) for each class comparison and class probabilities  
colnames(incpcd.all) <- c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                          "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                          "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                          "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                          "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                          "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se",
                          "p1","p2","p3","p4")

incpcd <- matrix(NA,1,(latent.classes-1)*mediation.effects+latent.classes)
colnames(incpcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                    "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

pcd.all <- matrix(NA,1,(class.comp*mediation.effects+latent.classes+6))

#labels represent the mediation effects (total, tie and pde) for each class comparison, 6 flags for potential issues and class probabilities 
colnames(pcd.all) <-  c("tot.1v4","tot.1v3","tot.1v2","tot.2v4","tot.2v3","tot.2v1","tot.3v4","tot.3v2","tot.3v1","tot.4v3","tot.4v2","tot.4v1",
                        "tie.1v4","tie.1v3","tie.1v2","tie.2v4","tie.2v3","tie.2v1","tie.3v4","tie.3v2","tie.3v1","tie.4v3","tie.4v2","tie.4v1",
                        "pde.1v4","pde.1v3","pde.1v2","pde.2v4","pde.2v3","pde.2v1","pde.3v4","pde.3v2","pde.3v1","pde.4v3","pde.4v2","pde.4v1",
                        "tot.1v4.se","tot.1v3.se","tot.1v2.se","tot.2v4.se","tot.2v3.se","tot.2v1.se","tot.3v4.se","tot.3v2.se","tot.3v1.se","tot.4v3.se","tot.4v2.se","tot.4v1.se",
                        "tie.1v4.se","tie.1v3.se","tie.1v2.se","tie.2v4.se","tie.2v3.se","tie.2v1.se","tie.3v4.se","tie.3v2.se","tie.3v1.se","tie.4v3.se","tie.4v2.se","tie.4v1.se",
                        "pde.1v4.se","pde.1v3.se","pde.1v2.se","pde.2v4.se","pde.2v3.se","pde.2v1.se","pde.3v4.se","pde.3v2.se","pde.3v1.se","pde.4v3.se","pde.4v2.se","pde.4v1.se", "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp",
                        "p1","p2","p3","p4")

pcd <- matrix(NA,1,(latent.classes-1)*mediation.effects+latent.classes)
colnames(pcd) <- c("tot.eop","tot.ao","tot.cl","tie.eop","tie.ao","tie.cl","pde.eop","pde.ao","pde.cl",
                    "tot.eop.se","tot.ao.se","tot.cl.se","tie.eop.se","tie.ao.se","tie.cl.se","pde.eop.se","pde.ao.se","pde.cl.se","p1","p2","p3","p4")

#create matrix to store the area under the trajectory for each class in the unconditional model (this allows us to know the order of the classes)
auc <- matrix(NA,1,latent.classes)
colnames(auc) <- c("auc1","auc2","auc3","auc4")

#create matrix to store the area under the trajectory for each class in the inclusive latent class model (this allows us to know the order of the classes)
auc.inc <- matrix(NA,1,latent.classes)
colnames(auc.inc) <- c("auc1","auc2","auc3","auc4")

############### READ IN THE DATA ##########################

#set working directory using file paths saved at the start
setwd(sim.data)

#read in data for applied example from mplus .dat file
data.original <- read.csv(file="applied.dat", sep=",", na.strings="-9999", header=FALSE, col.names = c("aln", #identification number
                                                                                                       "female", #confounder
                                                                                                       "rs_soc", #confounder
                                                                                                       "cd1","cd2","cd3","cd4","cd5","cd6", #latent class indicators
                                                                                                       "cd_modclass", #modal class assignment for the exposure
                                                                                                       "drugs_18", #mediator
                                                                                                       "int_18", #outcome
                                                                                                       "qlet" #identification letter
                                                                                                       ))
  
#subset data to include only those with complete data on outcome, mediator, confounders and at least one latent class indicator
data.cc <- subset(data.original, is.na(female)==FALSE & is.na(rs_soc)==FALSE & is.na(drugs_18)==FALSE & is.na(int_18)==FALSE & 
                      (is.na(cd1)==FALSE | is.na(cd2)==FALSE | is.na(cd3)==FALSE | is.na(cd4)==FALSE | is.na(cd5)==FALSE | is.na(cd6)==FALSE))

#data restricted to 3039
dim<-dim(data.cc)
#create an indicator of sample size
sample<-dim[1]

#replicate the original data and add empty columns to store imputed latent class membership for 4 classes (once generated) 
#and interactions between each latent class and the mediator
all.data <- data.cc
all.data$x1 <- NA #empty column to add dummy code for membership in latent class 1
all.data$x2 <- NA #empty column to add dummy code for membership in latent class 2
all.data$x3 <- NA #empty column to add dummy code for membership in latent class 3
all.data$x4 <- NA #empty column to add dummy code for membership in latent class 4
all.data$int1 <- NA  #empty column to add interaction between latent class 1 and mediator
all.data$int2 <- NA  #empty column to add interaction between latent class 2 and mediator
all.data$int3 <- NA  #empty column to add interaction between latent class 3 and mediator
all.data$int4 <- NA  #empty column to add interaction between latent class 4 and mediator
all.data$x <- NA  #empty column to add interaction between latent class 4 and mediator

############### UNCONDITIONAL LATENT CLASS MODEL ##########################

#set working directory using file paths saved at the start
setwd(sim.data)

#run the unconditional latent class model using sample with complete data on outcome, mediator and confounders (n=3,039)
runModels("4b uncond latent class.inp")

#read in the mplus output file and name it "model_output" to use later
model_output <- readModels("4b uncond latent class.out")

#save parameters from unconditional latent class model (within-class thresholds for latent class indicators and class intercepts) and their SEs
coef.se.original <- model_output$parameters$`unstandardized`[1:(indicators*latent.classes+(latent.classes-1)),3:4] 
#labels represent six within-class thresholds across 4 classes and 3 class intercepts)
rownames(coef.se.original) <- c("th1.c1","th2.c1","th3.c1","th4.c1","th5.c1","th6.c1", #within-class thresholds for class 1
                                "th1.c2","th2.c2","th3.c2","th4.c2","th5.c2","th6.c2", #within-class thresholds for class 2
                                "th1.c3","th2.c3","th3.c3","th4.c3","th5.c3","th6.c3", #within-class thresholds for class 3
                                "th1.c4","th2.c4","th3.c4","th4.c4","th5.c4","th6.c4", #within-class thresholds for class 4
                                "int.c1","int.c2","int.c3") #latent class intercepts

#record the class separation (entropy = 0.712)
entropy <- model_output$summaries$`Entropy`

#create a matrix of class probabilities (low = 65%, CL = 21%, EOP = 9%, AO = 4%)
p.original <- matrix(NA,latent.classes,1)
#these can be calculated using the 3 class intercepts we saved in "coef.se.original"
p.original[1,1] <- exp(coef.se.original["int.c1","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
p.original[2,1] <- exp(coef.se.original["int.c2","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))
p.original[3,1] <- exp(coef.se.original["int.c3","est"])/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))  
p.original[4,1] <- 1/(1 + exp(coef.se.original["int.c1","est"]) + exp(coef.se.original["int.c2","est"]) + exp(coef.se.original["int.c3","est"]))

#derive area under the trajectory parameters so we know the order of the latent classes in the unconditional model
#these can be calculated using the within-class thresholds we saved in "coef.se.original"
#from the auc we can see that the order of the classes is: low = 0.59, CL = 4.39, EOP = 16.14, AO = 10.96
auc[1,"auc1"] <- exp(-1*coef.se.original["th1.c1","est"])/(1+exp(-1*coef.se.original["th1.c1","est"]))+2*exp(-1*coef.se.original["th2.c1","est"])/(1+exp(-1*coef.se.original["th2.c1","est"]))+3*exp(-1*coef.se.original["th3.c1","est"])/(1+exp(-1*coef.se.original["th3.c1","est"]))+4*exp(-1*coef.se.original["th4.c1","est"])/(1+exp(-1*coef.se.original["th4.c1","est"]))+5*exp(-1*coef.se.original["th5.c1","est"])/(1+exp(-1*coef.se.original["th5.c1","est"]))+6*exp(-1*coef.se.original["th6.c1","est"])/(1+exp(-1*coef.se.original["th6.c1","est"]))
auc[1,"auc2"] <- exp(-1*coef.se.original["th1.c2","est"])/(1+exp(-1*coef.se.original["th1.c2","est"]))+2*exp(-1*coef.se.original["th2.c2","est"])/(1+exp(-1*coef.se.original["th2.c2","est"]))+3*exp(-1*coef.se.original["th3.c2","est"])/(1+exp(-1*coef.se.original["th3.c2","est"]))+4*exp(-1*coef.se.original["th4.c2","est"])/(1+exp(-1*coef.se.original["th4.c2","est"]))+5*exp(-1*coef.se.original["th5.c2","est"])/(1+exp(-1*coef.se.original["th5.c2","est"]))+6*exp(-1*coef.se.original["th6.c2","est"])/(1+exp(-1*coef.se.original["th6.c2","est"]))
auc[1,"auc3"] <- exp(-1*coef.se.original["th1.c3","est"])/(1+exp(-1*coef.se.original["th1.c3","est"]))+2*exp(-1*coef.se.original["th2.c3","est"])/(1+exp(-1*coef.se.original["th2.c3","est"]))+3*exp(-1*coef.se.original["th3.c3","est"])/(1+exp(-1*coef.se.original["th3.c3","est"]))+4*exp(-1*coef.se.original["th4.c3","est"])/(1+exp(-1*coef.se.original["th4.c3","est"]))+5*exp(-1*coef.se.original["th5.c3","est"])/(1+exp(-1*coef.se.original["th5.c3","est"]))+6*exp(-1*coef.se.original["th6.c3","est"])/(1+exp(-1*coef.se.original["th6.c3","est"]))
auc[1,"auc4"] <- exp(-1*coef.se.original["th1.c4","est"])/(1+exp(-1*coef.se.original["th1.c4","est"]))+2*exp(-1*coef.se.original["th2.c4","est"])/(1+exp(-1*coef.se.original["th2.c4","est"]))+3*exp(-1*coef.se.original["th3.c4","est"])/(1+exp(-1*coef.se.original["th3.c4","est"]))+4*exp(-1*coef.se.original["th4.c4","est"])/(1+exp(-1*coef.se.original["th4.c4","est"]))+5*exp(-1*coef.se.original["th5.c4","est"])/(1+exp(-1*coef.se.original["th5.c4","est"]))+6*exp(-1*coef.se.original["th6.c4","est"])/(1+exp(-1*coef.se.original["th6.c4","est"]))

#import tech3 as a covariance matrix to allow us to use the covariance between parameters for perturbing in updated PCD approach
#############################################################
#this outputs paramCov as a matrix object in R 
tech3 <- model_output$tech3$paramCov
#this creates a full, symmetrical covariance matrix
upperTriangle(tech3) <- lowerTriangle(tech3, byrow=TRUE)
cov <- tech3

#address any fixed parameters: within class thresholds that have been fixed at 15 or -15 (representing a probability of 0 or 100%) do not have (co)variances
#the code below means that fixed parameters will still get perturbed in updated PCD but only a very small amount
#in this applied example there are no fixed parameters, therefore this section of code is not needed
#replace missing (999) in covariance matrix with 0 to represent no covariance for the fixed parameters
cov[cov==999] <- 0
#change the variance for fixed parameters to be very small (0.000000001)
for(j in 1:(indicators*latent.classes+(latent.classes-1))) {
  if (cov[j,j]==0) cov[j,j]<- 0.000000001}

###################### ONESTEP MODEL ##########################

#run the onestep latent class mediation model
runModels("4c onestep mediation.inp")

#read in parameters from one-step latent class model    
est.onestep <- readModels("4c onestep mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects, area under the trajectory for each class, and class probabilities) and their SEs in the matrix created at the start "onestep.all"
onestep.all[1, (match("tot.1v4",colnames(onestep.all))):(match("pde.4v1",colnames(onestep.all)))] <- est.onestep[126:161,"est"] # TOT_1V4 to PDE_4V1
onestep.all[1, (match("tot.1v4.se",colnames(onestep.all))):(match("pde.4v1.se",colnames(onestep.all)))] <- est.onestep[126:161,"se"] # TOT_1V4 to PDE_4V1
onestep.all[1, (match("auc1",colnames(onestep.all))):(match("auc4",colnames(onestep.all)))] <- est.onestep[162:165,"est"] # AUC
onestep.all[1, (match("p1",colnames(onestep.all))):(match("p4",colnames(onestep.all)))] <- est.onestep[98:101,"est"]  # P_X1 to P_X4

#saving results for comparisons of interest (e.g. with low class as the reference group)

#move over class probabilities
onestep[,"p1"]<-onestep.all[,"p1"]
onestep[,"p2"]<-onestep.all[,"p2"]
onestep[,"p3"]<-onestep.all[,"p3"]
onestep[,"p4"]<-onestep.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in mplus output - determined using the area under the trajectory parameters (auc1-auc4))
  #1234 (eop,ao,cl,low)
  if(onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"])
    onestep[1,1:18] <- onestep.all[1,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"])
    onestep[1,1:18] <- onestep.all[1,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"])
    onestep[1,1:18] <- onestep.all[1,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc4"]>onestep.all[1,"auc3"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"])
    onestep[1,1:18] <- onestep.all[1,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc3"]>onestep.all[1,"auc4"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"])
    onestep[1,1:18] <- onestep.all[1,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc4"]>onestep.all[1,"auc3"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"])
    onestep[1,1:18] <- onestep.all[1,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"])
    onestep[1,1:18] <- onestep.all[1,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"])
    onestep[1,1:18] <- onestep.all[1,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"])
    onestep[1,1:18] <- onestep.all[1,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"]
     & onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"])
    onestep[1,1:18] <- onestep.all[1,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"])
    onestep[1,1:18] <- onestep.all[1,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"]
     & onestep.all[1,"auc1"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"])
    onestep[1,1:18] <- onestep.all[1,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"])
    onestep[1,1:18] <- onestep.all[1,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"])
    onestep[1,1:18] <- onestep.all[1,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"] & onestep.all[1,"auc1"]>onestep.all[1,"auc4"])
    onestep[1,1:18] <- onestep.all[1,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"]
     & onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc1"]>onestep.all[1,"auc3"])
    onestep[1,1:18] <- onestep.all[1,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc2"])
    onestep[1,1:18] <- onestep.all[1,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"]
     & onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc1"]>onestep.all[1,"auc2"])
    onestep[1,1:18] <- onestep.all[1,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"] & onestep.all[1,"auc4"]>onestep.all[1,"auc1"])
    onestep[1,1:18] <- onestep.all[1,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"] & onestep.all[1,"auc3"]>onestep.all[1,"auc1"])
    onestep[1,1:18] <- onestep.all[1,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc4"] & onestep.all[1,"auc4"]>onestep.all[1,"auc1"])
    onestep[1,1:18] <- onestep.all[1,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"]
     & onestep.all[1,"auc2"]>onestep.all[1,"auc1"] & onestep.all[1,"auc2"]>onestep.all[1,"auc3"] & onestep.all[1,"auc3"]>onestep.all[1,"auc1"])
    onestep[1,1:18] <- onestep.all[1,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc3"]>onestep.all[1,"auc4"]
     & onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc2"]>onestep.all[1,"auc1"])
    onestep[1,1:18] <- onestep.all[1,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(onestep.all[1,"auc4"]>onestep.all[1,"auc1"] & onestep.all[1,"auc4"]>onestep.all[1,"auc2"] & onestep.all[1,"auc4"]>onestep.all[1,"auc3"]
     & onestep.all[1,"auc3"]>onestep.all[1,"auc1"] & onestep.all[1,"auc3"]>onestep.all[1,"auc2"] & onestep.all[1,"auc2"]>onestep.all[1,"auc1"])
    onestep[1,1:18] <- onestep.all[1,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]

######################## BCH ##########################

#bch can only be used with multiple latent classes if they are combined into one class
#use one 8 class model for X and M and set up mediation model using nom-nom-cat approach (because bch needs to include an auxiliary variable)

#read in bch weights and modal class assignment which were exported out of the unconditional latent class model in "bch.txt"
data.bch <- read.table("bch.txt", sep="", na.strings="*", header=FALSE)
#select the weights for each class and the modal class assignment
data.bch <- data.bch[,c(7:10,15)]
#labels representing one weight for each class and the modal class assignment
colnames(data.bch) <- c("bch1","bch2","bch3","bch4","modal")
#combine the bch weights with the mediator and outcome from the original data
data.bch <- cbind(data.cc[,c("female","rs_soc","drugs_18","int_18")], data.bch)

#create weights for 8 class model by multiplying the bch weight for each class with the observed data on the mediator
data.bch$bch11 <- data.bch$bch1*data.bch$drugs_18
data.bch$bch12 <- data.bch$bch1*(1-data.bch$drugs_18)
data.bch$bch21 <- data.bch$bch2*data.bch$drugs_18
data.bch$bch22 <- data.bch$bch2*(1-data.bch$drugs_18)
data.bch$bch31 <- data.bch$bch3*data.bch$drugs_18
data.bch$bch32 <- data.bch$bch3*(1-data.bch$drugs_18)
data.bch$bch41 <- data.bch$bch4*data.bch$drugs_18
data.bch$bch42 <- data.bch$bch4*(1-data.bch$drugs_18)

#prepare mplus .dat file using "data.bch"  
prepareMplusData(data.bch,"bch.dat")
#run mediation model using bch method
runModels("4d bch mediation.inp")

#read in parameters from bch latent class mediation model   
est.bch <- readModels("4d bch mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "bch.all"
bch.all[1,(match("tot.1v4",colnames(bch.all))):(match("pde.4v1",colnames(bch.all)))] <- est.bch[70:105,"est"] # TOT_1V4 to PDE_4V1
bch.all[1,(match("tot.1v4.se",colnames(bch.all))):(match("pde.4v1.se",colnames(bch.all)))] <- est.bch[70:105,"se"] # TOT_1V4 to PDE_4V1
bch.all[1,(match("p1",colnames(bch.all))):(match("p4",colnames(bch.all)))] <- est.bch[42:45,"est"]  # P_X1 to P_X4

#move over class probabilities
bch[,"p1"]<-bch.all[,"p1"]
bch[,"p2"]<-bch.all[,"p2"]
bch[,"p3"]<-bch.all[,"p3"]
bch[,"p4"]<-bch.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in mplus output - determined using the area under the trajectory matrix (auc))
  #1234 (eop,ao,cl,low)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc4"])
    bch[1,1:18] <- bch.all[1,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc3"])
    bch[1,1:18] <- bch.all[1,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc2"]>auc[1,"auc4"])
    bch[1,1:18] <- bch.all[1,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc3"])
    bch[1,1:18] <- bch.all[1,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc2"])
    bch[1,1:18] <- bch.all[1,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc2"])
    bch[1,1:18] <- bch.all[1,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc4"])
    bch[1,1:18] <- bch.all[1,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc3"])
    bch[1,1:18] <- bch.all[1,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc2"]>auc[1,"auc4"])
    bch[1,1:18] <- bch.all[1,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc3"])
    bch[1,1:18] <- bch.all[1,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc2"])
    bch[1,1:18] <- bch.all[1,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc2"])
    bch[1,1:18] <- bch.all[1,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc1"]>auc[1,"auc4"])
    bch[1,1:18] <- bch.all[1,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc3"])
    bch[1,1:18] <- bch.all[1,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc1"]>auc[1,"auc4"])
    bch[1,1:18] <- bch.all[1,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc3"])
    bch[1,1:18] <- bch.all[1,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc2"])
    bch[1,1:18] <- bch.all[1,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc2"])
    bch[1,1:18] <- bch.all[1,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc1"])
    bch[1,1:18] <- bch.all[1,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc1"])
    bch[1,1:18] <- bch.all[1,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc1"])
    bch[1,1:18] <- bch.all[1,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc1"])
    bch[1,1:18] <- bch.all[1,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc1"])
    bch[1,1:18] <- bch.all[1,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc1"])
    bch[1,1:18] <- bch.all[1,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]


######################## MODAL CLASS ASSIGNMENT ##########################

#run mediation model using modal class assignment
runModels("4e modal mediation.inp")

#read in parameters from modal class mediation model  
est.modal <- readModels("4e modal mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "modal.all"
modal.all[1,(match("tot.1v4",colnames(modal.all))):(match("pde.4v1",colnames(modal.all)))] <- est.modal[102:137,"est"] # TOT_1V4 to PDE_4V1
modal.all[1,(match("tot.1v4.se",colnames(bch.all))):(match("pde.4v1.se",colnames(bch.all)))] <- est.modal[102:137,"se"] # TOT_1V4 to PDE_4V1
modal.all[1,(match("p1",colnames(bch.all))):(match("p4",colnames(bch.all)))] <- est.modal[74:77,"est"]  # P_X1 to P_X4

#move over class probabilities
modal[,"p1"]<-modal.all[,"p1"]
modal[,"p2"]<-modal.all[,"p2"]
modal[,"p3"]<-modal.all[,"p3"]
modal[,"p4"]<-modal.all[,"p4"]
  
#move over desired class comparisons (this will differ depending on order of the classes in mplus output - determined using the area under the trajectory matrix (auc))
  #1234 (eop,ao,cl,low)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc4"])
    modal[1,1:18] <- modal.all[1,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc3"])
    modal[1,1:18] <- modal.all[1,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc2"]>auc[1,"auc4"])
    modal[1,1:18] <- modal.all[1,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc3"])
    modal[1,1:18] <- modal.all[1,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc2"])
    modal[1,1:18] <- modal.all[1,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc2"])
    modal[1,1:18] <- modal.all[1,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc4"])
    modal[1,1:18] <- modal.all[1,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc3"])
    modal[1,1:18] <- modal.all[1,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc2"]>auc[1,"auc4"])
    modal[1,1:18] <- modal.all[1,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc3"])
    modal[1,1:18] <- modal.all[1,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc2"])
    modal[1,1:18] <- modal.all[1,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc2"])
    modal[1,1:18] <- modal.all[1,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc1"]>auc[1,"auc4"])
    modal[1,1:18] <- modal.all[1,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc3"])
    modal[1,1:18] <- modal.all[1,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc1"]>auc[1,"auc4"])
    modal[1,1:18] <- modal.all[1,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc3"])
    modal[1,1:18] <- modal.all[1,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc2"])
    modal[1,1:18] <- modal.all[1,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc2"])
    modal[1,1:18] <- modal.all[1,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc1"])
    modal[1,1:18] <- modal.all[1,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc1"])
    modal[1,1:18] <- modal.all[1,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc1"])
    modal[1,1:18] <- modal.all[1,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc1"])
    modal[1,1:18] <- modal.all[1,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc1"])
    modal[1,1:18] <- modal.all[1,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc1"])
    modal[1,1:18] <- modal.all[1,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]

######################## NON-INCLUSIVE PCD ##########################
  
#set seed for non-inclusive pcds
#this will set the seed for each session in R (not each time run loops)
#to make sure you get the same results each time, run the loop only once within each R session
set.seed(80)
  
#calculate P(X=x|U) - this will give us the probability of class membership (cprobs) from unconditional latent class model
  
#for each class, multiply the individual data (responses to 6 binary indicators: cd1 to cd6) with within-class thresholds
#number of latent class indicators
#class 1
#this creates a matrix with within-class thresholds for class 1 repeated for every individual in the dataset (e.g. repeated 3039 times)
theta <- matrix(rep(coef.se.original[(match("th1.c1",rownames(coef.se.original))):(match("th6.c1",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
#if latent class indicator is present (-1)*threshold is used, if latent class indicator is absent (1)*threshold is used
P1<-exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta))
#class 2
#this creates a matrix with within-class thresholds for class 2 repeated for every individual in the dataset (e.g. repeated 3039 times)
theta <- matrix(rep(coef.se.original[(match("th1.c2",rownames(coef.se.original))):(match("th6.c2",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
P2<-exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta))
#class 3
#this creates a matrix with within-class thresholds for class 3 repeated for every individual in the dataset (e.g. repeated 3039 times)
theta <- matrix(rep(coef.se.original[(match("th1.c3",rownames(coef.se.original))):(match("th6.c3",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
P3<-exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta))
#class 4
#this creates a matrix with within-class thresholds for class 4 repeated for every individual in the dataset (e.g. repeated 3039 times)
theta <- matrix(rep(coef.se.original[(match("th1.c4",rownames(coef.se.original))):(match("th6.c4",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data)) 
P4<-exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta))

#need to replace missing data (NA) from indicators with "1" to just take the product of f(U|X) over the observed U's
#this is making MAR assumption conditional on observed U's (same as FIML)
#however, missingness must not depend on unobserved U's or latent class
P1[is.na(P1)] <- 1
P2[is.na(P2)] <- 1
P3[is.na(P3)] <- 1
P4[is.na(P4)] <- 1
  
#multiply the 6 probabilities (for each latent class indicator) in each row of P* then multiply by class probabilities that were saved earlier in "p.original"  
N1<-apply(P1,1,prod)*p.original[1,1]
N2<-apply(P2,1,prod)*p.original[2,1]
N3<-apply(P3,1,prod)*p.original[3,1]
N4<-apply(P4,1,prod)*p.original[4,1]    
  
#this gives us the probability of class membership for each person in the dataset 
P<-cbind(N1/(N1+N2+N3+N4),N2/(N1+N2+N3+N4),N3/(N1+N2+N3+N4),N4/(N1+N2+N3+N4))

colnames(P) <-  c("cprob1","cprob2","cprob3","cprob4")

#check against cprobs exported from mplus (they should match except for rounding)
cprobs <- read.table("bch.txt", sep="", na.strings="*", header=FALSE)
cprobs <- cprobs[,c(11:14)]
colnames(cprobs) <- c("cprob1","cprob2","cprob3","cprob4")

#we will now impute class membership 40 times for each person using their probability of class membership
#we will create 40 imputed datasets (chosen to keep Monte Carlo error at less than 10% of standard error for parameters in regression model for Y)
imp.n <- 40
#create a matrix to store imputed class membership for each person (sample*imp.n=3039*40)
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
  imp.subset <- cbind(data.cc,subset(imp, imp[,1]==l))
  prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))
  } 
  
#create "imp.txt" file for Mplus to call imputed datasets
imp.txt <- matrix(NA,imp.n,1)
for(l in 1:imp.n) {
  imp.txt[l,1] <- paste0("imp_", l, ".dat")
  } 
write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
#run mediation model using non-inclusive PCD
runModels("4f npcd mediation.inp")
  
#read in parameters from npcd mediation model 
est.npcd <- readModels("4f npcd mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "npcd.all"
npcd.all[1,(match("tot.1v4",colnames(npcd.all))):(match("pde.4v1",colnames(npcd.all)))] <- est.npcd[102:137,"est"] # TOT_1V4 to PDE_4V1
npcd.all[1,(match("tot.1v4.se",colnames(npcd.all))):(match("pde.4v1.se",colnames(npcd.all)))] <- est.npcd[102:137,"se"] # TOT_1V4 to PDE_4V1
npcd.all[1,(match("p1",colnames(npcd.all))):(match("p4",colnames(npcd.all)))] <- est.npcd[74:77,"est"]  # P_X1 to P_X4
  
#move over class probabilities
npcd[,"p1"]<-npcd.all[,"p1"]
npcd[,"p2"]<-npcd.all[,"p2"]
npcd[,"p3"]<-npcd.all[,"p3"]
npcd[,"p4"]<-npcd.all[,"p4"]
  
#move over desired class comparisons (this will differ depending on order of the classes in Mplus output - determined using the area under the trajectory matrix (auc))
  #1234 (eop,ao,cl,low) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc4"])
    npcd[1,1:18] <- npcd.all[1,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc3"])
    npcd[1,1:18] <- npcd.all[1,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc2"]>auc[1,"auc4"])
    npcd[1,1:18] <- npcd.all[1,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc3"])
    npcd[1,1:18] <- npcd.all[1,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc2"])
    npcd[1,1:18] <- npcd.all[1,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc2"])
    npcd[1,1:18] <- npcd.all[1,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc4"])
    npcd[1,1:18] <- npcd.all[1,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc3"])
    npcd[1,1:18] <- npcd.all[1,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc2"]>auc[1,"auc4"])
    npcd[1,1:18] <- npcd.all[1,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc3"])
    npcd[1,1:18] <- npcd.all[1,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc2"])
    npcd[1,1:18] <- npcd.all[1,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc2"])
    npcd[1,1:18] <- npcd.all[1,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc1"]>auc[1,"auc4"])
    npcd[1,1:18] <- npcd.all[1,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc3"])
    npcd[1,1:18] <- npcd.all[1,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc1"]>auc[1,"auc4"])
    npcd[1,1:18] <- npcd.all[1,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc3"])
    npcd[1,1:18] <- npcd.all[1,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc2"])
    npcd[1,1:18] <- npcd.all[1,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc2"])
    npcd[1,1:18] <- npcd.all[1,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc1"])
    npcd[1,1:18] <- npcd.all[1,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc1"])
    npcd[1,1:18] <- npcd.all[1,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc1"])
    npcd[1,1:18] <- npcd.all[1,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc1"])
    npcd[1,1:18] <- npcd.all[1,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc1"])
    npcd[1,1:18] <- npcd.all[1,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc1"])
    npcd[1,1:18] <- npcd.all[1,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]

######################## INCLUSIVE PCD ##########################

#set seed for inclusive pcd
set.seed(81)

#set working directory using file paths saved at the start
setwd(sim.data)

#run inclusive latent class model and export out cprobs
runModels("4g inc latent class.inp")

#read in auc parameters from inclusive latent class model
est.inc <- readModels("4g inc latent class.out", what="parameters")$parameters$`unstandardized`
auc.inc[1,(match("auc1",colnames(auc.inc))):(match("auc4",colnames(auc.inc)))] <- est.inc[40:43,"est"]  # AUC1 to AUC4

#read in class probabilities from "inc.txt" (as we did with the bch weights exported from the unconditional latent class model) 
cprobs <- read.table("inc.txt", sep="", na.strings="*", header=FALSE)
#we could also extract modal classes from here (v12) to use the inclusive modal approach
P <- cprobs[,c(11:14)]
colnames(P) <- c("cprob1","cprob2","cprob3","cprob4")

#we will now impute class membership 40 times for each person using their probability of class membership
#create a matrix to store imputed class membership for each person (sample*imp.n=3039*40)
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
  imp.subset <- cbind(data.cc,subset(imp, imp[,1]==l))
  prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))
} 

#create "imp.txt" file for Mplus to call imputed datasets
imp.txt <- matrix(NA,imp.n,1)
for(l in 1:imp.n) {
  imp.txt[l,1] <- paste0("imp_", l, ".dat")
} 
write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#run mediation model using non-inclusive PCD
runModels("4h incpcd mediation.inp")

#read in parameters from incpcd mediation model 
est.incpcd <- readModels("4h incpcd mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "npcd.all"
incpcd.all[1,(match("tot.1v4",colnames(incpcd.all))):(match("pde.4v1",colnames(incpcd.all)))] <- est.incpcd[102:137,"est"] # TOT_1V4 to PDE_4V1
incpcd.all[1,(match("tot.1v4.se",colnames(incpcd.all))):(match("pde.4v1.se",colnames(incpcd.all)))] <- est.incpcd[102:137,"se"] # TOT_1V4 to PDE_4V1
incpcd.all[1,(match("p1",colnames(incpcd.all))):(match("p4",colnames(incpcd.all)))] <- est.incpcd[74:77,"est"]  # P_X1 to P_X4

#move over class probabilities
incpcd[,"p1"]<-incpcd.all[,"p1"]
incpcd[,"p2"]<-incpcd.all[,"p2"]
incpcd[,"p3"]<-incpcd.all[,"p3"]
incpcd[,"p4"]<-incpcd.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in mplus output - determined using the area under the trajectory matrix (auc))
  #1234 (eop,ao,cl,low) 
  if(auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"])
    incpcd[1,1:18] <- incpcd.all[1,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"])
    incpcd[1,1:18] <- incpcd.all[1,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"])
    incpcd[1,1:18] <- incpcd.all[1,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc4"]>auc.inc[1,"auc3"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"])
    incpcd[1,1:18] <- incpcd.all[1,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc3"]>auc.inc[1,"auc4"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"])
    incpcd[1,1:18] <- incpcd.all[1,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc4"]>auc.inc[1,"auc3"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"])
    incpcd[1,1:18] <- incpcd.all[1,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"])
    incpcd[1,1:18] <- incpcd.all[1,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"])
    incpcd[1,1:18] <- incpcd.all[1,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"])
    incpcd[1,1:18] <- incpcd.all[1,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"]
     & auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"])
    incpcd[1,1:18] <- incpcd.all[1,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"])
    incpcd[1,1:18] <- incpcd.all[1,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"]
     & auc.inc[1,"auc1"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"])
    incpcd[1,1:18] <- incpcd.all[1,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"])
    incpcd[1,1:18] <- incpcd.all[1,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"])
    incpcd[1,1:18] <- incpcd.all[1,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"] & auc.inc[1,"auc1"]>auc.inc[1,"auc4"])
    incpcd[1,1:18] <- incpcd.all[1,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"]
     & auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc1"]>auc.inc[1,"auc3"])
    incpcd[1,1:18] <- incpcd.all[1,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc2"])
    incpcd[1,1:18] <- incpcd.all[1,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"]
     & auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc1"]>auc.inc[1,"auc2"])
    incpcd[1,1:18] <- incpcd.all[1,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"] & auc.inc[1,"auc4"]>auc.inc[1,"auc1"])
    incpcd[1,1:18] <- incpcd.all[1,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"] & auc.inc[1,"auc3"]>auc.inc[1,"auc1"])
    incpcd[1,1:18] <- incpcd.all[1,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc4"] & auc.inc[1,"auc4"]>auc.inc[1,"auc1"])
    incpcd[1,1:18] <- incpcd.all[1,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"]
     & auc.inc[1,"auc2"]>auc.inc[1,"auc1"] & auc.inc[1,"auc2"]>auc.inc[1,"auc3"] & auc.inc[1,"auc3"]>auc.inc[1,"auc1"])
    incpcd[1,1:18] <- incpcd.all[1,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc3"]>auc.inc[1,"auc4"]
     & auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc2"]>auc.inc[1,"auc1"])
    incpcd[1,1:18] <- incpcd.all[1,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc.inc[1,"auc4"]>auc.inc[1,"auc1"] & auc.inc[1,"auc4"]>auc.inc[1,"auc2"] & auc.inc[1,"auc4"]>auc.inc[1,"auc3"]
     & auc.inc[1,"auc3"]>auc.inc[1,"auc1"] & auc.inc[1,"auc3"]>auc.inc[1,"auc2"] & auc.inc[1,"auc2"]>auc.inc[1,"auc1"])
    incpcd[1,1:18] <- incpcd.all[1,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]

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
model.y <- logistf(int_18~drugs_18+female+rs_soc, data=all.data)
#regression model for mediator m 
model.m <- glm(drugs_18~female+rs_soc, family=binomial, data=all.data)

#perturbing the coefficients once around coefficients in models using variance-covariance matrix 
#initially use zero for exposure coefficients (and exposure-mediator interactions in outcome model)
#once latent class exposure has been imputed in subsequent runs, these will become coefficients from models
beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y)),0,0,0,0,0,0)
beta.m<-c(rmnorm(1,coef(model.m),vcov(model.m)),0,0,0)

#we will create 60 imputed datasets for class membership (chosen to keep Monte Carlo error at 10% of SE or less)
imp.n <- 60
#we will allow 20 iterations between saving out imputed class membership
cycles <- 20
#we will allow a burn in of 100 iterations before starting to save out imputed class membership
burnin <- 100

#create a matrix to save the results from every iteration (all results are recorded to assess convergence later)
#rows=iterations=(cycles*imp.n+burnin)
#columns=iteration number, coefficients from regression models, cell sizes and 2 flags for issues=(6+(latent.classes-1)*3+latent.classes*4+confounders*2)=35
results <- matrix(NA,cycles*imp.n+burnin,6+(latent.classes-1)*3+latent.classes*4+confounders*2) 

colnames(results) <- c("iteration", #iteration number
                       "b0.y","b1.y","b2.y","b3.y","b4.y","b5.y","b6.y","b7.y","b8.y","b9.y", #beta coefficients from regression model for the outcome
                       "b0.m","b1.m","b2.m","b3.m","b4.m","b5.m", #beta coefficients from regression model for the mediator
                       "x1.m0.y0","x1.m1.y0","x2.m0.y0","x2.m1.y0","x3.m0.y0","x3.m1.y0","x4.m0.y0","x4.m1.y0","x1.m0.y1","x1.m1.y1","x2.m0.y1","x2.m1.y1","x3.m0.y1","x3.m1.y1","x4.m0.y1","x4.m1.y1", #cell sizes from crosstabs for classes by mediator by outcome
                       "zero.cell", #flag for presence of zero cells in crosstabs
                       "large.th" #flag for a within-class threshold in unconditional latent class model that was out of bounds after perturbing
                        ) 

#the first column is simply an indicator of iteration number (range from 1 to 1300)
results[,1]<-1:(cycles*imp.n+burnin)

#create a matrix to store cell sizes from crosstabs for classes by mediator by outcome 
xmy<-matrix(NA,1,latent.classes*4)
colnames(xmy) <-  c("x1m0y0","x1m1y0","x2m0y0","x2m1y0","x3m0y0","x3m1y0","x4m0y0","x4m1y0","x1m0y1","x1m1y1","x2m0y1","x2m1y1","x3m0y1","x3m1y1","x4m0y1","x4m1y1")
#create a matrix to store presence of zero cells in this crosstabs
zero.cell<-matrix(NA,1,1)
#create a matrix to store presence of a within-class threshold in unconditional latent class model that was out of bounds after perturbing (e.g., not corresponding to 0 to 100% probability) 
large.th<-matrix(NA,1,1)

#create a matrix to store imputed class membership for each person (sample*imp.n=3039*60)
imp <- matrix(NA,sample*imp.n,latent.classes+1)
colnames(imp) <- c("imp","x1","x2","x3","x4") 
#first column is simply an indicator for imputation number (range from 1 to 60)
for(h in 0:(imp.n-1)) {
  imp [c(sample*h+1:sample),1]<-h+1}

#we will now create a loop to repeat steps 2 and 3 below 1300 (imp.n*cycles+burnin) times and save estimates after every 20 iterations (after a 100 iteration burn in)
for(j in 1:(imp.n*cycles+burnin)) {
  results[j,2:11]<-beta.y #save the perturbed coefficients from regression model for the outcome "beta.y"
  results[j,12:17]<-beta.m #save the perturbed coefficients from regression model for the mediator "beta.m"    
  results[j,18:33]<-xmy #save cell sizes from crosstabs for classes by mediator by outcome "xmy"
  results[j,34]<-zero.cell #save a flag for presence of zero cells in this crosstabs "zero.cell"
  results[j,35]<-large.th #save a flag for a within-class threshold in unconditional latent class model that was out of bounds after perturbing "large.th"
  
  #this saves each individual's imputed class membership 60 times (generated later in the script) into the matrix "imp" we created earlier     
  for(l in 1:imp.n) {
    if (j==(cycles*l+burnin)) imp[c(sample*(l-1)+1:sample),2:(latent.classes+1)]<-x
  } 
  
  #STEP 2a:
  #we will now perturb the parameters (within-class thresholds and class intercepts) that we saved earlier "coef.se.original" from the unconditional latent class model 
  #NB we do not need ro reorder the parameters as we did in the simulation model as when there is missing data in the indicators (as there is here) Mplus numbers the parameters differently (see tech1)
  coef<-rmnorm(1,coef.se.original[,"est"],cov)

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
  rownames(coef) <- c("th1.c1","th2.c1","th3.c1","th4.c1","th5.c1","th6.c1",
                      "th1.c2","th2.c2","th3.c2","th4.c2","th5.c2","th6.c1",
                      "th1.c3","th2.c3","th3.c3","th4.c3","th5.c3","th6.c1",
                      "th1.c4","th2.c4","th3.c4","th4.c4","th5.c4","th6.c1",
                      "int.c1","int.c2","int.c3")
  colnames(coef) <- c("est")
        
  #now calculate P(X=x|U) - without perturbing this give us the probability of class membership (cprobs) from unconditional latent class model
  
  #for each class, multiply the individual data (responses to 6 binary indicators: cd1 to cd6) with within-class thresholds
  #class 1
  #this creates a matrix with within-class thresholds for class 1 repeated for every individual in the dataset (e.g. repeated 3039 times)
  theta <- matrix(rep(coef[(match("th1.c1",rownames(coef.se.original))):(match("th6.c1",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
  #if latent class indicator is present (-1)*threshold is used, if latent class indicator is absent (1)*threshold is used
  P1<-exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta))
  #class 2
  theta <- matrix(rep(coef[(match("th1.c2",rownames(coef.se.original))):(match("th6.c2",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
  P2<-exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta))
  #class 3
  theta <- matrix(rep(coef[(match("th1.c3",rownames(coef.se.original))):(match("th6.c3",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
  P3<-exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta))
  #class 4
  theta <- matrix(rep(coef[(match("th1.c4",rownames(coef.se.original))):(match("th6.c4",rownames(coef.se.original))),"est"],each=nrow(all.data)),ncol=indicators,nrow=nrow(all.data))
  P4<-exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta)/(1+exp((-1)^all.data[,(match("cd1",colnames(all.data))):(match("cd6",colnames(all.data)))]*theta))
  
  #need to replace missing data (NA) from indicators with "1" to just take the product of f(U|X) over the observed U's
  #this is making MAR assumption conditional on observed U's (same as FIML)
  #however, missingness must not depend on unobserved U's or latent class
  P1[is.na(P1)] <- 1
  P2[is.na(P2)] <- 1
  P3[is.na(P3)] <- 1
  P4[is.na(P4)] <- 1
    
  #we need to create a new matrix of class probabilities using the perturbed class intercepts from the unconditional latent class model
  p <- matrix(NA,sample,latent.classes)
  p[,1] <- exp(coef["int.c1","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
  p[,2] <- exp(coef["int.c2","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))
  p[,3] <- exp(coef["int.c3","est"])/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))  
  p[,4] <- 1/(1 + exp(coef["int.c1","est"]) + exp(coef["int.c2","est"]) + exp(coef["int.c3","est"]))

  #in the first round this is multiplying P1 etc by class probabilities (which will be the same for each person)
  #in subsequent rounds it uses the propensity score for classes given confounders (which differs depending on a persons confounder values)
  #this propensity score is derived later in the code
  if (j==1) P1.ps <- cbind(P1,p[,1])
  else P1.ps <- cbind(P1,ps[,1])
  if (j==1) P2.ps <- cbind(P2,p[,2])
  else P2.ps <- cbind(P2,ps[,2])
  if (j==1) P3.ps <- cbind(P3,p[,3])
  else P3.ps <- cbind(P3,ps[,3])
  if (j==1) P4.ps <- cbind(P4,p[,4])
  else P4.ps <- cbind(P4,ps[,4])  
    
  #this gives us the probability of class membership for each person in the dataset
  #these will differ slightly to "cprobs" that can be exported from the unconditional latent class model due to perturbing
  N1<-apply(P1.ps,1,prod)
  N2<-apply(P2.ps,1,prod)
  N3<-apply(P3.ps,1,prod)
  N4<-apply(P4.ps,1,prod)    
    
  #STEP 2b - Bayes rule:
  
  #combine P(X=x|U) along with coefficients from regression models for outcome and mediator
  #this will calculate probabilities P(X=x|Y,M,U) for each class (x = 1,..k) for each individual
  #confounders (Z) need to be used to calculate P(Y=y|M,X,Z) and P(M=m|X,Z)
  #we will do this via a number a steps below
  
  #P(Y=1|M,X,Z): calculate probability that the outcome (Y) = 1 given mediator (M), exposure (X, latent classes) , and confounders(Z)
  #this uses the coefficients from the regression model for the outcome (intercept (b0), coef for M (b1), coefs for Z (b2-b3), coefs for X (b4-b6) and coefs for XM interaction (b7-b9))
  #class 1
  num1 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$drugs_18+results[j,"b2.y"]*all.data$female+results[j,"b3.y"]*all.data$rs_soc+results[j,"b4.y"]+results[j,"b7.y"]*all.data$drugs_18)
  Y11 <- num1/(1+num1)  
  #class 2
  num2 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$drugs_18+results[j,"b2.y"]*all.data$female+results[j,"b3.y"]*all.data$rs_soc+results[j,"b5.y"]+results[j,"b8.y"]*all.data$drugs_18)
  Y21 <- num2/(1+num2)  
  #class 3
  num3 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$drugs_18+results[j,"b2.y"]*all.data$female+results[j,"b3.y"]*all.data$rs_soc+results[j,"b6.y"]+results[j,"b9.y"]*all.data$drugs_18)
  Y31 <- num3/(1+num3)
  #class 4
  num4 <- exp(results[j,"b0.y"]+results[j,"b1.y"]*all.data$drugs_18+results[j,"b2.y"]*all.data$female+results[j,"b3.y"]*all.data$rs_soc)
  Y41 <- num4/(1+num4)

  #P(Y=0|M,X,Z): calculate probability that the outcome (Y) = 0 given mediator (M), exposure (X, latent classes), and confounders (Z) 
  Y10 <- 1/(1+num1)
  Y20 <- 1/(1+num2)
  Y30 <- 1/(1+num3)
  Y40 <- 1/(1+num4)
    
  #P(M=1|X,Z): calculate probability that the mediator (M) = 1 given exposure (X, latent classes) and confounders (Z)
  #this uses the coefficients from the regression model for the mediator (intercept (b0), coefs for Z (b1-b2) and coefs for X (b3-b5))
  #class 1 
  num1 <- exp(results[j,"b0.m"]+results[j,"b1.m"]*all.data$female+results[j,"b2.m"]*all.data$rs_soc+results[j,"b3.m"])
  M11 <- num1/(1+num1)
  #class 2
  num2 <- exp(results[j,"b0.m"]+results[j,"b1.m"]*all.data$female+results[j,"b2.m"]*all.data$rs_soc+results[j,"b4.m"])
  M21 <- num2/(1+num2) 
  #class 3
  num3 <- exp(results[j,"b0.m"]+results[j,"b1.m"]*all.data$female+results[j,"b2.m"]*all.data$rs_soc+results[j,"b5.m"])
  M31 <- num3/(1+num3)   
  #class 4
  num4 <- exp(results[j,"b0.m"]+results[j,"b1.m"]*all.data$female+results[j,"b2.m"]*all.data$rs_soc)
  M41 <- num4/(1+num4) 
  
  #P(M=0|X,Z)    
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
  #using each individual's observed data on the mediator and outcome
  Q<-Q11
  #probabilities for those with mediator and outcome absent
  Q[all.data$int_18==0 & all.data$drugs_18==0,]<-Q00[all.data$int_18==0 & all.data$drugs_18==0,]
  #probabilities for those with mediator absent and outcome present
  Q[all.data$int_18==1 & all.data$drugs_18==0,]<-Q01[all.data$int_18==1 & all.data$drugs_18==0,] 
  #probabilities for those with mediator present and outcome absent
  Q[all.data$int_18==0 & all.data$drugs_18==1,]<-Q10[all.data$int_18==0 & all.data$drugs_18==1,]
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
  all.data[,"int1"] <- all.data[,"x1"]*all.data[,"drugs_18"]
  all.data[,"int2"] <- all.data[,"x2"]*all.data[,"drugs_18"]
  all.data[,"int3"] <- all.data[,"x3"]*all.data[,"drugs_18"]
  all.data[,"int4"] <- all.data[,"x4"]*all.data[,"drugs_18"]
  
  ###################################################
  #we will now perform some checks on the cell sizes from crosstabs for classes by mediator by outcome 
  #we will add these cell sizes into "results" to make traceplots to assess convergence
  
  #create a subset of the data for Y=0 
  no.y <- subset(all.data, all.data$int_18==0)
  #create a subset of the data for Y=1
  yes.y <- subset(all.data, all.data$int_18==1)
  
  #crosstabs for x1 and m for those with y=0
  x1m.no.y<-table(no.y$x1,no.y$drugs_18)
  
  #this is to make sure matrix is 2 by 2 even when there are zero cells
  if (nrow(x1m.no.y)==1) x1m.no.y <- rbind(x1m.no.y,matrix(0,1,2))
  #crosstabs for x2 and m for those with y=0 
  x2m.no.y<-table(no.y$x2,no.y$drugs_18)
  if (nrow(x2m.no.y)==1) x2m.no.y <- rbind(x2m.no.y,matrix(0,1,2))
  #crosstabs for x3 and m for those with y=0 
  x3m.no.y<-table(no.y$x3,no.y$drugs_18)
  if (nrow(x3m.no.y)==1) x3m.no.y <- rbind(x3m.no.y,matrix(0,1,2))
  #crosstabs for x4 and m for those with y=0 
  x4m.no.y<-table(no.y$x4,no.y$drugs_18)
  if (nrow(x4m.no.y)==1) x4m.no.y <- rbind(x4m.no.y,matrix(0,1,2))
  #crosstabs for x1 and m for those with y=1 
  x1m.yes.y<-table(yes.y$x1,yes.y$drugs_18)
  if (nrow(x1m.yes.y)==1) x1m.yes.y <- rbind(x1m.yes.y,matrix(0,1,2))    
  x2m.yes.y<-table(yes.y$x2,yes.y$drugs_18)
  if (nrow(x2m.yes.y)==1) x2m.yes.y <- rbind(x2m.yes.y,matrix(0,1,2))    
  x3m.yes.y<-table(yes.y$x3,yes.y$drugs_18)
  if (nrow(x3m.yes.y)==1) x3m.yes.y <- rbind(x3m.yes.y,matrix(0,1,2))
  x4m.yes.y<-table(yes.y$x4,yes.y$drugs_18)
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
    
  #fit the firth logistic regression model for P(Y|X,M,Z) to obtain updated parameter estimates now classes have been imputed
  if (zero.cell<2) model.y <- logistf(int_18~drugs_18+female+rs_soc+x1+x2+x3+int1+int2+int3, data=all.data)
  #if there is more than 1 zero cell, even firth logistic regression does not converge, therefore it is necessary to remove XM interactions from regression model
  if (zero.cell>1) model.y <- logistf(int_18~drugs_18+female+rs_soc+x1+x2+x3, data=all.data)
  
  #fit the logistic regression model for P(M|X,Z) to obtain updated parameter estimates
  model.m <- glm(drugs_18~female+rs_soc+x1+x2+x3, family=binomial, data=all.data)
  
  #perturbing the coefficients once around coefficients in models using variance-covariance matrix 
  if (zero.cell>1) beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y)),0,0,0)
  if (zero.cell<2) beta.y<-c(rmnorm(1,coef(model.y),vcov(model.y))) 
  beta.m<-c(rmnorm(1,coef(model.m),vcov(model.m)))
  
  ##Propensity scores
  #################################################
  #because including confounders, class probabilities will need to be replaced by propensity scores
  #these will be created by regressing the generated class variable on the confounders
  #first combine x1-x4 dummy variables into nominal x
  
  for(k in 1:sample) {
  if (all.data[k,"x1"]==1) all.data[k,22]<-1
  if (all.data[k,"x2"]==1) all.data[k,22]<-2
  if (all.data[k,"x3"]==1) all.data[k,22]<-3
  if (all.data[k,"x4"]==1) all.data[k,22]<-4}
  
  #next, use "detect_separation" with a generalised linear model to test for the occurrence of complete or quasi-complete separation
  #this will detect if any of the parameters have infinite maximum likelihood estimates - this is quite common with multinomial logistic regression
    
  #class order is low, CL, EOP, AO therefore in our binomial models we want x1 to represent the reference group
  #this code will return a "0" for any parameter estimates that are finite and return "inf" for any that are infinite 
  data.x2 <- subset(all.data, x3==0 & x4==0)
  model.x2 <- glm(x2~female+rs_soc, family=binomial, data=data.x2, method = "detect_separation")
  x2 <- coef(model.x2)
  data.x3 <- subset(all.data, x2==0 & x4==0)
  model.x3 <- glm(x3~female+rs_soc, family=binomial, data=data.x3, method = "detect_separation")
  x3 <- coef(model.x3)
  data.x4 <- subset(all.data, x2==0 & x3==0)
  model.x4 <- glm(x4~female+rs_soc, family=binomial, data=data.x4, method = "detect_separation")
  x4 <- coef(model.x4)
    
  check <- matrix(NA,1,9)
  check[1,1:3] <- x2
  check[1,4:6] <- x3
  check[1,7:9] <- x4
    
  #here we will flag the occurrence of any infinite parameter estimates in the model
  inf <- matrix(NA,1,1)
  inf <- 0
  for(k in 1:9) {
    if (check[k]!=0)  inf <- 1}
    
  #if there is an infinite parameter estimate we can use "brmultinom" to address this
  #if not, a standard multinomial regression can be used which is much quicker to run
  #the adjusted score approach to bias reduction that brmultinom implements (type = "AS_mean") is an alternative to 
  #maximum likelihood that results in estimates with smaller asymptotic bias that are also *always* finite, 
  #even in cases of complete or quasi-complete separation.
  if (inf==1) model.x <- brmultinom(x ~ female+rs_soc, data = all.data, type = "AS_mean")
  if (inf==0) model.x <- multinom(x ~ female+rs_soc, data = all.data, Hess=T)
    
  coef.model.x <- coef(model.x)
  coef.x <- matrix(NA,1,9)
  coef.x[1,1]<-coef.model.x[1,1]
  coef.x[1,2]<-coef.model.x[1,2]
  coef.x[1,3]<-coef.model.x[1,3]
  coef.x[1,4]<-coef.model.x[2,1]
  coef.x[1,5]<-coef.model.x[2,2]
  coef.x[1,6]<-coef.model.x[2,3]
  coef.x[1,7]<-coef.model.x[3,1]
  coef.x[1,8]<-coef.model.x[3,2]
  coef.x[1,9]<-coef.model.x[3,3]
    
  #perturb the estimates using their (co)variance matrix (code differs slightly depending which type of multinomial model is used)
  if (inf==1) beta.x<-c(rmnorm(1,coef.x,vcov(model.x)))
  if (inf==0) beta.x<-c(rmnorm(1,coef.x,solve(model.x$Hessian)))
    
  #create the propensity score
  num1 <- exp(beta.x[1]+beta.x[2]*all.data$female+beta.x[3]*all.data$rs_soc)
  num2 <- exp(beta.x[4]+beta.x[5]*all.data$female+beta.x[6]*all.data$rs_soc)
  num3 <- exp(beta.x[7]+beta.x[8]*all.data$female+beta.x[9]*all.data$rs_soc)
  denom <- 1+ num1 + num2 + num3
  ps1 <- 1/denom
  ps2 <- num1/denom
  ps3 <- num2/denom
  ps4 <- num3/denom
  ps <- matrix(NA,sample,4)
  ps[,1] <- ps1
  ps[,2] <- ps2
  ps[,3] <- ps3  
  ps[,4] <- ps4
  
} #end of iteration loop

################################################
#assessing convergence

#create plot of parameters and cell sizes across each iteration
results.df <- as.data.frame(results)

#remove the first iteration when no data for X (latent classes) 
results.df <- subset(results.df, iteration>1)

#to check autocorrelation for each beta after burn in of 100 iterations
b4.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b4.y"]
b5.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b5.y"]
b6.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b6.y"]
b7.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b7.y"]
b8.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b8.y"]
b9.y <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b9.y"]
b3.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b3.m"]
b4.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b4.m"]
b5.m <- results.df[(burnin-1):((cycles*imp.n+burnin)-1),"b5.m"]

#can use a separate folder for updated PCD because of large number of files created (e.g., 60 imputed datasets)
setwd(upcd.files)

#create a pdf showing traceplots, histograms and autocorrelation plots for each cefficient, and traceplots and histograms for cell sizes across all iterations

  pdf(file="trace_plot.pdf") 
  plot(results.df$iteration, results.df$b4.y,
       xlab = "iteration",
       ylab = "b4.y (y on class 1 vs class 4)")
  lines(results.df$iteration, results.df$b4.y)
  acf(b4.y, lag.max=50)
  plot(results.df$iteration, results.df$b5.y,
       xlab = "iteration",
       ylab = "b5.y (y on class 2 vs class 4)")
  lines(results.df$iteration, results.df$b5.y)
  acf(b5.y, lag.max=50)
  plot(results.df$iteration, results.df$b6.y,
       xlab = "iteration",
       ylab = "b6.y (y on class 3 vs class 4)")
  lines(results.df$iteration, results.df$b6.y)
  acf(b6.y, lag.max=50)
  plot(results.df$iteration, results.df$b7.y,
       xlab = "iteration",
       ylab = "b7.y (y on class 1 x m interaction)")
  lines(results.df$iteration, results.df$b7.y)
  acf(b7.y, lag.max=50)
  plot(results.df$iteration, results.df$b8.y,
       xlab = "iteration",
       ylab = "b8.y (y on class 2 x m interaction)")
  lines(results.df$iteration, results.df$b8.y)
  acf(b8.y, lag.max=50)
  plot(results.df$iteration, results.df$b9.y,
       xlab = "iteration",
       ylab = "b9.y (y on class 3 x m interaction)")
  lines(results.df$iteration, results.df$b9.y)
  acf(b9.y, lag.max=50)
  plot(results.df$iteration, results.df$b3.m,
       xlab = "iteration",
       ylab = "b3.m (m on class 1 vs class 4)")
  lines(results.df$iteration, results.df$b3.m)
  acf(b3.m, lag.max=50)
  plot(results.df$iteration, results.df$b4.m,
       xlab = "iteration",
       ylab = "b4.m (m on class 2 vs class 4)")
  lines(results.df$iteration, results.df$b4.m)
  acf(b4.m, lag.max=50)
  plot(results.df$iteration, results.df$b5.m,
       xlab = "iteration",
       ylab = "b5.m on class 3 vs class 4")
  lines(results.df$iteration, results.df$b5.m)
  acf(b5.m, lag.max=50)
  plot(results.df$iteration, results.df$x1.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=1 m=0 y=0")
  lines(results.df$iteration, results.df$x1.m0.y0)
  plot(results.df$iteration, results.df$x1.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=1 m=1 y=0")
  lines(results.df$iteration, results.df$x1.m1.y0)
  plot(results.df$iteration, results.df$x2.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=2 m=0 y=0")
  lines(results.df$iteration, results.df$x2.m0.y0)
  plot(results.df$iteration, results.df$x2.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=2 m=1 y=0")
  lines(results.df$iteration, results.df$x2.m1.y0)
  plot(results.df$iteration, results.df$x3.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=3 m=0 y=0")
  lines(results.df$iteration, results.df$x3.m0.y0)
  plot(results.df$iteration, results.df$x3.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=3 m=1 y=0")
  lines(results.df$iteration, results.df$x3.m1.y0)
  plot(results.df$iteration, results.df$x4.m0.y0,
       xlab = "iteration",
       ylab = "cell size x=4 m=0 y=0")
  lines(results.df$iteration, results.df$x4.m0.y0)
  plot(results.df$iteration, results.df$x4.m1.y0,
       xlab = "iteration",
       ylab = "cell size x=4 m=1 y=0")
  lines(results.df$iteration, results.df$x4.m1.y0)
  plot(results.df$iteration, results.df$x1.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=1 m=0 y=1")
  lines(results.df$iteration, results.df$x1.m0.y1)
  plot(results.df$iteration, results.df$x1.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=1 m=1 y=1")
  lines(results.df$iteration, results.df$x1.m1.y1)
  plot(results.df$iteration, results.df$x2.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=2 m=0 y=1")
  lines(results.df$iteration, results.df$x2.m0.y1)
  plot(results.df$iteration, results.df$x2.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=2 m=1 y=1")
  lines(results.df$iteration, results.df$x2.m1.y1)
  plot(results.df$iteration, results.df$x3.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=3 m=0 y=1")
  lines(results.df$iteration, results.df$x3.m0.y1)
  plot(results.df$iteration, results.df$x3.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=3 m=1 y=1")
  lines(results.df$iteration, results.df$x3.m1.y1)
  plot(results.df$iteration, results.df$x4.m0.y1,
       xlab = "iteration",
       ylab = "cell size x=4 m=0 y=1")
  lines(results.df$iteration, results.df$x4.m0.y1)
  plot(results.df$iteration, results.df$x4.m1.y1,
       xlab = "iteration",
       ylab = "cell size x=4 m=1 y=1")
  lines(results.df$iteration, results.df$x4.m1.y1)
  dev.off()

#analysis
  
#combine imputed class membership stored in "imp.n" with original data
#prepare 1 mplus .dat file for each imputed dataset to run mediation model (60 .dat files should be created)
for(l in 1:imp.n) {
  imp.subset <- cbind(data.cc,subset(imp, imp[,1]==l))
  prepareMplusData(imp.subset, file=paste0("imp_", l, ".dat"))}

#create "imp.txt" file for mplus to call imputed datasets
imp.txt <- matrix(NA,imp.n,1)
for(l in 1:imp.n) {
  imp.txt[l,1] <- paste0("imp_", l, ".dat")}
write.table(imp.txt, file="imp.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
#run mediation model using updated PCD
runModels("4i pcd mediation.inp")  
    
#read in parameters from updated PCD mediation model 
est.pcd <- readModels("4i pcd mediation.out", what="parameters")$parameters$`unstandardized`
#store only the parameters (mediation effects and class probabilities) and their SEs in the matrix created at the start "upcd.all"
pcd.all[1,(match("tot.1v4",colnames(pcd.all))):(match("pde.4v1",colnames(pcd.all)))] <- est.pcd[102:137,"est"] # TOT_1V4 to PDE_4V1
pcd.all[1,(match("tot.1v4.se",colnames(pcd.all))):(match("pde.4v1.se",colnames(pcd.all)))] <- est.pcd[102:137,"se"] # TOT_1V4 to PDE_4V1
pcd.all[1,(match("p1",colnames(pcd.all))):(match("p4",colnames(pcd.all)))] <- est.pcd[74:77,"est"]  # P_X1 to P_X4

#flags for potential issues: "fixed.th","large.th","largest.th","largest.th.se","zero.cell","zero.cell.imp"

#flag if fixed threshold in unconditional model which results in missing value in tech3 covariance matrix of parameters
pcd.all[1,"fixed.th"] <- 0 #set flag to 0 to start
for(k in 1:(indicators*latent.classes+(latent.classes-1))*(indicators*latent.classes+(latent.classes-1))) { #size of tech3 matrix is 27x27 (729) 
  if (tech3[k]==999) pcd.all[i,"fixed.th"] <- 1}

#when perturbing within-class thresholds, those with a large standard error can go out of bounds (e.g., corresponding to a probability that is not between 0 and 100%)
#we created a flag for this in "results" so here we will move this into "pcd.all" and record whether this was the case in any iteration
pcd.all[1,"large.th"] <- 0 #set flag to 0 to start
results[1,"large.th"]<-0 
for(k in 1:(cycles*imp.n+burnin)) {
  if (results[k,"large.th"]==1) pcd.all[i,"large.th"] <- 1}

#we also want a flag so that we know the largest standard error for a threshold in the unconditional model (and the threshold this SE corresponds to)
#this flag captures the largest SE
pcd.all[1,"largest.th.se"] <- max(coef.se.original[,"se"])
#this flag records the threshold that corresponds to largest SE
pcd.all[1,"largest.th"] <- coef.se.original[which.max(coef.se.original[,"se"]),"est"]

#we created a flag for number of zero cells in the crosstabs for classes by mediator by outcome in "results"
#we will move this flag to "pcd.all" and record the largest number of zero cells in any iteration
results[1,"zero.cell"]<-0 #change row 1 in "results" to 0 as latent classes not imputed in first iteration
pcd.all[1,"zero.cell"] <- max(results[,"zero.cell"])

#we will also flag if there was a zero cell (in the crosstabs for classes by mediator by outcome) in any of the iterations that we saved as one of the 60 imputed datasets
zero.cell.imp <- matrix(NA,1,imp.n)
for(l in 1:imp.n) {
  zero.cell.imp[1,l] <- results[cycles*l+burnin,"zero.cell"]}

pcd.all[1,"zero.cell.imp"] <- 0 #set flag to 0 to start
for(k in 1:imp.n) {
  if (zero.cell.imp[1,k]>0) pcd.all[i,"zero.cell.imp"] <- 1}

#############################################################

#move over class probabilities
pcd[,"p1"]<-pcd.all[,"p1"]
pcd[,"p2"]<-pcd.all[,"p2"]
pcd[,"p3"]<-pcd.all[,"p3"]
pcd[,"p4"]<-pcd.all[,"p4"]

#move over desired class comparisons (this will differ depending on order of the classes in mplus output - determined using the area under the trajectory matrix (auc))
  #1234 (eop,ao,cl,low) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc4"])
    pcd[1,1:18] <- pcd.all[1,c(1,4,7,13,16,19,25,28,31,37,40,43,49,52,55,61,64,67)]
  #1243 (eop,ao,low,cl)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc3"])
    pcd[1,1:18] <- pcd.all[1,c(2,5,10,14,17,22,26,29,34,38,41,46,50,53,58,62,65,70)]
  #1324 (eop,cl,ao,low) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc2"]>auc[1,"auc4"])
    pcd[1,1:18] <- pcd.all[1,c(1,7,4,13,19,16,25,31,28,37,43,40,49,55,52,61,67,64)]
  #1342 (eop,cl,low,ao)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc3"])
    pcd[1,1:18] <- pcd.all[1,c(2,10,5,14,22,17,26,34,29,38,46,41,50,58,53,62,70,65)]
  #1423 (eop,low,ao,cl) 
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc2"])
    pcd[1,1:18] <- pcd.all[1,c(3,8,11,15,20,23,27,32,35,39,44,47,51,56,59,63,68,71)] 
  #1432 (eop,low,cl,ao)
  if(auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc2"])
    pcd[1,1:18] <- pcd.all[1,c(3,11,8,15,23,20,27,35,32,39,47,44,51,59,56,63,71,68)] 
  #2134 (ao,eop,cl,low)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc3"]>auc[1,"auc4"])
    pcd[1,1:18] <- pcd.all[1,c(4,1,7,16,13,19,28,25,31,40,37,43,52,49,55,64,61,67)]
  #2143 (ao,eop,low,cl)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc3"])
    pcd[1,1:18] <- pcd.all[1,c(5,2,10,17,14,22,29,26,34,41,38,46,53,50,58,65,62,70)]
  #2314 (ao,cl,eop,low)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc2"]>auc[1,"auc4"])
    pcd[1,1:18] <- pcd.all[1,c(7,1,4,19,13,16,31,25,28,43,37,40,55,49,52,67,61,64)]
  #2341 (ao,cl,low,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc3"])
    pcd[1,1:18] <- pcd.all[1,c(10,2,5,22,14,17,34,26,29,46,38,41,58,50,53,70,62,65)]
  #2413 (ao,low,eop,cl)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc2"])
    pcd[1,1:18] <- pcd.all[1,c(8,3,11,20,15,23,32,27,35,44,39,47,56,51,59,68,63,71)]
  #2431 (ao,low,cl,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc1"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc2"])
    pcd[1,1:18] <- pcd.all[1,c(11,3,8,23,15,20,35,27,32,47,39,44,59,51,56,71,63,68)]
  #3124 (cl,eop,ao,low)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc1"]>auc[1,"auc4"])
    pcd[1,1:18] <- pcd.all[1,c(4,7,1,16,19,13,28,31,25,40,43,37,52,55,49,64,67,61)]
  #3142 (cl,eop,low,ao)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc3"])
    pcd[1,1:18] <- pcd.all[1,c(5,10,2,17,22,14,29,34,26,41,46,38,53,58,50,65,70,62)]
  #3214 (cl,ao,eop,low)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc1"]>auc[1,"auc4"])
    pcd[1,1:18] <- pcd.all[1,c(7,4,1,19,16,13,31,28,25,43,40,37,55,52,49,67,64,61)]
  #3241 (cl,ao,low,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc1"]>auc[1,"auc3"])
    pcd[1,1:18] <- pcd.all[1,c(10,5,2,22,17,14,34,29,26,46,41,38,58,53,50,70,65,62)]
  #3412 (cl,low,eop,ao)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc2"])
    pcd[1,1:18] <- pcd.all[1,c(8,11,3,20,23,15,32,35,27,44,47,39,56,59,51,68,71,63)]
  #3421 (cl,low,ao,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc1"]>auc[1,"auc2"])
    pcd[1,1:18] <- pcd.all[1,c(11,8,3,23,20,15,35,32,27,47,44,39,59,56,51,71,68,63)]
  #4123 (low,eop,ao,cl)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc1"])
    pcd[1,1:18] <- pcd.all[1,c(6,9,12,18,21,24,30,33,36,42,45,48,54,57,60,66,69,72)]
  #4132 (low,eop,cl,ao)
  if(auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc2"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc1"])
    pcd[1,1:18] <- pcd.all[1,c(6,12,9,18,24,21,30,36,33,42,48,45,54,60,57,66,72,69)]
  #4213 (low,ao,eop,cl)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc4"] & auc[1,"auc4"]>auc[1,"auc1"])
    pcd[1,1:18] <- pcd.all[1,c(9,6,12,21,18,24,33,30,36,45,42,48,57,54,60,69,66,72)]
  #4231 (low,ao,cl,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc2"]>auc[1,"auc1"] & auc[1,"auc2"]>auc[1,"auc3"] & auc[1,"auc3"]>auc[1,"auc1"])
    pcd[1,1:18] <- pcd.all[1,c(12,6,9,24,18,21,36,30,33,48,42,45,60,54,57,72,66,69)]
  #4312 (low,cl,eop,ao)
  if(auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc3"]>auc[1,"auc4"]
     & auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc1"])
    pcd[1,1:18] <- pcd.all[1,c(9,12,6,21,24,18,33,36,30,45,48,42,57,60,54,69,72,66)]
  #4321 (low,cl,ao,eop)
  if(auc[1,"auc4"]>auc[1,"auc1"] & auc[1,"auc4"]>auc[1,"auc2"] & auc[1,"auc4"]>auc[1,"auc3"]
     & auc[1,"auc3"]>auc[1,"auc1"] & auc[1,"auc3"]>auc[1,"auc2"] & auc[1,"auc2"]>auc[1,"auc1"])
    pcd[1,1:18] <- pcd.all[1,c(12,9,6,24,21,18,36,33,30,48,45,42,60,57,54,72,69,66)]

##################################################

#now combine estimates from all methods in one table to export
estimates <- rbind(onestep,bch,modal,npcd,incpcd,pcd)
rownames(estimates) <- c("onestep","bch","modal","npcd","incpcd","upcd")

#put results in a table and export to a results folder
setwd(results.files)
write.table(estimates, file="table.estimates.txt", sep = ",")
write.table(onestep.all, file="table.onestep.txt", sep = ",")
write.table(bch.all, file="table.bch.txt", sep = ",")
write.table(modal.all, file="table.modal.txt", sep = ",")
write.table(npcd.all, file="table.npcd.txt", sep = ",")
write.table(incpcd.all, file="table.incpcd.txt", sep = ",")
write.table(pcd.all, file="table.upcd.txt", sep = ",")
 

#create a figure to show results
##################################################

dt.total.eop <- data.frame(method = c("onestep","bch","modal","npcd","incpcd","upcd"),
                       estimand = c("total","total","total","total","total","total"),
                       comparison = c("eop","eop","eop","eop","eop","eop"),
                       logrr = estimates[,"tot.eop"],
                       lci = estimates[,"tot.eop"]-1.96*estimates[,"tot.eop.se"],
                       uci = estimates[,"tot.eop"]+1.96*estimates[,"tot.eop.se"])

dt.total.ao <- data.frame(method = c("onestep","bch","modal","npcd","incpcd","upcd"),
                           estimand = c("total","total","total","total","total","total"),
                           comparison = c("ao","ao","ao","ao","ao","ao"),
                          logrr = estimates[,"tot.ao"],
                          lci = estimates[,"tot.ao"]-1.96*estimates[,"tot.ao.se"],
                          uci = estimates[,"tot.ao"]+1.96*estimates[,"tot.ao.se"])

dt.total.cl <- data.frame(method = c("onestep","bch","modal","npcd","incpcd","upcd"),
                          estimand = c("total","total","total","total","total","total"),
                          comparison = c("cl","cl","cl","cl","cl","cl"),                       
                          logrr = estimates[,"tot.cl"],
                          lci = estimates[,"tot.cl"]-1.96*estimates[,"tot.cl.se"],
                          uci = estimates[,"tot.cl"]+1.96*estimates[,"tot.cl.se"])

dt.tie.eop <- data.frame(method = c("onestep","bch","modal","npcd","incpcd","upcd"),
                     estimand = c("tie","tie","tie","tie","tie","tie"),
                     comparison = c("eop","eop","eop","eop","eop","eop"),
                     logrr = estimates[,"tie.eop"],
                     lci = estimates[,"tie.eop"]-1.96*estimates[,"tie.eop.se"],
                     uci = estimates[,"tie.eop"]+1.96*estimates[,"tie.eop.se"])
                     
dt.tie.ao <- data.frame(method = c("onestep","bch","modal","npcd","incpcd","upcd"),
                     estimand = c("tie","tie","tie","tie","tie","tie"),
                     comparison = c("ao","ao","ao","ao","ao","ao"),                     
                     logrr = estimates[,"tie.ao"],
                     lci = estimates[,"tie.ao"]-1.96*estimates[,"tie.ao.se"],
                     uci = estimates[,"tie.ao"]+1.96*estimates[,"tie.ao.se"])

dt.tie.cl <- data.frame(method = c("onestep","bch","modal","npcd","incpcd","upcd"),
                     estimand = c("tie","tie","tie","tie","tie","tie"),
                     comparison = c("cl","cl","cl","cl","cl","cl"),  
                     logrr = estimates[,"tie.cl"],
                     lci = estimates[,"tie.cl"]-1.96*estimates[,"tie.cl.se"],
                     uci = estimates[,"tie.cl"]+1.96*estimates[,"tie.cl.se"])

dt.pde.eop <- data.frame(method = c("onestep","bch","modal","npcd","incpcd","upcd"),
                     estimand = c("pde","pde","pde","pde","pde","pde"),
                     comparison = c("eop","eop","eop","eop","eop","eop"),
                     logrr = estimates[,"pde.eop"],
                     lci = estimates[,"pde.eop"]-1.96*estimates[,"pde.eop.se"],
                     uci = estimates[,"pde.eop"]+1.96*estimates[,"pde.eop.se"])


dt.pde.ao <- data.frame(method = c("onestep","bch","modal","npcd","incpcd","upcd"),
                         estimand = c("pde","pde","pde","pde","pde","pde"),
                         comparison = c("ao","ao","ao","ao","ao","ao"),
                         logrr = estimates[,"pde.ao"],
                         lci = estimates[,"pde.ao"]-1.96*estimates[,"pde.ao.se"],
                         uci = estimates[,"pde.ao"]+1.96*estimates[,"pde.ao.se"])

dt.pde.cl <- data.frame(method = c("onestep","bch","modal","npcd","incpcd","upcd"),
                        estimand = c("pde","pde","pde","pde","pde","pde"),
                        comparison = c("cl","cl","cl","cl","cl","cl"),
                        logrr = estimates[,"pde.cl"],
                        lci = estimates[,"pde.cl"]-1.96*estimates[,"pde.cl.se"],
                        uci = estimates[,"pde.cl"]+1.96*estimates[,"pde.cl.se"])

dt.total <- rbind(dt.total.eop,dt.total.ao,dt.total.cl)
dt.tie <- rbind(dt.tie.eop,dt.tie.ao,dt.tie.cl)
dt.pde <- rbind(dt.pde.eop,dt.pde.ao,dt.pde.cl)

#uses package forcats
comparison.labs <- c(eop = "eop", ao = "ao", cl = "cl")

#total
pl1 <- ggplot(data = dt.total) +
  geom_point(aes(x=fct_inorder(method), y=logrr), color= "black") +
  geom_errorbar(aes(x=method, ymin=lci, ymax= uci), width = 0.4, color ="black", size = 1) +
  geom_text(aes(x=method, y=lci, label = round(lci,2)), size= 3, vjust = 1) +
  geom_text(aes(x=method, y=uci, label = round(uci,2)), size= 3, vjust = -1) +
  theme_bw() +
  labs(y = "log risk ratio", x = "method") +
  theme(text = element_text(size = 12))  +
  facet_wrap(~fct_inorder(comparison), 
             labeller = labeller(comparison = comparison.labs),
             nrow=1) +
  ylim(-0.5, 1.75)
pl1

#indirect
pl2 <- ggplot(data = dt.tie) +
  geom_point(aes(x=fct_inorder(method), y=logrr), color= "black") +
  geom_errorbar(aes(x=method, ymin=lci, ymax= uci), width = 0.4, color ="black", size = 1) +
  geom_text(aes(x=method, y=lci, label = round(lci,2)), size= 3, vjust = 1) +
  geom_text(aes(x=method, y=uci, label = round(uci,2)), size= 3, vjust = -1) +
  theme_bw() +
  labs(y = "log risk ratio", x = "method") +
  theme(text = element_text(size = 12))  +
  facet_wrap(~fct_inorder(comparison), 
             labeller = labeller(comparison = comparison.labs),
             nrow=1) +
  ylim(-0.3, 0.4) 
pl2

#direct
pl3 <- ggplot(data = dt.pde) +
  geom_point(aes(x=fct_inorder(method), y=logrr), color= "black") +
  geom_errorbar(aes(x=method, ymin=lci, ymax= uci), width = 0.4, color ="black", size = 1) +
  geom_text(aes(x=method, y=lci, label = round(lci,2)), size= 3, vjust = 1) +
  geom_text(aes(x=method, y=uci, label = round(uci,2)), size= 3, vjust = -1) +
  theme_bw() +
  labs(y = "log risk ratio", x = "method") +
  theme(text = element_text(size = 12))  +
  facet_wrap(~fct_inorder(comparison), 
             labeller = labeller(comparison = comparison.labs),
             nrow=1) +
  ylim(-0.5, 1.75) 
pl3

  ################################################end of script#####################################################################