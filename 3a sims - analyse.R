#R script for analysing simulation results
#script uses vignettes for rsimsum

#this script needs to be run after "2a sims.R"
#############################################################

############### INSTALL AND LOAD PACKAGES ###################

#install the required packages (if not already installed)
list.of.packages <- c("dplyr", "scales", "ggplot2", "tidyr", "rsimsum", "knitr", "forcats")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] 
if(length(new.packages)) install.packages(new.packages)

#call all required packages
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages)

#check versions of required packages
getNamespaceVersion("dplyr") #1.0.7
getNamespaceVersion("scales") #1.2.0     
getNamespaceVersion("ggplot2") #3.3.5
getNamespaceVersion("tidyr") #1.1.3
getNamespaceVersion("rsimsum") #0.11.2
getNamespaceVersion("knitr") #1.33
getNamespaceVersion("forcats")

##################################
#true mediation effects
betas <- c(0.356674944,0.251314428,0.133531393,0.05702914,0.080874497,0.011555227,0.299645804,0.170439931,0.121976166)

sims <- 500
##################################

#first bring estimates data back in
#we want to create a table in long format with all 500 sims, 3 dgms and 5 methods
#we want one estimand (with SE) per table (3 mediation effects for 3 comparisons = 9 tables)

setwd("file location/all sims/results/good entropy")

onestep.good <- read.table(file="table.onestep.txt", sep = ",", header=TRUE)
bch.good <- read.table(file="table.bch.txt", sep = ",", header=TRUE)
modal.good <- read.table(file="table.modal.txt", sep = ",", header=TRUE)
npcd.good <- read.table(file="table.npcd.txt", sep = ",", header=TRUE)
incpcd.good <- read.table(file="table.incpcd.txt", sep = ",", header=TRUE)
pcd.good <- read.table(file="table.pcd.txt", sep = ",", header=TRUE)

#create exclusion flag based on the largest standard error for the within-class thresholds in unconditional model
#this causes non-convergence (seen in traceplots for cell sizes)
#traceplots checked for sims from highest standard error until convergence is seen
#this seems to correspond to SEs less than 2*mean of all SEs
pcd.good$exclude <- 0
se.max <- 2*mean(pcd.good$largest.th.se)
for(k in 1:sims) {
  if (pcd.good[k,"largest.th.se"]>se.max) pcd.good[k,"exclude"]<-1}

table(pcd.good$exclude)
#we will use the same rule for SEs for within-class thresholds in onestep model
onestep.good$exclude <- 0
se.max <- 2*mean(onestep.good$largest.th.se)
for(k in 1:sims) {
  if (onestep.good[k,"largest.th.se"]>se.max) onestep.good[k,"exclude"]<-1}

#43 to exclude in total (16 for pcd, 23 for onestep and 4 for both)
table(pcd.good$exclude,onestep.good$exclude)

######################################################

setwd("file location/all sims/results/medium entropy")

onestep.med <- read.table(file="table.onestep.txt", sep = ",", header=TRUE)
bch.med <- read.table(file="table.bch.txt", sep = ",", header=TRUE)
modal.med <- read.table(file="table.modal.txt", sep = ",", header=TRUE)
npcd.med <- read.table(file="table.npcd.txt", sep = ",", header=TRUE)
incpcd.med <- read.table(file="table.incpcd.txt", sep = ",", header=TRUE)
pcd.med <- read.table(file="table.pcd.txt", sep = ",", header=TRUE)

#create exclusion flag based on the largest standard error for the within-class thresholds in unconditional model
#this causes non-convergence (seen in traceplots for cell sizes)
#traceplots checked for sims from highest standard error until convergence is seen
#this seems to correspond to SEs less than 2*mean of all SEs
pcd.med$exclude <- 0
se.max <- 2*mean(pcd.med$largest.th.se)
for(k in 1:sims) {
  if (pcd.med[k,"largest.th.se"]>se.max) pcd.med[k,"exclude"]<-1}

table(pcd.med$exclude)
#we will use the same rule for SEs for within-class thresholds in onestep model
onestep.med$exclude <- 0
se.max <- 2*mean(onestep.med$largest.th.se)
for(k in 1:sims) {
  if (onestep.med[k,"largest.th.se"]>se.max) onestep.med[k,"exclude"]<-1}

#66 to exclude in total (32 for pcd, 26 for onestep and 8 for both)
table(pcd.med$exclude,onestep.med$exclude)

######################################################

setwd("file location/all sims/results/poor entropy")

onestep.poor <- read.table(file="table.onestep.txt", sep = ",", header=TRUE)
bch.poor <- read.table(file="table.bch.txt", sep = ",", header=TRUE)
modal.poor <- read.table(file="table.modal.txt", sep = ",", header=TRUE)
npcd.poor <- read.table(file="table.npcd.txt", sep = ",", header=TRUE)
incpcd.poor <- read.table(file="table.incpcd.txt", sep = ",", header=TRUE)
pcd.poor <- read.table(file="table.pcd.txt", sep = ",", header=TRUE)

#create exclusion flag based on the largest standard error for the within-class thresholds in unconditional model
#this sometimes causes non-convergence (seen in traceplots for cell sizes)
#traceplots checked for sims from highest standard error until convergence is seen
#with poor entropy, sims with $SEs more than 2*mean of SEs across all sims do not always show non-convergence on the trace plot 
pcd.poor$exclude <- 0
se.max <- 2*mean(pcd.poor$largest.th.se)
for(k in 1:sims) {
  if (pcd.poor[k,"largest.th.se"]>se.max) pcd.poor[k,"exclude"]<-1}

table(pcd.poor$exclude)
#we will use the same rule for SEs for within-class thresholds in onestep model
onestep.poor$exclude <- 0
se.max <- 2*mean(onestep.poor$largest.th.se)
for(k in 1:sims) {
  if (onestep.poor[k,"largest.th.se"]>se.max) onestep.poor[k,"exclude"]<-1}

#76 to exclude in total (35 for pcd, 33 for onestep and 8 for both)
table(pcd.poor$exclude,onestep.poor$exclude)

######################################################

#create dataframes for each estimand
dgm <- 3
method <- 6

#create a matrix with exclusion criteria (from within-class thresholds in unconditional model or onestep)
exclude <- matrix(NA,sims,dgm)
colnames(exclude) <- c("good","medium","poor")
exclude[,"good"] <- pcd.good[,"exclude"]
for(k in 1:sims) {
if (onestep.good[k,"exclude"]==1) exclude[k,"good"]=1 }
exclude[,"medium"] <- pcd.med[,"exclude"]
for(k in 1:sims) {
  if (onestep.med[k,"exclude"]==1) exclude[k,"medium"]=1 }
exclude[,"poor"] <- pcd.poor[,"exclude"]
for(k in 1:sims) {
  if (onestep.poor[k,"exclude"]==1) exclude[k,"poor"]=1 }

te.eop <- matrix(NA,sims*dgm*method,6)
colnames(te.eop) <- c("rep","dgm","method","se","theta","exclude")
te.eop <- data.frame(te.eop)
te.ao <- matrix(NA,sims*dgm*method,6)
colnames(te.ao) <- c("rep","dgm","method","se","theta","exclude")
te.ao <- data.frame(te.ao)
te.cl <- matrix(NA,sims*dgm*method,6)
colnames(te.cl) <- c("rep","dgm","method","se","theta","exclude")
te.cl <- data.frame(te.cl)

tnie.eop <- matrix(NA,sims*dgm*method,6)
colnames(tnie.eop) <- c("rep","dgm","method","se","theta","exclude")
tnie.eop <- data.frame(tnie.eop)
tnie.ao <- matrix(NA,sims*dgm*method,6)
colnames(tnie.ao) <- c("rep","dgm","method","se","theta","exclude")
tnie.ao <- data.frame(tnie.ao)
tnie.cl <- matrix(NA,sims*dgm*method,6)
colnames(tnie.cl) <- c("rep","dgm","method","se","theta","exclude")
tnie.cl <- data.frame(tnie.cl)

pnde.eop <- matrix(NA,sims*dgm*method,6)
colnames(pnde.eop) <- c("rep","dgm","method","se","theta","exclude")
pnde.eop <- data.frame(pnde.eop)
pnde.ao <- matrix(NA,sims*dgm*method,6)
colnames(pnde.ao) <- c("rep","dgm","method","se","theta","exclude")
pnde.ao <- data.frame(pnde.ao)
pnde.cl <- matrix(NA,sims*dgm*method,6)
colnames(pnde.cl) <- c("rep","dgm","method","se","theta","exclude")
pnde.cl <- data.frame(pnde.cl)

#total effect for eop vs low
te.eop[1:sims,"rep"] <- 1:sims
te.eop[1:sims,"dgm"] <- "good"
te.eop[1:sims,"method"] <- "onestep"
te.eop[1:sims,"se"] <- onestep.good[1:sims,"tot.eop.se"]
te.eop[1:sims,"theta"] <- onestep.good[1:sims,"tot.eop"]
te.eop[1:sims,"exclude"] <- exclude[1:sims,"good"]
te.eop[(sims+1):(sims*2),"rep"] <- 1:sims
te.eop[(sims+1):(sims*2),"dgm"] <- "good"
te.eop[(sims+1):(sims*2),"method"] <- "bch"
te.eop[(sims+1):(sims*2),"se"] <- bch.good[1:sims,"tot.eop.se"]
te.eop[(sims+1):(sims*2),"theta"] <- bch.good[1:sims,"tot.eop"]
te.eop[(sims+1):(sims*2),"exclude"] <- exclude[1:sims,"good"]
te.eop[(sims*2+1):(sims*3),"rep"] <- 1:sims
te.eop[(sims*2+1):(sims*3),"dgm"] <- "good"
te.eop[(sims*2+1):(sims*3),"method"] <- "modal"
te.eop[(sims*2+1):(sims*3),"se"] <- modal.good[1:sims,"tot.eop.se"]
te.eop[(sims*2+1):(sims*3),"theta"] <- modal.good[1:sims,"tot.eop"]
te.eop[(sims*2+1):(sims*3),"exclude"] <- exclude[1:sims,"good"]
te.eop[(sims*3+1):(sims*4),"rep"] <- 1:sims
te.eop[(sims*3+1):(sims*4),"dgm"] <- "good"
te.eop[(sims*3+1):(sims*4),"method"] <- "npcd"
te.eop[(sims*3+1):(sims*4),"se"] <- npcd.good[1:sims,"tot.eop.se"]
te.eop[(sims*3+1):(sims*4),"theta"] <- npcd.good[1:sims,"tot.eop"]
te.eop[(sims*3+1):(sims*4),"exclude"] <- exclude[1:sims,"good"]
te.eop[(sims*4+1):(sims*5),"rep"] <- 1:sims
te.eop[(sims*4+1):(sims*5),"dgm"] <- "good"
te.eop[(sims*4+1):(sims*5),"method"] <- "incpcd"
te.eop[(sims*4+1):(sims*5),"se"] <- incpcd.good[1:sims,"tot.eop.se"]
te.eop[(sims*4+1):(sims*5),"theta"] <- incpcd.good[1:sims,"tot.eop"]
te.eop[(sims*4+1):(sims*5),"exclude"] <- exclude[1:sims,"good"]
te.eop[(sims*5+1):(sims*6),"rep"] <- 1:sims
te.eop[(sims*5+1):(sims*6),"dgm"] <- "good"
te.eop[(sims*5+1):(sims*6),"method"] <- "upcd"
te.eop[(sims*5+1):(sims*6),"se"] <- pcd.good[1:sims,"tot.eop.se"]
te.eop[(sims*5+1):(sims*6),"theta"] <- pcd.good[1:sims,"tot.eop"]
te.eop[(sims*5+1):(sims*6),"exclude"] <- exclude[1:sims,"good"]
te.eop[(sims*6+1):(sims*7),"rep"] <- 1:sims
te.eop[(sims*6+1):(sims*7),"dgm"] <- "medium"
te.eop[(sims*6+1):(sims*7),"method"] <- "onestep"
te.eop[(sims*6+1):(sims*7),"se"] <- onestep.med[1:sims,"tot.eop.se"]
te.eop[(sims*6+1):(sims*7),"theta"] <- onestep.med[1:sims,"tot.eop"]
te.eop[(sims*6+1):(sims*7),"exclude"] <- exclude[1:sims,"medium"]
te.eop[(sims*7+1):(sims*8),"rep"] <- 1:sims
te.eop[(sims*7+1):(sims*8),"dgm"] <- "medium"
te.eop[(sims*7+1):(sims*8),"method"] <- "bch"
te.eop[(sims*7+1):(sims*8),"se"] <- bch.med[1:sims,"tot.eop.se"]
te.eop[(sims*7+1):(sims*8),"theta"] <- bch.med[1:sims,"tot.eop"]
te.eop[(sims*7+1):(sims*8),"exclude"] <- exclude[1:sims,"medium"]
te.eop[(sims*8+1):(sims*9),"rep"] <- 1:sims
te.eop[(sims*8+1):(sims*9),"dgm"] <- "medium"
te.eop[(sims*8+1):(sims*9),"method"] <- "modal"
te.eop[(sims*8+1):(sims*9),"se"] <- modal.med[1:sims,"tot.eop.se"]
te.eop[(sims*8+1):(sims*9),"theta"] <- modal.med[1:sims,"tot.eop"]
te.eop[(sims*8+1):(sims*9),"exclude"] <- exclude[1:sims,"medium"]
te.eop[(sims*9+1):(sims*10),"rep"] <- 1:sims
te.eop[(sims*9+1):(sims*10),"dgm"] <- "medium"
te.eop[(sims*9+1):(sims*10),"method"] <- "npcd"
te.eop[(sims*9+1):(sims*10),"se"] <- npcd.med[1:sims,"tot.eop.se"]
te.eop[(sims*9+1):(sims*10),"theta"] <- npcd.med[1:sims,"tot.eop"]
te.eop[(sims*9+1):(sims*10),"exclude"] <- exclude[1:sims,"medium"]
te.eop[(sims*10+1):(sims*11),"rep"] <- 1:sims
te.eop[(sims*10+1):(sims*11),"dgm"] <- "medium"
te.eop[(sims*10+1):(sims*11),"method"] <- "incpcd"
te.eop[(sims*10+1):(sims*11),"se"] <- incpcd.med[1:sims,"tot.eop.se"]
te.eop[(sims*10+1):(sims*11),"theta"] <- incpcd.med[1:sims,"tot.eop"]
te.eop[(sims*10+1):(sims*11),"exclude"] <- exclude[1:sims,"medium"]
te.eop[(sims*11+1):(sims*12),"rep"] <- 1:sims
te.eop[(sims*11+1):(sims*12),"dgm"] <- "medium"
te.eop[(sims*11+1):(sims*12),"method"] <- "upcd"
te.eop[(sims*11+1):(sims*12),"se"] <- pcd.med[1:sims,"tot.eop.se"]
te.eop[(sims*11+1):(sims*12),"theta"] <- pcd.med[1:sims,"tot.eop"]
te.eop[(sims*11+1):(sims*12),"exclude"] <- exclude[1:sims,"medium"]
te.eop[(sims*12+1):(sims*13),"rep"] <- 1:sims
te.eop[(sims*12+1):(sims*13),"dgm"] <- "poor"
te.eop[(sims*12+1):(sims*13),"method"] <- "onestep"
te.eop[(sims*12+1):(sims*13),"se"] <- onestep.poor[1:sims,"tot.eop.se"]
te.eop[(sims*12+1):(sims*13),"theta"] <- onestep.poor[1:sims,"tot.eop"]
te.eop[(sims*12+1):(sims*13),"exclude"] <- exclude[1:sims,"poor"]
te.eop[(sims*13+1):(sims*14),"rep"] <- 1:sims
te.eop[(sims*13+1):(sims*14),"dgm"] <- "poor"
te.eop[(sims*13+1):(sims*14),"method"] <- "bch"
te.eop[(sims*13+1):(sims*14),"se"] <- bch.poor[1:sims,"tot.eop.se"]
te.eop[(sims*13+1):(sims*14),"theta"] <- bch.poor[1:sims,"tot.eop"]
te.eop[(sims*13+1):(sims*14),"exclude"] <- exclude[1:sims,"poor"]
te.eop[(sims*14+1):(sims*15),"rep"] <- 1:sims
te.eop[(sims*14+1):(sims*15),"dgm"] <- "poor"
te.eop[(sims*14+1):(sims*15),"method"] <- "modal"
te.eop[(sims*14+1):(sims*15),"se"] <- modal.poor[1:sims,"tot.eop.se"]
te.eop[(sims*14+1):(sims*15),"theta"] <- modal.poor[1:sims,"tot.eop"]
te.eop[(sims*14+1):(sims*15),"exclude"] <- exclude[1:sims,"poor"]
te.eop[(sims*15+1):(sims*16),"rep"] <- 1:sims
te.eop[(sims*15+1):(sims*16),"dgm"] <- "poor"
te.eop[(sims*15+1):(sims*16),"method"] <- "npcd"
te.eop[(sims*15+1):(sims*16),"se"] <- npcd.poor[1:sims,"tot.eop.se"]
te.eop[(sims*15+1):(sims*16),"theta"] <- npcd.poor[1:sims,"tot.eop"]
te.eop[(sims*15+1):(sims*16),"exclude"] <- exclude[1:sims,"poor"]
te.eop[(sims*16+1):(sims*17),"rep"] <- 1:sims
te.eop[(sims*16+1):(sims*17),"dgm"] <- "poor"
te.eop[(sims*16+1):(sims*17),"method"] <- "incpcd"
te.eop[(sims*16+1):(sims*17),"se"] <- incpcd.poor[1:sims,"tot.eop.se"]
te.eop[(sims*16+1):(sims*17),"theta"] <- incpcd.poor[1:sims,"tot.eop"]
te.eop[(sims*16+1):(sims*17),"exclude"] <- exclude[1:sims,"poor"]
te.eop[(sims*17+1):(sims*18),"rep"] <- 1:sims
te.eop[(sims*17+1):(sims*18),"dgm"] <- "poor"
te.eop[(sims*17+1):(sims*18),"method"] <- "upcd"
te.eop[(sims*17+1):(sims*18),"se"] <- pcd.poor[1:sims,"tot.eop.se"]
te.eop[(sims*17+1):(sims*18),"theta"] <- pcd.poor[1:sims,"tot.eop"]
te.eop[(sims*17+1):(sims*18),"exclude"] <- exclude[1:sims,"poor"]

#total effect for ao vs low
te.ao[1:sims,"rep"] <- 1:sims
te.ao[1:sims,"dgm"] <- "good"
te.ao[1:sims,"method"] <- "onestep"
te.ao[1:sims,"se"] <- onestep.good[1:sims,"tot.ao.se"]
te.ao[1:sims,"theta"] <- onestep.good[1:sims,"tot.ao"]
te.ao[1:sims,"exclude"] <- exclude[1:sims,"good"]
te.ao[(sims+1):(sims*2),"rep"] <- 1:sims
te.ao[(sims+1):(sims*2),"dgm"] <- "good"
te.ao[(sims+1):(sims*2),"method"] <- "bch"
te.ao[(sims+1):(sims*2),"se"] <- bch.good[1:sims,"tot.ao.se"]
te.ao[(sims+1):(sims*2),"theta"] <- bch.good[1:sims,"tot.ao"]
te.ao[(sims+1):(sims*2),"exclude"] <- exclude[1:sims,"good"]
te.ao[(sims*2+1):(sims*3),"rep"] <- 1:sims
te.ao[(sims*2+1):(sims*3),"dgm"] <- "good"
te.ao[(sims*2+1):(sims*3),"method"] <- "modal"
te.ao[(sims*2+1):(sims*3),"se"] <- modal.good[1:sims,"tot.ao.se"]
te.ao[(sims*2+1):(sims*3),"theta"] <- modal.good[1:sims,"tot.ao"]
te.ao[(sims*2+1):(sims*3),"exclude"] <- exclude[1:sims,"good"]
te.ao[(sims*3+1):(sims*4),"rep"] <- 1:sims
te.ao[(sims*3+1):(sims*4),"dgm"] <- "good"
te.ao[(sims*3+1):(sims*4),"method"] <- "npcd"
te.ao[(sims*3+1):(sims*4),"se"] <- npcd.good[1:sims,"tot.ao.se"]
te.ao[(sims*3+1):(sims*4),"theta"] <- npcd.good[1:sims,"tot.ao"]
te.ao[(sims*3+1):(sims*4),"exclude"] <- exclude[1:sims,"good"]
te.ao[(sims*4+1):(sims*5),"rep"] <- 1:sims
te.ao[(sims*4+1):(sims*5),"dgm"] <- "good"
te.ao[(sims*4+1):(sims*5),"method"] <- "incpcd"
te.ao[(sims*4+1):(sims*5),"se"] <- incpcd.good[1:sims,"tot.ao.se"]
te.ao[(sims*4+1):(sims*5),"theta"] <- incpcd.good[1:sims,"tot.ao"]
te.ao[(sims*4+1):(sims*5),"exclude"] <- exclude[1:sims,"good"]
te.ao[(sims*5+1):(sims*6),"rep"] <- 1:sims
te.ao[(sims*5+1):(sims*6),"dgm"] <- "good"
te.ao[(sims*5+1):(sims*6),"method"] <- "upcd"
te.ao[(sims*5+1):(sims*6),"se"] <- pcd.good[1:sims,"tot.ao.se"]
te.ao[(sims*5+1):(sims*6),"theta"] <- pcd.good[1:sims,"tot.ao"]
te.ao[(sims*5+1):(sims*6),"exclude"] <- exclude[1:sims,"good"]
te.ao[(sims*6+1):(sims*7),"rep"] <- 1:sims
te.ao[(sims*6+1):(sims*7),"dgm"] <- "medium"
te.ao[(sims*6+1):(sims*7),"method"] <- "onestep"
te.ao[(sims*6+1):(sims*7),"se"] <- onestep.med[1:sims,"tot.ao.se"]
te.ao[(sims*6+1):(sims*7),"theta"] <- onestep.med[1:sims,"tot.ao"]
te.ao[(sims*6+1):(sims*7),"exclude"] <- exclude[1:sims,"medium"]
te.ao[(sims*7+1):(sims*8),"rep"] <- 1:sims
te.ao[(sims*7+1):(sims*8),"dgm"] <- "medium"
te.ao[(sims*7+1):(sims*8),"method"] <- "bch"
te.ao[(sims*7+1):(sims*8),"se"] <- bch.med[1:sims,"tot.ao.se"]
te.ao[(sims*7+1):(sims*8),"theta"] <- bch.med[1:sims,"tot.ao"]
te.ao[(sims*7+1):(sims*8),"exclude"] <- exclude[1:sims,"medium"]
te.ao[(sims*8+1):(sims*9),"rep"] <- 1:sims
te.ao[(sims*8+1):(sims*9),"dgm"] <- "medium"
te.ao[(sims*8+1):(sims*9),"method"] <- "modal"
te.ao[(sims*8+1):(sims*9),"se"] <- modal.med[1:sims,"tot.ao.se"]
te.ao[(sims*8+1):(sims*9),"theta"] <- modal.med[1:sims,"tot.ao"]
te.ao[(sims*8+1):(sims*9),"exclude"] <- exclude[1:sims,"medium"]
te.ao[(sims*9+1):(sims*10),"rep"] <- 1:sims
te.ao[(sims*9+1):(sims*10),"dgm"] <- "medium"
te.ao[(sims*9+1):(sims*10),"method"] <- "npcd"
te.ao[(sims*9+1):(sims*10),"se"] <- npcd.med[1:sims,"tot.ao.se"]
te.ao[(sims*9+1):(sims*10),"theta"] <- npcd.med[1:sims,"tot.ao"]
te.ao[(sims*9+1):(sims*10),"exclude"] <- exclude[1:sims,"medium"]
te.ao[(sims*10+1):(sims*11),"rep"] <- 1:sims
te.ao[(sims*10+1):(sims*11),"dgm"] <- "medium"
te.ao[(sims*10+1):(sims*11),"method"] <- "incpcd"
te.ao[(sims*10+1):(sims*11),"se"] <- incpcd.med[1:sims,"tot.ao.se"]
te.ao[(sims*10+1):(sims*11),"theta"] <- incpcd.med[1:sims,"tot.ao"]
te.ao[(sims*10+1):(sims*11),"exclude"] <- exclude[1:sims,"medium"]
te.ao[(sims*11+1):(sims*12),"rep"] <- 1:sims
te.ao[(sims*11+1):(sims*12),"dgm"] <- "medium"
te.ao[(sims*11+1):(sims*12),"method"] <- "upcd"
te.ao[(sims*11+1):(sims*12),"se"] <- pcd.med[1:sims,"tot.ao.se"]
te.ao[(sims*11+1):(sims*12),"theta"] <- pcd.med[1:sims,"tot.ao"]
te.ao[(sims*11+1):(sims*12),"exclude"] <- exclude[1:sims,"medium"]
te.ao[(sims*12+1):(sims*13),"rep"] <- 1:sims
te.ao[(sims*12+1):(sims*13),"dgm"] <- "poor"
te.ao[(sims*12+1):(sims*13),"method"] <- "onestep"
te.ao[(sims*12+1):(sims*13),"se"] <- onestep.poor[1:sims,"tot.ao.se"]
te.ao[(sims*12+1):(sims*13),"theta"] <- onestep.poor[1:sims,"tot.ao"]
te.ao[(sims*12+1):(sims*13),"exclude"] <- exclude[1:sims,"poor"]
te.ao[(sims*13+1):(sims*14),"rep"] <- 1:sims
te.ao[(sims*13+1):(sims*14),"dgm"] <- "poor"
te.ao[(sims*13+1):(sims*14),"method"] <- "bch"
te.ao[(sims*13+1):(sims*14),"se"] <- bch.poor[1:sims,"tot.ao.se"]
te.ao[(sims*13+1):(sims*14),"theta"] <- bch.poor[1:sims,"tot.ao"]
te.ao[(sims*13+1):(sims*14),"exclude"] <- exclude[1:sims,"poor"]
te.ao[(sims*14+1):(sims*15),"rep"] <- 1:sims
te.ao[(sims*14+1):(sims*15),"dgm"] <- "poor"
te.ao[(sims*14+1):(sims*15),"method"] <- "modal"
te.ao[(sims*14+1):(sims*15),"se"] <- modal.poor[1:sims,"tot.ao.se"]
te.ao[(sims*14+1):(sims*15),"theta"] <- modal.poor[1:sims,"tot.ao"]
te.ao[(sims*14+1):(sims*15),"exclude"] <- exclude[1:sims,"poor"]
te.ao[(sims*15+1):(sims*16),"rep"] <- 1:sims
te.ao[(sims*15+1):(sims*16),"dgm"] <- "poor"
te.ao[(sims*15+1):(sims*16),"method"] <- "npcd"
te.ao[(sims*15+1):(sims*16),"se"] <- npcd.poor[1:sims,"tot.ao.se"]
te.ao[(sims*15+1):(sims*16),"theta"] <- npcd.poor[1:sims,"tot.ao"]
te.ao[(sims*15+1):(sims*16),"exclude"] <- exclude[1:sims,"poor"]
te.ao[(sims*16+1):(sims*17),"rep"] <- 1:sims
te.ao[(sims*16+1):(sims*17),"dgm"] <- "poor"
te.ao[(sims*16+1):(sims*17),"method"] <- "incpcd"
te.ao[(sims*16+1):(sims*17),"se"] <- incpcd.poor[1:sims,"tot.ao.se"]
te.ao[(sims*16+1):(sims*17),"theta"] <- incpcd.poor[1:sims,"tot.ao"]
te.ao[(sims*16+1):(sims*17),"exclude"] <- exclude[1:sims,"poor"]
te.ao[(sims*17+1):(sims*18),"rep"] <- 1:sims
te.ao[(sims*17+1):(sims*18),"dgm"] <- "poor"
te.ao[(sims*17+1):(sims*18),"method"] <- "upcd"
te.ao[(sims*17+1):(sims*18),"se"] <- pcd.poor[1:sims,"tot.ao.se"]
te.ao[(sims*17+1):(sims*18),"theta"] <- pcd.poor[1:sims,"tot.ao"]
te.ao[(sims*17+1):(sims*18),"exclude"] <- exclude[1:sims,"poor"]

#total effect for cl vs low
te.cl[1:sims,"rep"] <- 1:sims
te.cl[1:sims,"dgm"] <- "good"
te.cl[1:sims,"method"] <- "onestep"
te.cl[1:sims,"se"] <- onestep.good[1:sims,"tot.cl.se"]
te.cl[1:sims,"theta"] <- onestep.good[1:sims,"tot.cl"]
te.cl[1:sims,"exclude"] <- exclude[1:sims,"good"]
te.cl[(sims+1):(sims*2),"rep"] <- 1:sims
te.cl[(sims+1):(sims*2),"dgm"] <- "good"
te.cl[(sims+1):(sims*2),"method"] <- "bch"
te.cl[(sims+1):(sims*2),"se"] <- bch.good[1:sims,"tot.cl.se"]
te.cl[(sims+1):(sims*2),"theta"] <- bch.good[1:sims,"tot.cl"]
te.cl[(sims+1):(sims*2),"exclude"] <- exclude[1:sims,"good"]
te.cl[(sims*2+1):(sims*3),"rep"] <- 1:sims
te.cl[(sims*2+1):(sims*3),"dgm"] <- "good"
te.cl[(sims*2+1):(sims*3),"method"] <- "modal"
te.cl[(sims*2+1):(sims*3),"se"] <- modal.good[1:sims,"tot.cl.se"]
te.cl[(sims*2+1):(sims*3),"theta"] <- modal.good[1:sims,"tot.cl"]
te.cl[(sims*2+1):(sims*3),"exclude"] <- exclude[1:sims,"good"]
te.cl[(sims*3+1):(sims*4),"rep"] <- 1:sims
te.cl[(sims*3+1):(sims*4),"dgm"] <- "good"
te.cl[(sims*3+1):(sims*4),"method"] <- "npcd"
te.cl[(sims*3+1):(sims*4),"se"] <- npcd.good[1:sims,"tot.cl.se"]
te.cl[(sims*3+1):(sims*4),"theta"] <- npcd.good[1:sims,"tot.cl"]
te.cl[(sims*3+1):(sims*4),"exclude"] <- exclude[1:sims,"good"]
te.cl[(sims*4+1):(sims*5),"rep"] <- 1:sims
te.cl[(sims*4+1):(sims*5),"dgm"] <- "good"
te.cl[(sims*4+1):(sims*5),"method"] <- "incpcd"
te.cl[(sims*4+1):(sims*5),"se"] <- incpcd.good[1:sims,"tot.cl.se"]
te.cl[(sims*4+1):(sims*5),"theta"] <- incpcd.good[1:sims,"tot.cl"]
te.cl[(sims*4+1):(sims*5),"exclude"] <- exclude[1:sims,"good"]
te.cl[(sims*5+1):(sims*6),"rep"] <- 1:sims
te.cl[(sims*5+1):(sims*6),"dgm"] <- "good"
te.cl[(sims*5+1):(sims*6),"method"] <- "upcd"
te.cl[(sims*5+1):(sims*6),"se"] <- pcd.good[1:sims,"tot.cl.se"]
te.cl[(sims*5+1):(sims*6),"theta"] <- pcd.good[1:sims,"tot.cl"]
te.cl[(sims*5+1):(sims*6),"exclude"] <- exclude[1:sims,"good"]
te.cl[(sims*6+1):(sims*7),"rep"] <- 1:sims
te.cl[(sims*6+1):(sims*7),"dgm"] <- "medium"
te.cl[(sims*6+1):(sims*7),"method"] <- "onestep"
te.cl[(sims*6+1):(sims*7),"se"] <- onestep.med[1:sims,"tot.cl.se"]
te.cl[(sims*6+1):(sims*7),"theta"] <- onestep.med[1:sims,"tot.cl"]
te.cl[(sims*6+1):(sims*7),"exclude"] <- exclude[1:sims,"medium"]
te.cl[(sims*7+1):(sims*8),"rep"] <- 1:sims
te.cl[(sims*7+1):(sims*8),"dgm"] <- "medium"
te.cl[(sims*7+1):(sims*8),"method"] <- "bch"
te.cl[(sims*7+1):(sims*8),"se"] <- bch.med[1:sims,"tot.cl.se"]
te.cl[(sims*7+1):(sims*8),"theta"] <- bch.med[1:sims,"tot.cl"]
te.cl[(sims*7+1):(sims*8),"exclude"] <- exclude[1:sims,"medium"]
te.cl[(sims*8+1):(sims*9),"rep"] <- 1:sims
te.cl[(sims*8+1):(sims*9),"dgm"] <- "medium"
te.cl[(sims*8+1):(sims*9),"method"] <- "modal"
te.cl[(sims*8+1):(sims*9),"se"] <- modal.med[1:sims,"tot.cl.se"]
te.cl[(sims*8+1):(sims*9),"theta"] <- modal.med[1:sims,"tot.cl"]
te.cl[(sims*8+1):(sims*9),"exclude"] <- exclude[1:sims,"medium"]
te.cl[(sims*9+1):(sims*10),"rep"] <- 1:sims
te.cl[(sims*9+1):(sims*10),"dgm"] <- "medium"
te.cl[(sims*9+1):(sims*10),"method"] <- "npcd"
te.cl[(sims*9+1):(sims*10),"se"] <- npcd.med[1:sims,"tot.cl.se"]
te.cl[(sims*9+1):(sims*10),"theta"] <- npcd.med[1:sims,"tot.cl"]
te.cl[(sims*9+1):(sims*10),"exclude"] <- exclude[1:sims,"medium"]
te.cl[(sims*10+1):(sims*11),"rep"] <- 1:sims
te.cl[(sims*10+1):(sims*11),"dgm"] <- "medium"
te.cl[(sims*10+1):(sims*11),"method"] <- "incpcd"
te.cl[(sims*10+1):(sims*11),"se"] <- incpcd.med[1:sims,"tot.cl.se"]
te.cl[(sims*10+1):(sims*11),"theta"] <- incpcd.med[1:sims,"tot.cl"]
te.cl[(sims*10+1):(sims*11),"exclude"] <- exclude[1:sims,"medium"]
te.cl[(sims*11+1):(sims*12),"rep"] <- 1:sims
te.cl[(sims*11+1):(sims*12),"dgm"] <- "medium"
te.cl[(sims*11+1):(sims*12),"method"] <- "upcd"
te.cl[(sims*11+1):(sims*12),"se"] <- pcd.med[1:sims,"tot.cl.se"]
te.cl[(sims*11+1):(sims*12),"theta"] <- pcd.med[1:sims,"tot.cl"]
te.cl[(sims*11+1):(sims*12),"exclude"] <- exclude[1:sims,"medium"]
te.cl[(sims*12+1):(sims*13),"rep"] <- 1:sims
te.cl[(sims*12+1):(sims*13),"dgm"] <- "poor"
te.cl[(sims*12+1):(sims*13),"method"] <- "onestep"
te.cl[(sims*12+1):(sims*13),"se"] <- onestep.poor[1:sims,"tot.cl.se"]
te.cl[(sims*12+1):(sims*13),"theta"] <- onestep.poor[1:sims,"tot.cl"]
te.cl[(sims*12+1):(sims*13),"exclude"] <- exclude[1:sims,"poor"]
te.cl[(sims*13+1):(sims*14),"rep"] <- 1:sims
te.cl[(sims*13+1):(sims*14),"dgm"] <- "poor"
te.cl[(sims*13+1):(sims*14),"method"] <- "bch"
te.cl[(sims*13+1):(sims*14),"se"] <- bch.poor[1:sims,"tot.cl.se"]
te.cl[(sims*13+1):(sims*14),"theta"] <- bch.poor[1:sims,"tot.cl"]
te.cl[(sims*13+1):(sims*14),"exclude"] <- exclude[1:sims,"poor"]
te.cl[(sims*14+1):(sims*15),"rep"] <- 1:sims
te.cl[(sims*14+1):(sims*15),"dgm"] <- "poor"
te.cl[(sims*14+1):(sims*15),"method"] <- "modal"
te.cl[(sims*14+1):(sims*15),"se"] <- modal.poor[1:sims,"tot.cl.se"]
te.cl[(sims*14+1):(sims*15),"theta"] <- modal.poor[1:sims,"tot.cl"]
te.cl[(sims*14+1):(sims*15),"exclude"] <- exclude[1:sims,"poor"]
te.cl[(sims*15+1):(sims*16),"rep"] <- 1:sims
te.cl[(sims*15+1):(sims*16),"dgm"] <- "poor"
te.cl[(sims*15+1):(sims*16),"method"] <- "npcd"
te.cl[(sims*15+1):(sims*16),"se"] <- npcd.poor[1:sims,"tot.cl.se"]
te.cl[(sims*15+1):(sims*16),"theta"] <- npcd.poor[1:sims,"tot.cl"]
te.cl[(sims*15+1):(sims*16),"exclude"] <- exclude[1:sims,"poor"]
te.cl[(sims*16+1):(sims*17),"rep"] <- 1:sims
te.cl[(sims*16+1):(sims*17),"dgm"] <- "poor"
te.cl[(sims*16+1):(sims*17),"method"] <- "incpcd"
te.cl[(sims*16+1):(sims*17),"se"] <- incpcd.poor[1:sims,"tot.cl.se"]
te.cl[(sims*16+1):(sims*17),"theta"] <- incpcd.poor[1:sims,"tot.cl"]
te.cl[(sims*16+1):(sims*17),"exclude"] <- exclude[1:sims,"poor"]
te.cl[(sims*17+1):(sims*18),"rep"] <- 1:sims
te.cl[(sims*17+1):(sims*18),"dgm"] <- "poor"
te.cl[(sims*17+1):(sims*18),"method"] <- "upcd"
te.cl[(sims*17+1):(sims*18),"se"] <- pcd.poor[1:sims,"tot.cl.se"]
te.cl[(sims*17+1):(sims*18),"theta"] <- pcd.poor[1:sims,"tot.cl"]
te.cl[(sims*17+1):(sims*18),"exclude"] <- exclude[1:sims,"poor"]

#tnie effect for eop vs low
tnie.eop[1:sims,"rep"] <- 1:sims
tnie.eop[1:sims,"dgm"] <- "good"
tnie.eop[1:sims,"method"] <- "onestep"
tnie.eop[1:sims,"se"] <- onestep.good[1:sims,"tie.eop.se"]
tnie.eop[1:sims,"theta"] <- onestep.good[1:sims,"tie.eop"]
tnie.eop[1:sims,"exclude"] <- exclude[1:sims,"good"]
tnie.eop[(sims+1):(sims*2),"rep"] <- 1:sims
tnie.eop[(sims+1):(sims*2),"dgm"] <- "good"
tnie.eop[(sims+1):(sims*2),"method"] <- "bch"
tnie.eop[(sims+1):(sims*2),"se"] <- bch.good[1:sims,"tie.eop.se"]
tnie.eop[(sims+1):(sims*2),"theta"] <- bch.good[1:sims,"tie.eop"]
tnie.eop[(sims+1):(sims*2),"exclude"] <- exclude[1:sims,"good"]
tnie.eop[(sims*2+1):(sims*3),"rep"] <- 1:sims
tnie.eop[(sims*2+1):(sims*3),"dgm"] <- "good"
tnie.eop[(sims*2+1):(sims*3),"method"] <- "modal"
tnie.eop[(sims*2+1):(sims*3),"se"] <- modal.good[1:sims,"tie.eop.se"]
tnie.eop[(sims*2+1):(sims*3),"theta"] <- modal.good[1:sims,"tie.eop"]
tnie.eop[(sims*2+1):(sims*3),"exclude"] <- exclude[1:sims,"good"]
tnie.eop[(sims*3+1):(sims*4),"rep"] <- 1:sims
tnie.eop[(sims*3+1):(sims*4),"dgm"] <- "good"
tnie.eop[(sims*3+1):(sims*4),"method"] <- "npcd"
tnie.eop[(sims*3+1):(sims*4),"se"] <- npcd.good[1:sims,"tie.eop.se"]
tnie.eop[(sims*3+1):(sims*4),"theta"] <- npcd.good[1:sims,"tie.eop"]
tnie.eop[(sims*3+1):(sims*4),"exclude"] <- exclude[1:sims,"good"]
tnie.eop[(sims*4+1):(sims*5),"rep"] <- 1:sims
tnie.eop[(sims*4+1):(sims*5),"dgm"] <- "good"
tnie.eop[(sims*4+1):(sims*5),"method"] <- "incpcd"
tnie.eop[(sims*4+1):(sims*5),"se"] <- incpcd.good[1:sims,"tie.eop.se"]
tnie.eop[(sims*4+1):(sims*5),"theta"] <- incpcd.good[1:sims,"tie.eop"]
tnie.eop[(sims*4+1):(sims*5),"exclude"] <- exclude[1:sims,"good"]
tnie.eop[(sims*5+1):(sims*6),"rep"] <- 1:sims
tnie.eop[(sims*5+1):(sims*6),"dgm"] <- "good"
tnie.eop[(sims*5+1):(sims*6),"method"] <- "upcd"
tnie.eop[(sims*5+1):(sims*6),"se"] <- pcd.good[1:sims,"tie.eop.se"]
tnie.eop[(sims*5+1):(sims*6),"theta"] <- pcd.good[1:sims,"tie.eop"]
tnie.eop[(sims*5+1):(sims*6),"exclude"] <- exclude[1:sims,"good"]
tnie.eop[(sims*6+1):(sims*7),"rep"] <- 1:sims
tnie.eop[(sims*6+1):(sims*7),"dgm"] <- "medium"
tnie.eop[(sims*6+1):(sims*7),"method"] <- "onestep"
tnie.eop[(sims*6+1):(sims*7),"se"] <- onestep.med[1:sims,"tie.eop.se"]
tnie.eop[(sims*6+1):(sims*7),"theta"] <- onestep.med[1:sims,"tie.eop"]
tnie.eop[(sims*6+1):(sims*7),"exclude"] <- exclude[1:sims,"medium"]
tnie.eop[(sims*7+1):(sims*8),"rep"] <- 1:sims
tnie.eop[(sims*7+1):(sims*8),"dgm"] <- "medium"
tnie.eop[(sims*7+1):(sims*8),"method"] <- "bch"
tnie.eop[(sims*7+1):(sims*8),"se"] <- bch.med[1:sims,"tie.eop.se"]
tnie.eop[(sims*7+1):(sims*8),"theta"] <- bch.med[1:sims,"tie.eop"]
tnie.eop[(sims*7+1):(sims*8),"exclude"] <- exclude[1:sims,"medium"]
tnie.eop[(sims*8+1):(sims*9),"rep"] <- 1:sims
tnie.eop[(sims*8+1):(sims*9),"dgm"] <- "medium"
tnie.eop[(sims*8+1):(sims*9),"method"] <- "modal"
tnie.eop[(sims*8+1):(sims*9),"se"] <- modal.med[1:sims,"tie.eop.se"]
tnie.eop[(sims*8+1):(sims*9),"theta"] <- modal.med[1:sims,"tie.eop"]
tnie.eop[(sims*8+1):(sims*9),"exclude"] <- exclude[1:sims,"medium"]
tnie.eop[(sims*9+1):(sims*10),"rep"] <- 1:sims
tnie.eop[(sims*9+1):(sims*10),"dgm"] <- "medium"
tnie.eop[(sims*9+1):(sims*10),"method"] <- "npcd"
tnie.eop[(sims*9+1):(sims*10),"se"] <- npcd.med[1:sims,"tie.eop.se"]
tnie.eop[(sims*9+1):(sims*10),"theta"] <- npcd.med[1:sims,"tie.eop"]
tnie.eop[(sims*9+1):(sims*10),"exclude"] <- exclude[1:sims,"medium"]
tnie.eop[(sims*10+1):(sims*11),"rep"] <- 1:sims
tnie.eop[(sims*10+1):(sims*11),"dgm"] <- "medium"
tnie.eop[(sims*10+1):(sims*11),"method"] <- "incpcd"
tnie.eop[(sims*10+1):(sims*11),"se"] <- incpcd.med[1:sims,"tie.eop.se"]
tnie.eop[(sims*10+1):(sims*11),"theta"] <- incpcd.med[1:sims,"tie.eop"]
tnie.eop[(sims*10+1):(sims*11),"exclude"] <- exclude[1:sims,"medium"]
tnie.eop[(sims*11+1):(sims*12),"rep"] <- 1:sims
tnie.eop[(sims*11+1):(sims*12),"dgm"] <- "medium"
tnie.eop[(sims*11+1):(sims*12),"method"] <- "upcd"
tnie.eop[(sims*11+1):(sims*12),"se"] <- pcd.med[1:sims,"tie.eop.se"]
tnie.eop[(sims*11+1):(sims*12),"theta"] <- pcd.med[1:sims,"tie.eop"]
tnie.eop[(sims*11+1):(sims*12),"exclude"] <- exclude[1:sims,"medium"]
tnie.eop[(sims*12+1):(sims*13),"rep"] <- 1:sims
tnie.eop[(sims*12+1):(sims*13),"dgm"] <- "poor"
tnie.eop[(sims*12+1):(sims*13),"method"] <- "onestep"
tnie.eop[(sims*12+1):(sims*13),"se"] <- onestep.poor[1:sims,"tie.eop.se"]
tnie.eop[(sims*12+1):(sims*13),"theta"] <- onestep.poor[1:sims,"tie.eop"]
tnie.eop[(sims*12+1):(sims*13),"exclude"] <- exclude[1:sims,"poor"]
tnie.eop[(sims*13+1):(sims*14),"rep"] <- 1:sims
tnie.eop[(sims*13+1):(sims*14),"dgm"] <- "poor"
tnie.eop[(sims*13+1):(sims*14),"method"] <- "bch"
tnie.eop[(sims*13+1):(sims*14),"se"] <- bch.poor[1:sims,"tie.eop.se"]
tnie.eop[(sims*13+1):(sims*14),"theta"] <- bch.poor[1:sims,"tie.eop"]
tnie.eop[(sims*13+1):(sims*14),"exclude"] <- exclude[1:sims,"poor"]
tnie.eop[(sims*14+1):(sims*15),"rep"] <- 1:sims
tnie.eop[(sims*14+1):(sims*15),"dgm"] <- "poor"
tnie.eop[(sims*14+1):(sims*15),"method"] <- "modal"
tnie.eop[(sims*14+1):(sims*15),"se"] <- modal.poor[1:sims,"tie.eop.se"]
tnie.eop[(sims*14+1):(sims*15),"theta"] <- modal.poor[1:sims,"tie.eop"]
tnie.eop[(sims*14+1):(sims*15),"exclude"] <- exclude[1:sims,"poor"]
tnie.eop[(sims*15+1):(sims*16),"rep"] <- 1:sims
tnie.eop[(sims*15+1):(sims*16),"dgm"] <- "poor"
tnie.eop[(sims*15+1):(sims*16),"method"] <- "npcd"
tnie.eop[(sims*15+1):(sims*16),"se"] <- npcd.poor[1:sims,"tie.eop.se"]
tnie.eop[(sims*15+1):(sims*16),"theta"] <- npcd.poor[1:sims,"tie.eop"]
tnie.eop[(sims*15+1):(sims*16),"exclude"] <- exclude[1:sims,"poor"]
tnie.eop[(sims*16+1):(sims*17),"rep"] <- 1:sims
tnie.eop[(sims*16+1):(sims*17),"dgm"] <- "poor"
tnie.eop[(sims*16+1):(sims*17),"method"] <- "incpcd"
tnie.eop[(sims*16+1):(sims*17),"se"] <- incpcd.poor[1:sims,"tie.eop.se"]
tnie.eop[(sims*16+1):(sims*17),"theta"] <- incpcd.poor[1:sims,"tie.eop"]
tnie.eop[(sims*16+1):(sims*17),"exclude"] <- exclude[1:sims,"poor"]
tnie.eop[(sims*17+1):(sims*18),"rep"] <- 1:sims
tnie.eop[(sims*17+1):(sims*18),"dgm"] <- "poor"
tnie.eop[(sims*17+1):(sims*18),"method"] <- "upcd"
tnie.eop[(sims*17+1):(sims*18),"se"] <- pcd.poor[1:sims,"tie.eop.se"]
tnie.eop[(sims*17+1):(sims*18),"theta"] <- pcd.poor[1:sims,"tie.eop"]
tnie.eop[(sims*17+1):(sims*18),"exclude"] <- exclude[1:sims,"poor"]

#tnie effect for ao vs low
tnie.ao[1:sims,"rep"] <- 1:sims
tnie.ao[1:sims,"dgm"] <- "good"
tnie.ao[1:sims,"method"] <- "onestep"
tnie.ao[1:sims,"se"] <- onestep.good[1:sims,"tie.ao.se"]
tnie.ao[1:sims,"theta"] <- onestep.good[1:sims,"tie.ao"]
tnie.ao[1:sims,"exclude"] <- exclude[1:sims,"good"]
tnie.ao[(sims+1):(sims*2),"rep"] <- 1:sims
tnie.ao[(sims+1):(sims*2),"dgm"] <- "good"
tnie.ao[(sims+1):(sims*2),"method"] <- "bch"
tnie.ao[(sims+1):(sims*2),"se"] <- bch.good[1:sims,"tie.ao.se"]
tnie.ao[(sims+1):(sims*2),"theta"] <- bch.good[1:sims,"tie.ao"]
tnie.ao[(sims+1):(sims*2),"exclude"] <- exclude[1:sims,"good"]
tnie.ao[(sims*2+1):(sims*3),"rep"] <- 1:sims
tnie.ao[(sims*2+1):(sims*3),"dgm"] <- "good"
tnie.ao[(sims*2+1):(sims*3),"method"] <- "modal"
tnie.ao[(sims*2+1):(sims*3),"se"] <- modal.good[1:sims,"tie.ao.se"]
tnie.ao[(sims*2+1):(sims*3),"theta"] <- modal.good[1:sims,"tie.ao"]
tnie.ao[(sims*2+1):(sims*3),"exclude"] <- exclude[1:sims,"good"]
tnie.ao[(sims*3+1):(sims*4),"rep"] <- 1:sims
tnie.ao[(sims*3+1):(sims*4),"dgm"] <- "good"
tnie.ao[(sims*3+1):(sims*4),"method"] <- "npcd"
tnie.ao[(sims*3+1):(sims*4),"se"] <- npcd.good[1:sims,"tie.ao.se"]
tnie.ao[(sims*3+1):(sims*4),"theta"] <- npcd.good[1:sims,"tie.ao"]
tnie.ao[(sims*3+1):(sims*4),"exclude"] <- exclude[1:sims,"good"]
tnie.ao[(sims*4+1):(sims*5),"rep"] <- 1:sims
tnie.ao[(sims*4+1):(sims*5),"dgm"] <- "good"
tnie.ao[(sims*4+1):(sims*5),"method"] <- "incpcd"
tnie.ao[(sims*4+1):(sims*5),"se"] <- incpcd.good[1:sims,"tie.ao.se"]
tnie.ao[(sims*4+1):(sims*5),"theta"] <- incpcd.good[1:sims,"tie.ao"]
tnie.ao[(sims*4+1):(sims*5),"exclude"] <- exclude[1:sims,"good"]
tnie.ao[(sims*5+1):(sims*6),"rep"] <- 1:sims
tnie.ao[(sims*5+1):(sims*6),"dgm"] <- "good"
tnie.ao[(sims*5+1):(sims*6),"method"] <- "upcd"
tnie.ao[(sims*5+1):(sims*6),"se"] <- pcd.good[1:sims,"tie.ao.se"]
tnie.ao[(sims*5+1):(sims*6),"theta"] <- pcd.good[1:sims,"tie.ao"]
tnie.ao[(sims*5+1):(sims*6),"exclude"] <- exclude[1:sims,"good"]
tnie.ao[(sims*6+1):(sims*7),"rep"] <- 1:sims
tnie.ao[(sims*6+1):(sims*7),"dgm"] <- "medium"
tnie.ao[(sims*6+1):(sims*7),"method"] <- "onestep"
tnie.ao[(sims*6+1):(sims*7),"se"] <- onestep.med[1:sims,"tie.ao.se"]
tnie.ao[(sims*6+1):(sims*7),"theta"] <- onestep.med[1:sims,"tie.ao"]
tnie.ao[(sims*6+1):(sims*7),"exclude"] <- exclude[1:sims,"medium"]
tnie.ao[(sims*7+1):(sims*8),"rep"] <- 1:sims
tnie.ao[(sims*7+1):(sims*8),"dgm"] <- "medium"
tnie.ao[(sims*7+1):(sims*8),"method"] <- "bch"
tnie.ao[(sims*7+1):(sims*8),"se"] <- bch.med[1:sims,"tie.ao.se"]
tnie.ao[(sims*7+1):(sims*8),"theta"] <- bch.med[1:sims,"tie.ao"]
tnie.ao[(sims*7+1):(sims*8),"exclude"] <- exclude[1:sims,"medium"]
tnie.ao[(sims*8+1):(sims*9),"rep"] <- 1:sims
tnie.ao[(sims*8+1):(sims*9),"dgm"] <- "medium"
tnie.ao[(sims*8+1):(sims*9),"method"] <- "modal"
tnie.ao[(sims*8+1):(sims*9),"se"] <- modal.med[1:sims,"tie.ao.se"]
tnie.ao[(sims*8+1):(sims*9),"theta"] <- modal.med[1:sims,"tie.ao"]
tnie.ao[(sims*8+1):(sims*9),"exclude"] <- exclude[1:sims,"medium"]
tnie.ao[(sims*9+1):(sims*10),"rep"] <- 1:sims
tnie.ao[(sims*9+1):(sims*10),"dgm"] <- "medium"
tnie.ao[(sims*9+1):(sims*10),"method"] <- "npcd"
tnie.ao[(sims*9+1):(sims*10),"se"] <- npcd.med[1:sims,"tie.ao.se"]
tnie.ao[(sims*9+1):(sims*10),"theta"] <- npcd.med[1:sims,"tie.ao"]
tnie.ao[(sims*9+1):(sims*10),"exclude"] <- exclude[1:sims,"medium"]
tnie.ao[(sims*10+1):(sims*11),"rep"] <- 1:sims
tnie.ao[(sims*10+1):(sims*11),"dgm"] <- "medium"
tnie.ao[(sims*10+1):(sims*11),"method"] <- "incpcd"
tnie.ao[(sims*10+1):(sims*11),"se"] <- incpcd.med[1:sims,"tie.ao.se"]
tnie.ao[(sims*10+1):(sims*11),"theta"] <- incpcd.med[1:sims,"tie.ao"]
tnie.ao[(sims*10+1):(sims*11),"exclude"] <- exclude[1:sims,"medium"]
tnie.ao[(sims*11+1):(sims*12),"rep"] <- 1:sims
tnie.ao[(sims*11+1):(sims*12),"dgm"] <- "medium"
tnie.ao[(sims*11+1):(sims*12),"method"] <- "upcd"
tnie.ao[(sims*11+1):(sims*12),"se"] <- pcd.med[1:sims,"tie.ao.se"]
tnie.ao[(sims*11+1):(sims*12),"theta"] <- pcd.med[1:sims,"tie.ao"]
tnie.ao[(sims*11+1):(sims*12),"exclude"] <- exclude[1:sims,"medium"]
tnie.ao[(sims*12+1):(sims*13),"rep"] <- 1:sims
tnie.ao[(sims*12+1):(sims*13),"dgm"] <- "poor"
tnie.ao[(sims*12+1):(sims*13),"method"] <- "onestep"
tnie.ao[(sims*12+1):(sims*13),"se"] <- onestep.poor[1:sims,"tie.ao.se"]
tnie.ao[(sims*12+1):(sims*13),"theta"] <- onestep.poor[1:sims,"tie.ao"]
tnie.ao[(sims*12+1):(sims*13),"exclude"] <- exclude[1:sims,"poor"]
tnie.ao[(sims*13+1):(sims*14),"rep"] <- 1:sims
tnie.ao[(sims*13+1):(sims*14),"dgm"] <- "poor"
tnie.ao[(sims*13+1):(sims*14),"method"] <- "bch"
tnie.ao[(sims*13+1):(sims*14),"se"] <- bch.poor[1:sims,"tie.ao.se"]
tnie.ao[(sims*13+1):(sims*14),"theta"] <- bch.poor[1:sims,"tie.ao"]
tnie.ao[(sims*13+1):(sims*14),"exclude"] <- exclude[1:sims,"poor"]
tnie.ao[(sims*14+1):(sims*15),"rep"] <- 1:sims
tnie.ao[(sims*14+1):(sims*15),"dgm"] <- "poor"
tnie.ao[(sims*14+1):(sims*15),"method"] <- "modal"
tnie.ao[(sims*14+1):(sims*15),"se"] <- modal.poor[1:sims,"tie.ao.se"]
tnie.ao[(sims*14+1):(sims*15),"theta"] <- modal.poor[1:sims,"tie.ao"]
tnie.ao[(sims*14+1):(sims*15),"exclude"] <- exclude[1:sims,"poor"]
tnie.ao[(sims*15+1):(sims*16),"rep"] <- 1:sims
tnie.ao[(sims*15+1):(sims*16),"dgm"] <- "poor"
tnie.ao[(sims*15+1):(sims*16),"method"] <- "npcd"
tnie.ao[(sims*15+1):(sims*16),"se"] <- npcd.poor[1:sims,"tie.ao.se"]
tnie.ao[(sims*15+1):(sims*16),"theta"] <- npcd.poor[1:sims,"tie.ao"]
tnie.ao[(sims*15+1):(sims*16),"exclude"] <- exclude[1:sims,"poor"]
tnie.ao[(sims*16+1):(sims*17),"rep"] <- 1:sims
tnie.ao[(sims*16+1):(sims*17),"dgm"] <- "poor"
tnie.ao[(sims*16+1):(sims*17),"method"] <- "incpcd"
tnie.ao[(sims*16+1):(sims*17),"se"] <- incpcd.poor[1:sims,"tie.ao.se"]
tnie.ao[(sims*16+1):(sims*17),"theta"] <- incpcd.poor[1:sims,"tie.ao"]
tnie.ao[(sims*16+1):(sims*17),"exclude"] <- exclude[1:sims,"poor"]
tnie.ao[(sims*17+1):(sims*18),"rep"] <- 1:sims
tnie.ao[(sims*17+1):(sims*18),"dgm"] <- "poor"
tnie.ao[(sims*17+1):(sims*18),"method"] <- "upcd"
tnie.ao[(sims*17+1):(sims*18),"se"] <- pcd.poor[1:sims,"tie.ao.se"]
tnie.ao[(sims*17+1):(sims*18),"theta"] <- pcd.poor[1:sims,"tie.ao"]
tnie.ao[(sims*17+1):(sims*18),"exclude"] <- exclude[1:sims,"poor"]

#tnie effect for cl vs low
tnie.cl[1:sims,"rep"] <- 1:sims
tnie.cl[1:sims,"dgm"] <- "good"
tnie.cl[1:sims,"method"] <- "onestep"
tnie.cl[1:sims,"se"] <- onestep.good[1:sims,"tie.cl.se"]
tnie.cl[1:sims,"theta"] <- onestep.good[1:sims,"tie.cl"]
tnie.cl[1:sims,"exclude"] <- exclude[1:sims,"good"]
tnie.cl[(sims+1):(sims*2),"rep"] <- 1:sims
tnie.cl[(sims+1):(sims*2),"dgm"] <- "good"
tnie.cl[(sims+1):(sims*2),"method"] <- "bch"
tnie.cl[(sims+1):(sims*2),"se"] <- bch.good[1:sims,"tie.cl.se"]
tnie.cl[(sims+1):(sims*2),"theta"] <- bch.good[1:sims,"tie.cl"]
tnie.cl[(sims+1):(sims*2),"exclude"] <- exclude[1:sims,"good"]
tnie.cl[(sims*2+1):(sims*3),"rep"] <- 1:sims
tnie.cl[(sims*2+1):(sims*3),"dgm"] <- "good"
tnie.cl[(sims*2+1):(sims*3),"method"] <- "modal"
tnie.cl[(sims*2+1):(sims*3),"se"] <- modal.good[1:sims,"tie.cl.se"]
tnie.cl[(sims*2+1):(sims*3),"theta"] <- modal.good[1:sims,"tie.cl"]
tnie.cl[(sims*2+1):(sims*3),"exclude"] <- exclude[1:sims,"good"]
tnie.cl[(sims*3+1):(sims*4),"rep"] <- 1:sims
tnie.cl[(sims*3+1):(sims*4),"dgm"] <- "good"
tnie.cl[(sims*3+1):(sims*4),"method"] <- "npcd"
tnie.cl[(sims*3+1):(sims*4),"se"] <- npcd.good[1:sims,"tie.cl.se"]
tnie.cl[(sims*3+1):(sims*4),"theta"] <- npcd.good[1:sims,"tie.cl"]
tnie.cl[(sims*3+1):(sims*4),"exclude"] <- exclude[1:sims,"good"]
tnie.cl[(sims*4+1):(sims*5),"rep"] <- 1:sims
tnie.cl[(sims*4+1):(sims*5),"dgm"] <- "good"
tnie.cl[(sims*4+1):(sims*5),"method"] <- "incpcd"
tnie.cl[(sims*4+1):(sims*5),"se"] <- incpcd.good[1:sims,"tie.cl.se"]
tnie.cl[(sims*4+1):(sims*5),"theta"] <- incpcd.good[1:sims,"tie.cl"]
tnie.cl[(sims*4+1):(sims*5),"exclude"] <- exclude[1:sims,"good"]
tnie.cl[(sims*5+1):(sims*6),"rep"] <- 1:sims
tnie.cl[(sims*5+1):(sims*6),"dgm"] <- "good"
tnie.cl[(sims*5+1):(sims*6),"method"] <- "upcd"
tnie.cl[(sims*5+1):(sims*6),"se"] <- pcd.good[1:sims,"tie.cl.se"]
tnie.cl[(sims*5+1):(sims*6),"theta"] <- pcd.good[1:sims,"tie.cl"]
tnie.cl[(sims*5+1):(sims*6),"exclude"] <- exclude[1:sims,"good"]
tnie.cl[(sims*6+1):(sims*7),"rep"] <- 1:sims
tnie.cl[(sims*6+1):(sims*7),"dgm"] <- "medium"
tnie.cl[(sims*6+1):(sims*7),"method"] <- "onestep"
tnie.cl[(sims*6+1):(sims*7),"se"] <- onestep.med[1:sims,"tie.cl.se"]
tnie.cl[(sims*6+1):(sims*7),"theta"] <- onestep.med[1:sims,"tie.cl"]
tnie.cl[(sims*6+1):(sims*7),"exclude"] <- exclude[1:sims,"medium"]
tnie.cl[(sims*7+1):(sims*8),"rep"] <- 1:sims
tnie.cl[(sims*7+1):(sims*8),"dgm"] <- "medium"
tnie.cl[(sims*7+1):(sims*8),"method"] <- "bch"
tnie.cl[(sims*7+1):(sims*8),"se"] <- bch.med[1:sims,"tie.cl.se"]
tnie.cl[(sims*7+1):(sims*8),"theta"] <- bch.med[1:sims,"tie.cl"]
tnie.cl[(sims*7+1):(sims*8),"exclude"] <- exclude[1:sims,"medium"]
tnie.cl[(sims*8+1):(sims*9),"rep"] <- 1:sims
tnie.cl[(sims*8+1):(sims*9),"dgm"] <- "medium"
tnie.cl[(sims*8+1):(sims*9),"method"] <- "modal"
tnie.cl[(sims*8+1):(sims*9),"se"] <- modal.med[1:sims,"tie.cl.se"]
tnie.cl[(sims*8+1):(sims*9),"theta"] <- modal.med[1:sims,"tie.cl"]
tnie.cl[(sims*8+1):(sims*9),"exclude"] <- exclude[1:sims,"medium"]
tnie.cl[(sims*9+1):(sims*10),"rep"] <- 1:sims
tnie.cl[(sims*9+1):(sims*10),"dgm"] <- "medium"
tnie.cl[(sims*9+1):(sims*10),"method"] <- "npcd"
tnie.cl[(sims*9+1):(sims*10),"se"] <- npcd.med[1:sims,"tie.cl.se"]
tnie.cl[(sims*9+1):(sims*10),"theta"] <- npcd.med[1:sims,"tie.cl"]
tnie.cl[(sims*9+1):(sims*10),"exclude"] <- exclude[1:sims,"medium"]
tnie.cl[(sims*10+1):(sims*11),"rep"] <- 1:sims
tnie.cl[(sims*10+1):(sims*11),"dgm"] <- "medium"
tnie.cl[(sims*10+1):(sims*11),"method"] <- "incpcd"
tnie.cl[(sims*10+1):(sims*11),"se"] <- incpcd.med[1:sims,"tie.cl.se"]
tnie.cl[(sims*10+1):(sims*11),"theta"] <- incpcd.med[1:sims,"tie.cl"]
tnie.cl[(sims*10+1):(sims*11),"exclude"] <- exclude[1:sims,"medium"]
tnie.cl[(sims*11+1):(sims*12),"rep"] <- 1:sims
tnie.cl[(sims*11+1):(sims*12),"dgm"] <- "medium"
tnie.cl[(sims*11+1):(sims*12),"method"] <- "upcd"
tnie.cl[(sims*11+1):(sims*12),"se"] <- pcd.med[1:sims,"tie.cl.se"]
tnie.cl[(sims*11+1):(sims*12),"theta"] <- pcd.med[1:sims,"tie.cl"]
tnie.cl[(sims*11+1):(sims*12),"exclude"] <- exclude[1:sims,"medium"]
tnie.cl[(sims*12+1):(sims*13),"rep"] <- 1:sims
tnie.cl[(sims*12+1):(sims*13),"dgm"] <- "poor"
tnie.cl[(sims*12+1):(sims*13),"method"] <- "onestep"
tnie.cl[(sims*12+1):(sims*13),"se"] <- onestep.poor[1:sims,"tie.cl.se"]
tnie.cl[(sims*12+1):(sims*13),"theta"] <- onestep.poor[1:sims,"tie.cl"]
tnie.cl[(sims*12+1):(sims*13),"exclude"] <- exclude[1:sims,"poor"]
tnie.cl[(sims*13+1):(sims*14),"rep"] <- 1:sims
tnie.cl[(sims*13+1):(sims*14),"dgm"] <- "poor"
tnie.cl[(sims*13+1):(sims*14),"method"] <- "bch"
tnie.cl[(sims*13+1):(sims*14),"se"] <- bch.poor[1:sims,"tie.cl.se"]
tnie.cl[(sims*13+1):(sims*14),"theta"] <- bch.poor[1:sims,"tie.cl"]
tnie.cl[(sims*13+1):(sims*14),"exclude"] <- exclude[1:sims,"poor"]
tnie.cl[(sims*14+1):(sims*15),"rep"] <- 1:sims
tnie.cl[(sims*14+1):(sims*15),"dgm"] <- "poor"
tnie.cl[(sims*14+1):(sims*15),"method"] <- "modal"
tnie.cl[(sims*14+1):(sims*15),"se"] <- modal.poor[1:sims,"tie.cl.se"]
tnie.cl[(sims*14+1):(sims*15),"theta"] <- modal.poor[1:sims,"tie.cl"]
tnie.cl[(sims*14+1):(sims*15),"exclude"] <- exclude[1:sims,"poor"]
tnie.cl[(sims*15+1):(sims*16),"rep"] <- 1:sims
tnie.cl[(sims*15+1):(sims*16),"dgm"] <- "poor"
tnie.cl[(sims*15+1):(sims*16),"method"] <- "npcd"
tnie.cl[(sims*15+1):(sims*16),"se"] <- npcd.poor[1:sims,"tie.cl.se"]
tnie.cl[(sims*15+1):(sims*16),"theta"] <- npcd.poor[1:sims,"tie.cl"]
tnie.cl[(sims*15+1):(sims*16),"exclude"] <- exclude[1:sims,"poor"]
tnie.cl[(sims*16+1):(sims*17),"rep"] <- 1:sims
tnie.cl[(sims*16+1):(sims*17),"dgm"] <- "poor"
tnie.cl[(sims*16+1):(sims*17),"method"] <- "incpcd"
tnie.cl[(sims*16+1):(sims*17),"se"] <- incpcd.poor[1:sims,"tie.cl.se"]
tnie.cl[(sims*16+1):(sims*17),"theta"] <- incpcd.poor[1:sims,"tie.cl"]
tnie.cl[(sims*16+1):(sims*17),"exclude"] <- exclude[1:sims,"poor"]
tnie.cl[(sims*17+1):(sims*18),"rep"] <- 1:sims
tnie.cl[(sims*17+1):(sims*18),"dgm"] <- "poor"
tnie.cl[(sims*17+1):(sims*18),"method"] <- "upcd"
tnie.cl[(sims*17+1):(sims*18),"se"] <- pcd.poor[1:sims,"tie.cl.se"]
tnie.cl[(sims*17+1):(sims*18),"theta"] <- pcd.poor[1:sims,"tie.cl"]
tnie.cl[(sims*17+1):(sims*18),"exclude"] <- exclude[1:sims,"poor"]

#pnde effect for eop vs low
pnde.eop[1:sims,"rep"] <- 1:sims
pnde.eop[1:sims,"dgm"] <- "good"
pnde.eop[1:sims,"method"] <- "onestep"
pnde.eop[1:sims,"se"] <- onestep.good[1:sims,"pde.eop.se"]
pnde.eop[1:sims,"theta"] <- onestep.good[1:sims,"pde.eop"]
pnde.eop[1:sims,"exclude"] <- exclude[1:sims,"good"]
pnde.eop[(sims+1):(sims*2),"rep"] <- 1:sims
pnde.eop[(sims+1):(sims*2),"dgm"] <- "good"
pnde.eop[(sims+1):(sims*2),"method"] <- "bch"
pnde.eop[(sims+1):(sims*2),"se"] <- bch.good[1:sims,"pde.eop.se"]
pnde.eop[(sims+1):(sims*2),"theta"] <- bch.good[1:sims,"pde.eop"]
pnde.eop[(sims+1):(sims*2),"exclude"] <- exclude[1:sims,"good"]
pnde.eop[(sims*2+1):(sims*3),"rep"] <- 1:sims
pnde.eop[(sims*2+1):(sims*3),"dgm"] <- "good"
pnde.eop[(sims*2+1):(sims*3),"method"] <- "modal"
pnde.eop[(sims*2+1):(sims*3),"se"] <- modal.good[1:sims,"pde.eop.se"]
pnde.eop[(sims*2+1):(sims*3),"theta"] <- modal.good[1:sims,"pde.eop"]
pnde.eop[(sims*2+1):(sims*3),"exclude"] <- exclude[1:sims,"good"]
pnde.eop[(sims*3+1):(sims*4),"rep"] <- 1:sims
pnde.eop[(sims*3+1):(sims*4),"dgm"] <- "good"
pnde.eop[(sims*3+1):(sims*4),"method"] <- "npcd"
pnde.eop[(sims*3+1):(sims*4),"se"] <- npcd.good[1:sims,"pde.eop.se"]
pnde.eop[(sims*3+1):(sims*4),"theta"] <- npcd.good[1:sims,"pde.eop"]
pnde.eop[(sims*3+1):(sims*4),"exclude"] <- exclude[1:sims,"good"]
pnde.eop[(sims*4+1):(sims*5),"rep"] <- 1:sims
pnde.eop[(sims*4+1):(sims*5),"dgm"] <- "good"
pnde.eop[(sims*4+1):(sims*5),"method"] <- "incpcd"
pnde.eop[(sims*4+1):(sims*5),"se"] <- incpcd.good[1:sims,"pde.eop.se"]
pnde.eop[(sims*4+1):(sims*5),"theta"] <- incpcd.good[1:sims,"pde.eop"]
pnde.eop[(sims*4+1):(sims*5),"exclude"] <- exclude[1:sims,"good"]
pnde.eop[(sims*5+1):(sims*6),"rep"] <- 1:sims
pnde.eop[(sims*5+1):(sims*6),"dgm"] <- "good"
pnde.eop[(sims*5+1):(sims*6),"method"] <- "upcd"
pnde.eop[(sims*5+1):(sims*6),"se"] <- pcd.good[1:sims,"pde.eop.se"]
pnde.eop[(sims*5+1):(sims*6),"theta"] <- pcd.good[1:sims,"pde.eop"]
pnde.eop[(sims*5+1):(sims*6),"exclude"] <- exclude[1:sims,"good"]
pnde.eop[(sims*6+1):(sims*7),"rep"] <- 1:sims
pnde.eop[(sims*6+1):(sims*7),"dgm"] <- "medium"
pnde.eop[(sims*6+1):(sims*7),"method"] <- "onestep"
pnde.eop[(sims*6+1):(sims*7),"se"] <- onestep.med[1:sims,"pde.eop.se"]
pnde.eop[(sims*6+1):(sims*7),"theta"] <- onestep.med[1:sims,"pde.eop"]
pnde.eop[(sims*6+1):(sims*7),"exclude"] <- exclude[1:sims,"medium"]
pnde.eop[(sims*7+1):(sims*8),"rep"] <- 1:sims
pnde.eop[(sims*7+1):(sims*8),"dgm"] <- "medium"
pnde.eop[(sims*7+1):(sims*8),"method"] <- "bch"
pnde.eop[(sims*7+1):(sims*8),"se"] <- bch.med[1:sims,"pde.eop.se"]
pnde.eop[(sims*7+1):(sims*8),"theta"] <- bch.med[1:sims,"pde.eop"]
pnde.eop[(sims*7+1):(sims*8),"exclude"] <- exclude[1:sims,"medium"]
pnde.eop[(sims*8+1):(sims*9),"rep"] <- 1:sims
pnde.eop[(sims*8+1):(sims*9),"dgm"] <- "medium"
pnde.eop[(sims*8+1):(sims*9),"method"] <- "modal"
pnde.eop[(sims*8+1):(sims*9),"se"] <- modal.med[1:sims,"pde.eop.se"]
pnde.eop[(sims*8+1):(sims*9),"theta"] <- modal.med[1:sims,"pde.eop"]
pnde.eop[(sims*8+1):(sims*9),"exclude"] <- exclude[1:sims,"medium"]
pnde.eop[(sims*9+1):(sims*10),"rep"] <- 1:sims
pnde.eop[(sims*9+1):(sims*10),"dgm"] <- "medium"
pnde.eop[(sims*9+1):(sims*10),"method"] <- "npcd"
pnde.eop[(sims*9+1):(sims*10),"se"] <- npcd.med[1:sims,"pde.eop.se"]
pnde.eop[(sims*9+1):(sims*10),"theta"] <- npcd.med[1:sims,"pde.eop"]
pnde.eop[(sims*9+1):(sims*10),"exclude"] <- exclude[1:sims,"medium"]
pnde.eop[(sims*10+1):(sims*11),"rep"] <- 1:sims
pnde.eop[(sims*10+1):(sims*11),"dgm"] <- "medium"
pnde.eop[(sims*10+1):(sims*11),"method"] <- "incpcd"
pnde.eop[(sims*10+1):(sims*11),"se"] <- incpcd.med[1:sims,"pde.eop.se"]
pnde.eop[(sims*10+1):(sims*11),"theta"] <- incpcd.med[1:sims,"pde.eop"]
pnde.eop[(sims*10+1):(sims*11),"exclude"] <- exclude[1:sims,"medium"]
pnde.eop[(sims*11+1):(sims*12),"rep"] <- 1:sims
pnde.eop[(sims*11+1):(sims*12),"dgm"] <- "medium"
pnde.eop[(sims*11+1):(sims*12),"method"] <- "upcd"
pnde.eop[(sims*11+1):(sims*12),"se"] <- pcd.med[1:sims,"pde.eop.se"]
pnde.eop[(sims*11+1):(sims*12),"theta"] <- pcd.med[1:sims,"pde.eop"]
pnde.eop[(sims*11+1):(sims*12),"exclude"] <- exclude[1:sims,"medium"]
pnde.eop[(sims*12+1):(sims*13),"rep"] <- 1:sims
pnde.eop[(sims*12+1):(sims*13),"dgm"] <- "poor"
pnde.eop[(sims*12+1):(sims*13),"method"] <- "onestep"
pnde.eop[(sims*12+1):(sims*13),"se"] <- onestep.poor[1:sims,"pde.eop.se"]
pnde.eop[(sims*12+1):(sims*13),"theta"] <- onestep.poor[1:sims,"pde.eop"]
pnde.eop[(sims*12+1):(sims*13),"exclude"] <- exclude[1:sims,"poor"]
pnde.eop[(sims*13+1):(sims*14),"rep"] <- 1:sims
pnde.eop[(sims*13+1):(sims*14),"dgm"] <- "poor"
pnde.eop[(sims*13+1):(sims*14),"method"] <- "bch"
pnde.eop[(sims*13+1):(sims*14),"se"] <- bch.poor[1:sims,"pde.eop.se"]
pnde.eop[(sims*13+1):(sims*14),"theta"] <- bch.poor[1:sims,"pde.eop"]
pnde.eop[(sims*13+1):(sims*14),"exclude"] <- exclude[1:sims,"poor"]
pnde.eop[(sims*14+1):(sims*15),"rep"] <- 1:sims
pnde.eop[(sims*14+1):(sims*15),"dgm"] <- "poor"
pnde.eop[(sims*14+1):(sims*15),"method"] <- "modal"
pnde.eop[(sims*14+1):(sims*15),"se"] <- modal.poor[1:sims,"pde.eop.se"]
pnde.eop[(sims*14+1):(sims*15),"theta"] <- modal.poor[1:sims,"pde.eop"]
pnde.eop[(sims*14+1):(sims*15),"exclude"] <- exclude[1:sims,"poor"]
pnde.eop[(sims*15+1):(sims*16),"rep"] <- 1:sims
pnde.eop[(sims*15+1):(sims*16),"dgm"] <- "poor"
pnde.eop[(sims*15+1):(sims*16),"method"] <- "npcd"
pnde.eop[(sims*15+1):(sims*16),"se"] <- npcd.poor[1:sims,"pde.eop.se"]
pnde.eop[(sims*15+1):(sims*16),"theta"] <- npcd.poor[1:sims,"pde.eop"]
pnde.eop[(sims*15+1):(sims*16),"exclude"] <- exclude[1:sims,"poor"]
pnde.eop[(sims*16+1):(sims*17),"rep"] <- 1:sims
pnde.eop[(sims*16+1):(sims*17),"dgm"] <- "poor"
pnde.eop[(sims*16+1):(sims*17),"method"] <- "incpcd"
pnde.eop[(sims*16+1):(sims*17),"se"] <- incpcd.poor[1:sims,"pde.eop.se"]
pnde.eop[(sims*16+1):(sims*17),"theta"] <- incpcd.poor[1:sims,"pde.eop"]
pnde.eop[(sims*16+1):(sims*17),"exclude"] <- exclude[1:sims,"poor"]
pnde.eop[(sims*17+1):(sims*18),"rep"] <- 1:sims
pnde.eop[(sims*17+1):(sims*18),"dgm"] <- "poor"
pnde.eop[(sims*17+1):(sims*18),"method"] <- "upcd"
pnde.eop[(sims*17+1):(sims*18),"se"] <- pcd.poor[1:sims,"pde.eop.se"]
pnde.eop[(sims*17+1):(sims*18),"theta"] <- pcd.poor[1:sims,"pde.eop"]
pnde.eop[(sims*17+1):(sims*18),"exclude"] <- exclude[1:sims,"poor"]

#pnde effect for ao vs low
pnde.ao[1:sims,"rep"] <- 1:sims
pnde.ao[1:sims,"dgm"] <- "good"
pnde.ao[1:sims,"method"] <- "onestep"
pnde.ao[1:sims,"se"] <- onestep.good[1:sims,"pde.ao.se"]
pnde.ao[1:sims,"theta"] <- onestep.good[1:sims,"pde.ao"]
pnde.ao[1:sims,"exclude"] <- exclude[1:sims,"good"]
pnde.ao[(sims+1):(sims*2),"rep"] <- 1:sims
pnde.ao[(sims+1):(sims*2),"dgm"] <- "good"
pnde.ao[(sims+1):(sims*2),"method"] <- "bch"
pnde.ao[(sims+1):(sims*2),"se"] <- bch.good[1:sims,"pde.ao.se"]
pnde.ao[(sims+1):(sims*2),"theta"] <- bch.good[1:sims,"pde.ao"]
pnde.ao[(sims+1):(sims*2),"exclude"] <- exclude[1:sims,"good"]
pnde.ao[(sims*2+1):(sims*3),"rep"] <- 1:sims
pnde.ao[(sims*2+1):(sims*3),"dgm"] <- "good"
pnde.ao[(sims*2+1):(sims*3),"method"] <- "modal"
pnde.ao[(sims*2+1):(sims*3),"se"] <- modal.good[1:sims,"pde.ao.se"]
pnde.ao[(sims*2+1):(sims*3),"theta"] <- modal.good[1:sims,"pde.ao"]
pnde.ao[(sims*2+1):(sims*3),"exclude"] <- exclude[1:sims,"good"]
pnde.ao[(sims*3+1):(sims*4),"rep"] <- 1:sims
pnde.ao[(sims*3+1):(sims*4),"dgm"] <- "good"
pnde.ao[(sims*3+1):(sims*4),"method"] <- "npcd"
pnde.ao[(sims*3+1):(sims*4),"se"] <- npcd.good[1:sims,"pde.ao.se"]
pnde.ao[(sims*3+1):(sims*4),"theta"] <- npcd.good[1:sims,"pde.ao"]
pnde.ao[(sims*3+1):(sims*4),"exclude"] <- exclude[1:sims,"good"]
pnde.ao[(sims*4+1):(sims*5),"rep"] <- 1:sims
pnde.ao[(sims*4+1):(sims*5),"dgm"] <- "good"
pnde.ao[(sims*4+1):(sims*5),"method"] <- "incpcd"
pnde.ao[(sims*4+1):(sims*5),"se"] <- incpcd.good[1:sims,"pde.ao.se"]
pnde.ao[(sims*4+1):(sims*5),"theta"] <- incpcd.good[1:sims,"pde.ao"]
pnde.ao[(sims*4+1):(sims*5),"exclude"] <- exclude[1:sims,"good"]
pnde.ao[(sims*5+1):(sims*6),"rep"] <- 1:sims
pnde.ao[(sims*5+1):(sims*6),"dgm"] <- "good"
pnde.ao[(sims*5+1):(sims*6),"method"] <- "upcd"
pnde.ao[(sims*5+1):(sims*6),"se"] <- pcd.good[1:sims,"pde.ao.se"]
pnde.ao[(sims*5+1):(sims*6),"theta"] <- pcd.good[1:sims,"pde.ao"]
pnde.ao[(sims*5+1):(sims*6),"exclude"] <- exclude[1:sims,"good"]
pnde.ao[(sims*6+1):(sims*7),"rep"] <- 1:sims
pnde.ao[(sims*6+1):(sims*7),"dgm"] <- "medium"
pnde.ao[(sims*6+1):(sims*7),"method"] <- "onestep"
pnde.ao[(sims*6+1):(sims*7),"se"] <- onestep.med[1:sims,"pde.ao.se"]
pnde.ao[(sims*6+1):(sims*7),"theta"] <- onestep.med[1:sims,"pde.ao"]
pnde.ao[(sims*6+1):(sims*7),"exclude"] <- exclude[1:sims,"medium"]
pnde.ao[(sims*7+1):(sims*8),"rep"] <- 1:sims
pnde.ao[(sims*7+1):(sims*8),"dgm"] <- "medium"
pnde.ao[(sims*7+1):(sims*8),"method"] <- "bch"
pnde.ao[(sims*7+1):(sims*8),"se"] <- bch.med[1:sims,"pde.ao.se"]
pnde.ao[(sims*7+1):(sims*8),"theta"] <- bch.med[1:sims,"pde.ao"]
pnde.ao[(sims*7+1):(sims*8),"exclude"] <- exclude[1:sims,"medium"]
pnde.ao[(sims*8+1):(sims*9),"rep"] <- 1:sims
pnde.ao[(sims*8+1):(sims*9),"dgm"] <- "medium"
pnde.ao[(sims*8+1):(sims*9),"method"] <- "modal"
pnde.ao[(sims*8+1):(sims*9),"se"] <- modal.med[1:sims,"pde.ao.se"]
pnde.ao[(sims*8+1):(sims*9),"theta"] <- modal.med[1:sims,"pde.ao"]
pnde.ao[(sims*8+1):(sims*9),"exclude"] <- exclude[1:sims,"medium"]
pnde.ao[(sims*9+1):(sims*10),"rep"] <- 1:sims
pnde.ao[(sims*9+1):(sims*10),"dgm"] <- "medium"
pnde.ao[(sims*9+1):(sims*10),"method"] <- "npcd"
pnde.ao[(sims*9+1):(sims*10),"se"] <- npcd.med[1:sims,"pde.ao.se"]
pnde.ao[(sims*9+1):(sims*10),"theta"] <- npcd.med[1:sims,"pde.ao"]
pnde.ao[(sims*9+1):(sims*10),"exclude"] <- exclude[1:sims,"medium"]
pnde.ao[(sims*10+1):(sims*11),"rep"] <- 1:sims
pnde.ao[(sims*10+1):(sims*11),"dgm"] <- "medium"
pnde.ao[(sims*10+1):(sims*11),"method"] <- "incpcd"
pnde.ao[(sims*10+1):(sims*11),"se"] <- incpcd.med[1:sims,"pde.ao.se"]
pnde.ao[(sims*10+1):(sims*11),"theta"] <- incpcd.med[1:sims,"pde.ao"]
pnde.ao[(sims*10+1):(sims*11),"exclude"] <- exclude[1:sims,"medium"]
pnde.ao[(sims*11+1):(sims*12),"rep"] <- 1:sims
pnde.ao[(sims*11+1):(sims*12),"dgm"] <- "medium"
pnde.ao[(sims*11+1):(sims*12),"method"] <- "upcd"
pnde.ao[(sims*11+1):(sims*12),"se"] <- pcd.med[1:sims,"pde.ao.se"]
pnde.ao[(sims*11+1):(sims*12),"theta"] <- pcd.med[1:sims,"pde.ao"]
pnde.ao[(sims*11+1):(sims*12),"exclude"] <- exclude[1:sims,"medium"]
pnde.ao[(sims*12+1):(sims*13),"rep"] <- 1:sims
pnde.ao[(sims*12+1):(sims*13),"dgm"] <- "poor"
pnde.ao[(sims*12+1):(sims*13),"method"] <- "onestep"
pnde.ao[(sims*12+1):(sims*13),"se"] <- onestep.poor[1:sims,"pde.ao.se"]
pnde.ao[(sims*12+1):(sims*13),"theta"] <- onestep.poor[1:sims,"pde.ao"]
pnde.ao[(sims*12+1):(sims*13),"exclude"] <- exclude[1:sims,"poor"]
pnde.ao[(sims*13+1):(sims*14),"rep"] <- 1:sims
pnde.ao[(sims*13+1):(sims*14),"dgm"] <- "poor"
pnde.ao[(sims*13+1):(sims*14),"method"] <- "bch"
pnde.ao[(sims*13+1):(sims*14),"se"] <- bch.poor[1:sims,"pde.ao.se"]
pnde.ao[(sims*13+1):(sims*14),"theta"] <- bch.poor[1:sims,"pde.ao"]
pnde.ao[(sims*13+1):(sims*14),"exclude"] <- exclude[1:sims,"poor"]
pnde.ao[(sims*14+1):(sims*15),"rep"] <- 1:sims
pnde.ao[(sims*14+1):(sims*15),"dgm"] <- "poor"
pnde.ao[(sims*14+1):(sims*15),"method"] <- "modal"
pnde.ao[(sims*14+1):(sims*15),"se"] <- modal.poor[1:sims,"pde.ao.se"]
pnde.ao[(sims*14+1):(sims*15),"theta"] <- modal.poor[1:sims,"pde.ao"]
pnde.ao[(sims*14+1):(sims*15),"exclude"] <- exclude[1:sims,"poor"]
pnde.ao[(sims*15+1):(sims*16),"rep"] <- 1:sims
pnde.ao[(sims*15+1):(sims*16),"dgm"] <- "poor"
pnde.ao[(sims*15+1):(sims*16),"method"] <- "npcd"
pnde.ao[(sims*15+1):(sims*16),"se"] <- npcd.poor[1:sims,"pde.ao.se"]
pnde.ao[(sims*15+1):(sims*16),"theta"] <- npcd.poor[1:sims,"pde.ao"]
pnde.ao[(sims*15+1):(sims*16),"exclude"] <- exclude[1:sims,"poor"]
pnde.ao[(sims*16+1):(sims*17),"rep"] <- 1:sims
pnde.ao[(sims*16+1):(sims*17),"dgm"] <- "poor"
pnde.ao[(sims*16+1):(sims*17),"method"] <- "incpcd"
pnde.ao[(sims*16+1):(sims*17),"se"] <- incpcd.poor[1:sims,"pde.ao.se"]
pnde.ao[(sims*16+1):(sims*17),"theta"] <- incpcd.poor[1:sims,"pde.ao"]
pnde.ao[(sims*16+1):(sims*17),"exclude"] <- exclude[1:sims,"poor"]
pnde.ao[(sims*17+1):(sims*18),"rep"] <- 1:sims
pnde.ao[(sims*17+1):(sims*18),"dgm"] <- "poor"
pnde.ao[(sims*17+1):(sims*18),"method"] <- "upcd"
pnde.ao[(sims*17+1):(sims*18),"se"] <- pcd.poor[1:sims,"pde.ao.se"]
pnde.ao[(sims*17+1):(sims*18),"theta"] <- pcd.poor[1:sims,"pde.ao"]
pnde.ao[(sims*17+1):(sims*18),"exclude"] <- exclude[1:sims,"poor"]

#pnde effect for cl vs low
pnde.cl[1:sims,"rep"] <- 1:sims
pnde.cl[1:sims,"dgm"] <- "good"
pnde.cl[1:sims,"method"] <- "onestep"
pnde.cl[1:sims,"se"] <- onestep.good[1:sims,"pde.cl.se"]
pnde.cl[1:sims,"theta"] <- onestep.good[1:sims,"pde.cl"]
pnde.cl[1:sims,"exclude"] <- exclude[1:sims,"good"]
pnde.cl[(sims+1):(sims*2),"rep"] <- 1:sims
pnde.cl[(sims+1):(sims*2),"dgm"] <- "good"
pnde.cl[(sims+1):(sims*2),"method"] <- "bch"
pnde.cl[(sims+1):(sims*2),"se"] <- bch.good[1:sims,"pde.cl.se"]
pnde.cl[(sims+1):(sims*2),"theta"] <- bch.good[1:sims,"pde.cl"]
pnde.cl[(sims+1):(sims*2),"exclude"] <- exclude[1:sims,"good"]
pnde.cl[(sims*2+1):(sims*3),"rep"] <- 1:sims
pnde.cl[(sims*2+1):(sims*3),"dgm"] <- "good"
pnde.cl[(sims*2+1):(sims*3),"method"] <- "modal"
pnde.cl[(sims*2+1):(sims*3),"se"] <- modal.good[1:sims,"pde.cl.se"]
pnde.cl[(sims*2+1):(sims*3),"theta"] <- modal.good[1:sims,"pde.cl"]
pnde.cl[(sims*2+1):(sims*3),"exclude"] <- exclude[1:sims,"good"]
pnde.cl[(sims*3+1):(sims*4),"rep"] <- 1:sims
pnde.cl[(sims*3+1):(sims*4),"dgm"] <- "good"
pnde.cl[(sims*3+1):(sims*4),"method"] <- "npcd"
pnde.cl[(sims*3+1):(sims*4),"se"] <- npcd.good[1:sims,"pde.cl.se"]
pnde.cl[(sims*3+1):(sims*4),"theta"] <- npcd.good[1:sims,"pde.cl"]
pnde.cl[(sims*3+1):(sims*4),"exclude"] <- exclude[1:sims,"good"]
pnde.cl[(sims*4+1):(sims*5),"rep"] <- 1:sims
pnde.cl[(sims*4+1):(sims*5),"dgm"] <- "good"
pnde.cl[(sims*4+1):(sims*5),"method"] <- "incpcd"
pnde.cl[(sims*4+1):(sims*5),"se"] <- incpcd.good[1:sims,"pde.cl.se"]
pnde.cl[(sims*4+1):(sims*5),"theta"] <- incpcd.good[1:sims,"pde.cl"]
pnde.cl[(sims*4+1):(sims*5),"exclude"] <- exclude[1:sims,"good"]
pnde.cl[(sims*5+1):(sims*6),"rep"] <- 1:sims
pnde.cl[(sims*5+1):(sims*6),"dgm"] <- "good"
pnde.cl[(sims*5+1):(sims*6),"method"] <- "upcd"
pnde.cl[(sims*5+1):(sims*6),"se"] <- pcd.good[1:sims,"pde.cl.se"]
pnde.cl[(sims*5+1):(sims*6),"theta"] <- pcd.good[1:sims,"pde.cl"]
pnde.cl[(sims*5+1):(sims*6),"exclude"] <- exclude[1:sims,"good"]
pnde.cl[(sims*6+1):(sims*7),"rep"] <- 1:sims
pnde.cl[(sims*6+1):(sims*7),"dgm"] <- "medium"
pnde.cl[(sims*6+1):(sims*7),"method"] <- "onestep"
pnde.cl[(sims*6+1):(sims*7),"se"] <- onestep.med[1:sims,"pde.cl.se"]
pnde.cl[(sims*6+1):(sims*7),"theta"] <- onestep.med[1:sims,"pde.cl"]
pnde.cl[(sims*6+1):(sims*7),"exclude"] <- exclude[1:sims,"medium"]
pnde.cl[(sims*7+1):(sims*8),"rep"] <- 1:sims
pnde.cl[(sims*7+1):(sims*8),"dgm"] <- "medium"
pnde.cl[(sims*7+1):(sims*8),"method"] <- "bch"
pnde.cl[(sims*7+1):(sims*8),"se"] <- bch.med[1:sims,"pde.cl.se"]
pnde.cl[(sims*7+1):(sims*8),"theta"] <- bch.med[1:sims,"pde.cl"]
pnde.cl[(sims*7+1):(sims*8),"exclude"] <- exclude[1:sims,"medium"]
pnde.cl[(sims*8+1):(sims*9),"rep"] <- 1:sims
pnde.cl[(sims*8+1):(sims*9),"dgm"] <- "medium"
pnde.cl[(sims*8+1):(sims*9),"method"] <- "modal"
pnde.cl[(sims*8+1):(sims*9),"se"] <- modal.med[1:sims,"pde.cl.se"]
pnde.cl[(sims*8+1):(sims*9),"theta"] <- modal.med[1:sims,"pde.cl"]
pnde.cl[(sims*8+1):(sims*9),"exclude"] <- exclude[1:sims,"medium"]
pnde.cl[(sims*9+1):(sims*10),"rep"] <- 1:sims
pnde.cl[(sims*9+1):(sims*10),"dgm"] <- "medium"
pnde.cl[(sims*9+1):(sims*10),"method"] <- "npcd"
pnde.cl[(sims*9+1):(sims*10),"se"] <- npcd.med[1:sims,"pde.cl.se"]
pnde.cl[(sims*9+1):(sims*10),"theta"] <- npcd.med[1:sims,"pde.cl"]
pnde.cl[(sims*9+1):(sims*10),"exclude"] <- exclude[1:sims,"medium"]
pnde.cl[(sims*10+1):(sims*11),"rep"] <- 1:sims
pnde.cl[(sims*10+1):(sims*11),"dgm"] <- "medium"
pnde.cl[(sims*10+1):(sims*11),"method"] <- "incpcd"
pnde.cl[(sims*10+1):(sims*11),"se"] <- incpcd.med[1:sims,"pde.cl.se"]
pnde.cl[(sims*10+1):(sims*11),"theta"] <- incpcd.med[1:sims,"pde.cl"]
pnde.cl[(sims*10+1):(sims*11),"exclude"] <- exclude[1:sims,"medium"]
pnde.cl[(sims*11+1):(sims*12),"rep"] <- 1:sims
pnde.cl[(sims*11+1):(sims*12),"dgm"] <- "medium"
pnde.cl[(sims*11+1):(sims*12),"method"] <- "upcd"
pnde.cl[(sims*11+1):(sims*12),"se"] <- pcd.med[1:sims,"pde.cl.se"]
pnde.cl[(sims*11+1):(sims*12),"theta"] <- pcd.med[1:sims,"pde.cl"]
pnde.cl[(sims*11+1):(sims*12),"exclude"] <- exclude[1:sims,"medium"]
pnde.cl[(sims*12+1):(sims*13),"rep"] <- 1:sims
pnde.cl[(sims*12+1):(sims*13),"dgm"] <- "poor"
pnde.cl[(sims*12+1):(sims*13),"method"] <- "onestep"
pnde.cl[(sims*12+1):(sims*13),"se"] <- onestep.poor[1:sims,"pde.cl.se"]
pnde.cl[(sims*12+1):(sims*13),"theta"] <- onestep.poor[1:sims,"pde.cl"]
pnde.cl[(sims*12+1):(sims*13),"exclude"] <- exclude[1:sims,"poor"]
pnde.cl[(sims*13+1):(sims*14),"rep"] <- 1:sims
pnde.cl[(sims*13+1):(sims*14),"dgm"] <- "poor"
pnde.cl[(sims*13+1):(sims*14),"method"] <- "bch"
pnde.cl[(sims*13+1):(sims*14),"se"] <- bch.poor[1:sims,"pde.cl.se"]
pnde.cl[(sims*13+1):(sims*14),"theta"] <- bch.poor[1:sims,"pde.cl"]
pnde.cl[(sims*13+1):(sims*14),"exclude"] <- exclude[1:sims,"poor"]
pnde.cl[(sims*14+1):(sims*15),"rep"] <- 1:sims
pnde.cl[(sims*14+1):(sims*15),"dgm"] <- "poor"
pnde.cl[(sims*14+1):(sims*15),"method"] <- "modal"
pnde.cl[(sims*14+1):(sims*15),"se"] <- modal.poor[1:sims,"pde.cl.se"]
pnde.cl[(sims*14+1):(sims*15),"theta"] <- modal.poor[1:sims,"pde.cl"]
pnde.cl[(sims*14+1):(sims*15),"exclude"] <- exclude[1:sims,"poor"]
pnde.cl[(sims*15+1):(sims*16),"rep"] <- 1:sims
pnde.cl[(sims*15+1):(sims*16),"dgm"] <- "poor"
pnde.cl[(sims*15+1):(sims*16),"method"] <- "npcd"
pnde.cl[(sims*15+1):(sims*16),"se"] <- npcd.poor[1:sims,"pde.cl.se"]
pnde.cl[(sims*15+1):(sims*16),"theta"] <- npcd.poor[1:sims,"pde.cl"]
pnde.cl[(sims*15+1):(sims*16),"exclude"] <- exclude[1:sims,"poor"]
pnde.cl[(sims*16+1):(sims*17),"rep"] <- 1:sims
pnde.cl[(sims*16+1):(sims*17),"dgm"] <- "poor"
pnde.cl[(sims*16+1):(sims*17),"method"] <- "incpcd"
pnde.cl[(sims*16+1):(sims*17),"se"] <- incpcd.poor[1:sims,"pde.cl.se"]
pnde.cl[(sims*16+1):(sims*17),"theta"] <- incpcd.poor[1:sims,"pde.cl"]
pnde.cl[(sims*16+1):(sims*17),"exclude"] <- exclude[1:sims,"poor"]
pnde.cl[(sims*17+1):(sims*18),"rep"] <- 1:sims
pnde.cl[(sims*17+1):(sims*18),"dgm"] <- "poor"
pnde.cl[(sims*17+1):(sims*18),"method"] <- "upcd"
pnde.cl[(sims*17+1):(sims*18),"se"] <- pcd.poor[1:sims,"pde.cl.se"]
pnde.cl[(sims*17+1):(sims*18),"theta"] <- pcd.poor[1:sims,"pde.cl"]
pnde.cl[(sims*17+1):(sims*18),"exclude"] <- exclude[1:sims,"poor"]

#explore estimates data using summary, but be sure to do this separately for each DGM and method

#using total effect for eop vs low as an example- ALL 500 SIMS
################################################################################################

#summary statistics for theta
te.eop %>% 
  group_by(dgm, method) %>%
  summarise(
    Min = min(theta),
    Q1 = quantile(theta, 0.25),
    Median = median(theta),
    Mean = mean(theta),
    Q3 = quantile(theta, 0.75),
    Max = max(theta)
  )

tnie.eop %>% 
  group_by(dgm, method) %>%
  summarise(
    Min = min(theta),
    Q1 = quantile(theta, 0.25),
    Median = median(theta),
    Mean = mean(theta),
    Q3 = quantile(theta, 0.75),
    Max = max(theta)
  )

pnde.eop %>% 
  group_by(dgm, method) %>%
  summarise(
    Min = min(theta),
    Q1 = quantile(theta, 0.25),
    Median = median(theta),
    Mean = mean(theta),
    Q3 = quantile(theta, 0.75),
    Max = max(theta)
  )

#summary statistics for se:
te.eop %>% 
  group_by(dgm, method) %>%
  summarise(
    Min = min(se),
    Q1 = quantile(se, 0.25),
    Median = median(se),
    Mean = mean(se),
    Q3 = quantile(se, 0.75),
    Max = max(se)
  ) 

tnie.eop %>% 
  group_by(dgm, method) %>%
  summarise(
    Min = min(se),
    Q1 = quantile(se, 0.25),
    Median = median(se),
    Mean = mean(se),
    Q3 = quantile(se, 0.75),
    Max = max(se)
  ) 

pnde.eop %>% 
  group_by(dgm, method) %>%
  summarise(
    Min = min(se),
    Q1 = quantile(se, 0.25),
    Median = median(se),
    Mean = mean(se),
    Q3 = quantile(se, 0.75),
    Max = max(se)
  ) 

#there are a few sims that have a very large SE for onestep and upcds - this might be when there is a threshold with large SE in unconditional model

#Check that diff and se are complete
te.eop %>%
  group_by(dgm, method) %>%
  summarise(
    N.theta = sum(!is.na(theta)),
    P.theta = percent(mean(!is.na(theta))),
    N.se = sum(!is.na(se)),
    P.se = percent(mean(!is.na(se)))
  )

tnie.eop %>%
  group_by(dgm, method) %>%
  summarise(
    N.theta = sum(!is.na(theta)),
    P.theta = percent(mean(!is.na(theta))),
    N.se = sum(!is.na(se)),
    P.se = percent(mean(!is.na(se)))
  )

pnde.eop %>%
  group_by(dgm, method) %>%
  summarise(
    N.theta = sum(!is.na(theta)),
    P.theta = percent(mean(!is.na(theta))),
    N.se = sum(!is.na(se)),
    P.se = percent(mean(!is.na(se)))
  )

#total effect for ao vs low
################################################################################################

#summary statistics for theta
te.ao %>% 
  group_by(dgm, method) %>%
  summarise(
    Min = min(theta),
    Q1 = quantile(theta, 0.25),
    Median = median(theta),
    Mean = mean(theta),
    Q3 = quantile(theta, 0.75),
    Max = max(theta)
  )

#summary statistics for se:
te.ao %>%
  group_by(dgm, method) %>%
  summarise(
    Min = min(se),
    Q1 = quantile(se, 0.25),
    Median = median(se),
    Mean = mean(se),
    Q3 = quantile(se, 0.75),
    Max = max(se)
  )

#Check that diff and se are complete
te.ao %>%
  group_by(dgm, method) %>%
  summarise(
    N.theta = sum(!is.na(theta)),
    P.theta = percent(mean(!is.na(theta))),
    N.se = sum(!is.na(se)),
    P.se = percent(mean(!is.na(se)))
  )

#total effect for cl vs low
################################################################################################

#summary statistics for theta
te.cl %>% 
  group_by(dgm, method) %>%
  summarise(
    Min = min(theta),
    Q1 = quantile(theta, 0.25),
    Median = median(theta),
    Mean = mean(theta),
    Q3 = quantile(theta, 0.75),
    Max = max(theta)
  )

#summary statistics for se:
te.cl %>%
  group_by(dgm, method) %>%
  summarise(
    Min = min(se),
    Q1 = quantile(se, 0.25),
    Median = median(se),
    Mean = mean(se),
    Q3 = quantile(se, 0.75),
    Max = max(se)
  )

#Check that diff and se are complete
te.ao %>%
  group_by(dgm, method) %>%
  summarise(
    N.theta = sum(!is.na(theta)),
    P.theta = percent(mean(!is.na(theta))),
    N.se = sum(!is.na(se)),
    P.se = percent(mean(!is.na(se)))
  )

#Produce some plots of the results: e.g. histograms and scatterplots of SE against point estimate

#histograms - each resemble a normal distribution
ggplot(te.eop, aes(x = theta)) +
  geom_histogram(binwidth = 0.01) +
  theme_bw() +
  facet_grid(method ~ dgm, label = label_both) +
  theme(text = element_text(size = 9)) 

ggplot(tnie.eop, aes(x = theta)) +
  geom_histogram(binwidth = 0.01) +
  theme_bw() +
  facet_grid(method ~ dgm, label = label_both) +
  theme(text = element_text(size = 9)) 

ggplot(pnde.eop, aes(x = theta)) +
  geom_histogram(binwidth = 0.01) +
  theme_bw() +
  facet_grid(method ~ dgm, label = label_both) +
  theme(text = element_text(size = 9)) 


##############################
#APPLYING EXCLUSION CRITERIA #examining results for sims after applying exclusion criteria based on within-class threshold SE
##############################

#CHECKS ON ESTIMATES

#using package "forcats"
method.labs <- c(onestep = "onestep", bch = "bch", modal = "modal", npcd = "npcd", incpcd = "incpcd", upcd = "upcd")

#EOP TE removing exclusions
ggplot(subset(te.eop, exclude==0), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12)) +
  geom_point(size=0.1)

#EOP TE changing colour of exclusions
ggplot(subset(te.eop), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12)) +
  geom_point(size=0.1) +
  geom_point(aes(colour=exclude), show.legend = FALSE) +
  xlab("log risk-ratio for TE of EOP versus Low") +
  ylab("standard error") + 
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  geom_vline(xintercept=0.357, linetype="dashed", color = "red")

#EOP TNIE removing exclusions
ggplot(subset(tnie.eop, exclude==0), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12))

#EOP TNIE changing colour of exclusions
ggplot(subset(tnie.eop), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12)) +
  geom_point(size=0.1) +
  geom_point(aes(colour=exclude), show.legend = FALSE) +
  xlab("log risk-ratio for TNIE of EOP versus Low") +
  ylab("standard error") + 
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  geom_vline(xintercept=0.057, linetype="dashed", color = "red")

#EOP PNDE removing exclusions
ggplot(subset(pnde.eop, exclude==0), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12))

#EOP PNDE changing colour of exclusions
ggplot(subset(pnde.eop), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12)) +
  geom_point(size=0.1) +
  geom_point(aes(colour=exclude), show.legend = FALSE) +
  xlab("log risk-ratio for PNDE of EOP versus Low") +
  ylab("standard error") + 
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  geom_vline(xintercept=0.300, linetype="dashed", color = "red")

#AO TE  removing exclusions
ggplot(subset(te.ao, exclude==0), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12))

#AO TE changing colour of exclusions
ggplot(subset(te.ao), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12)) +
  geom_point(size=0.1) +
  geom_point(aes(colour=exclude), show.legend = FALSE) +
  xlab("log risk-ratio for TE of AO versus Low") +
  ylab("standard error") + 
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  geom_vline(xintercept=0.251, linetype="dashed", color = "red")

#AO TNIE  removing exclusions
ggplot(subset(tnie.ao, exclude==0), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12))

#AO TNIE changing colour of exclusions
ggplot(subset(tnie.ao), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12)) +
  geom_point(size=0.1) +
  geom_point(aes(colour=exclude), show.legend = FALSE) +
  xlab("log risk-ratio for TNIE of AO versus Low") +
  ylab("standard error") + 
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  geom_vline(xintercept=0.081, linetype="dashed", color = "red")

#AO PNDE  removing exclusions
ggplot(subset(pnde.ao, exclude==0), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12))

#AO PNDE changing colour of exclusions
ggplot(subset(pnde.ao), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 10)) +
  geom_point(size=0.1) +
  geom_point(aes(colour=exclude), show.legend = FALSE) +
  xlab("log risk-ratio for PNDE of AO versus Low") +
  ylab("standard error") + 
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  geom_vline(xintercept=0.170, linetype="dashed", color = "red")

#CL TE  removing exclusions
ggplot(subset(te.cl, exclude==0), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12))

#CL TE changing colour of exclusions
ggplot(subset(te.cl), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12)) +
  geom_point(size=0.1) +
  geom_point(aes(colour=exclude), show.legend = FALSE) +
  xlab("log risk-ratio for TE of CL versus Low") +
  ylab("standard error") + 
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  geom_vline(xintercept=0.134, linetype="dashed", color = "red")

#CL TNIE  removing exclusions
ggplot(subset(tnie.cl, exclude==0), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12))

#CL TNIE changing colour of exclusions
ggplot(subset(tnie.cl), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 10)) +
  geom_point(size=0.1) +
  geom_point(aes(colour=exclude), show.legend = FALSE) +
  xlab("log risk-ratio for TNIE of CL versus Low") +
  ylab("standard error") + 
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  geom_vline(xintercept=0.012, linetype="dashed", color = "red")

#CL PNDE  removing exclusions
ggplot(subset(pnde.cl, exclude==0), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12))

#CL PNDE changing colour of exclusions
ggplot(subset(pnde.cl), aes(x = theta, y = se)) +
  geom_point() +
  theme_bw() +
  facet_grid(fct_inorder(method) ~ dgm, labeller = labeller(method = method.labs)) +
  theme(text = element_text(size = 12)) +
  geom_point(size=0.1) +
  geom_point(aes(colour=exclude), show.legend = FALSE) +
  xlab("log risk-ratio for PNDE of CL versus Low") +
  ylab("standard error") + 
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  geom_vline(xintercept=0.122, linetype="dashed", color = "red")

#scatter plots comparing onestep and upcds
#dropping sims with SE > 2 times the average standard error (for within-class threshold)
te.eop.subset <- subset(te.eop, (te.eop[,"method"]=="onestep" | te.eop[,"method"]=="upcd") & te.eop[,"exclude"]==0)
s1.subset <- simsum(data = te.eop.subset, estvarname = "theta", se = "se", true = 0.356674944, methodvar = "method", by = "dgm", ref = "onestep", x = TRUE)

autoplot(s1.subset, type = "est") +
  ggplot2::theme_bw()
autoplot(s1.subset, type = "se") +
  ggplot2::theme_bw()

#scatter plots comparing upcd and bch
#dropping sims with SE > 2 times the average standard error (for within-class threshold)
te.eop.subset <- subset(te.eop, (te.eop[,"method"]=="bch" | te.eop[,"method"]=="upcd") & te.eop[,"exclude"]==0)
s1.subset <- simsum(data = te.eop.subset, estvarname = "theta", se = "se", true = 0.356674944, methodvar = "method", by = "dgm", x = TRUE)

autoplot(s1.subset, type = "est") +
  ggplot2::theme_bw()
autoplot(s1.subset, type = "se") +
  ggplot2::theme_bw()

#TABLES

te.eop.subset <- subset(te.eop, (te.eop[,"exclude"]==0))
tnie.eop.subset <- subset(tnie.eop, (tnie.eop[,"exclude"]==0))
pnde.eop.subset <- subset(pnde.eop, (pnde.eop[,"exclude"]==0))
te.ao.subset <- subset(te.ao, (te.ao[,"exclude"]==0))
tnie.ao.subset <- subset(tnie.ao, (tnie.ao[,"exclude"]==0))
pnde.ao.subset <- subset(pnde.ao, (pnde.ao[,"exclude"]==0))
te.cl.subset <- subset(te.cl, (te.cl[,"exclude"]==0))
tnie.cl.subset <- subset(tnie.cl, (tnie.cl[,"exclude"]==0))
pnde.cl.subset <- subset(pnde.cl, (pnde.cl[,"exclude"]==0))

te.eop.subset %>%
  simsum(estvarname = "theta", se = "se", methodvar = "method", true = 0.356674944, by = "dgm", ref = "onestep", x = TRUE) %>%
  summary() %>%
  tidy() %>%
  kable()

tnie.eop.subset %>%
  simsum(estvarname = "theta", se = "se", methodvar = "method", true = 0.05702914, by = "dgm", ref = "onestep", x = TRUE) %>%
  summary() %>%
  tidy() %>%
  kable()

pnde.eop.subset %>%
  simsum(estvarname = "theta", se = "se", methodvar = "method", true = 0.299645804, by = "dgm", ref = "onestep", x = TRUE) %>%
  summary() %>%
  tidy() %>%
  kable()

te.ao.subset %>%
  simsum(estvarname = "theta", se = "se", methodvar = "method", true = 0.251314428, by = "dgm", ref = "onestep", x = TRUE) %>%
  summary() %>%
  tidy() %>%
  kable()

te.cl.subset %>%
  simsum(estvarname = "theta", se = "se", methodvar = "method", true = 0.133531393, by = "dgm", ref = "onestep", x = TRUE) %>%
  summary() %>%
  tidy() %>%
  kable()

#PLOTS

#including multiple estimands at once
#https://cran.r-project.org/web/packages/rsimsum/vignettes/E-custom-inputs.html

#total effects
te.eop.subset$par <- "1.eop"
te.ao.subset$par <- "2.ao"
te.cl.subset$par <- "3.cl"

te <- rbind(te.eop.subset,te.ao.subset,te.cl.subset)
is.data.frame(te) #[1] TRUE

te$true <- 0.356674944

dim<-dim(te)
size<-dim[1]

for(j in 1:size) {
if(te[j,7] == "2.ao")
  te[j,8] <- 0.251314428
if(te[j,7] == "3.cl")
  te[j,8] <- 0.133531393
}
  
ms1 <- multisimsum(
  data = te,
  par = "par", true = "true",
  estvarname = "theta", se = "se", methodvar = "method",
  by = "dgm", ref = "onestep", x = TRUE
)

sms1 <- summary(ms1)

#bias
ggplot(tidy(sms1, stats = "bias"), aes(x = factor(method, level = c('onestep', 'bch', 'modal', 'npcd', 'incpcd', 'upcd')), y = est, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_point() +
  geom_errorbar(width = 1 / 3) +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(par ~ dgm) +
  labs(x = "Method", y = "Bias")

#bias-eliminated coverage
ggplot(tidy(sms1, stats = "becover"), aes(x = factor(method, level = c('onestep', 'bch', 'modal', 'npcd', 'incpcd', 'upcd')), y = est, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0.95, color = "red", lty = "dashed") +
  geom_point() +
  geom_errorbar(width = 1 / 3) +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(par ~ dgm) +
  labs(x = "Method", y = "Bias eliminated coverage")

#indirect effects
tnie.eop.subset$par <- "1.eop"
tnie.ao.subset$par <- "2.ao"
tnie.cl.subset$par <- "3.cl"

tnie <- rbind(tnie.eop.subset,tnie.ao.subset,tnie.cl.subset)

tnie$true <- 0.05702914

for(j in 1:size) {
  if(tnie[j,7] == "2.ao")
    tnie[j,8] <- 0.080874497
  if(tnie[j,7] == "3.cl")
    tnie[j,8] <- 0.011555227
}

ms1 <- multisimsum(
  data = tnie,
  par = "par", true = "true",
  estvarname = "theta", se = "se", methodvar = "method",
  by = "dgm", ref = "onestep", x = TRUE
)

sms1 <- summary(ms1)

#bias
ggplot(tidy(sms1, stats = "bias"), aes(x = factor(method, level = c('onestep', 'bch', 'modal', 'npcd', 'incpcd', 'upcd')), y = est, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_point() +
  geom_errorbar(width = 1 / 3) +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(par ~ dgm) +
  labs(x = "Method", y = "Bias")

#direct effects
pnde.eop.subset$par <- "1.eop"
pnde.ao.subset$par <- "2.ao"
pnde.cl.subset$par <- "3.cl"

pnde <- rbind(pnde.eop.subset,pnde.ao.subset,pnde.cl.subset)

pnde$true <- 0.299645804

for(j in 1:size) {
  if(pnde[j,7] == "2.ao")
    pnde[j,8] <- 0.170439931
  if(pnde[j,7] == "3.cl")
    pnde[j,8] <- 0.121976166
}

ms1 <- multisimsum(
  data = pnde,
  par = "par", true = "true",
  estvarname = "theta", se = "se", methodvar = "method",
  by = "dgm", ref = "onestep", x = TRUE
)

sms1 <- summary(ms1)

#bias
ggplot(tidy(sms1, stats = "bias"), aes(x = factor(method, level = c('onestep', 'bch', 'modal', 'npcd', 'incpcd', 'upcd')), y = est, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_point() +
  geom_errorbar(width = 1 / 3) +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(par ~ dgm) +
  labs(x = "Method", y = "Bias")

#########################################
#WITHOUT APPLYING EXCLUSIONS
#########################################
#PLOTS

#total effects
te.eop$par <- "1.eop"
te.ao$par <- "2.ao"
te.cl$par <- "3.cl"

te <- rbind(te.eop,te.ao,te.cl)
is.data.frame(te) #[1] TRUE

te$true <- 0.356674944

dim<-dim(te)
size<-dim[1]

for(j in 1:size) {
  if(te[j,7] == "2.ao")
    te[j,8] <- 0.251314428
  if(te[j,7] == "3.cl")
    te[j,8] <- 0.133531393
}

ms1 <- multisimsum(
  data = te,
  par = "par", true = "true",
  estvarname = "theta", se = "se", methodvar = "method",
  by = "dgm", ref = "onestep", x = TRUE
)

sms1 <- summary(ms1)

#bias
ggplot(tidy(sms1, stats = "bias"), aes(x = factor(method, level = c('onestep', 'bch', 'modal', 'npcd', 'incpcd', 'upcd')), y = est, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_point() +
  geom_errorbar(width = 1 / 3) +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(par ~ dgm) +
  labs(x = "Method", y = "Bias")

#indirect effects
tnie.eop$par <- "1.eop"
tnie.ao$par <- "2.ao"
tnie.cl$par <- "3.cl"

tnie <- rbind(tnie.eop,tnie.ao,tnie.cl)

tnie$true <- 0.05702914

for(j in 1:size) {
  if(tnie[j,7] == "2.ao")
    tnie[j,8] <- 0.080874497
  if(tnie[j,7] == "3.cl")
    tnie[j,8] <- 0.011555227
}

ms1 <- multisimsum(
  data = tnie,
  par = "par", true = "true",
  estvarname = "theta", se = "se", methodvar = "method",
  by = "dgm", ref = "onestep", x = TRUE
)

sms1 <- summary(ms1)

#bias
ggplot(tidy(sms1, stats = "bias"), aes(x = factor(method, level = c('onestep', 'bch', 'modal', 'npcd', 'incpcd', 'upcd')), y = est, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_point() +
  geom_errorbar(width = 1 / 3) +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(par ~ dgm) +
  labs(x = "Method", y = "Bias")

#direct effects
pnde.eop$par <- "1.eop"
pnde.ao$par <- "2.ao"
pnde.cl$par <- "3.cl"

pnde <- rbind(pnde.eop,pnde.ao,pnde.cl)

pnde$true <- 0.299645804

for(j in 1:size) {
  if(pnde[j,7] == "2.ao")
    pnde[j,8] <- 0.170439931
  if(pnde[j,7] == "3.cl")
    pnde[j,8] <- 0.121976166
}

ms1 <- multisimsum(
  data = pnde,
  par = "par", true = "true",
  estvarname = "theta", se = "se", methodvar = "method",
  by = "dgm", ref = "onestep", x = TRUE
)

sms1 <- summary(ms1)

#bias
ggplot(tidy(sms1, stats = "bias"), aes(x = factor(method, level = c('onestep', 'bch', 'modal', 'npcd', 'incpcd', 'upcd')), y = est, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  geom_point() +
  geom_errorbar(width = 1 / 3) +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(par ~ dgm) +
  labs(x = "Method", y = "Bias")

################################################end of script#####################################################################
