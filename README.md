# latentclass-mediation
R and Mplus code to perform counterfactual mediation with a latent class exposure using updated PCD (and comparing against existing methods)

Code here is to perform counterfactual mediation with a latent class exposure comparing a new method (updated pseudo class draws; uPCD) to existing methods (onestep, bias-adjusted three step (BCH), modal class assignment, non-inclusive PCD, and inclusive PCD). There are three Mplus input files to generate the simulated data across three different levels of class separation (entropy). There are two R files containing code to compare the six methods (one using simulated data, and one using ALSPAC data as an applied example) each with corresponding Mplus input files needed to run the code. There is an additional R file which analyses the dataset of estimates from the simulated data. Finally, there is a webpage generated using R markdown (upcd-for-sim1-poor-entropy) describing the code to run uPCD on the first simulated dataset with poor entropy.

1a. Mplus input file to generate 500 simulated datasets with poor entropy

1b. Mplus input file to generate 500 simulated datasets with medium entropy

1c. Mplus input file to generate 500 simulated datasets with good entropy

2a. R file to analyse simulated datasets across all entropy levels

2b. Mplus input file for unconditional latent class model

2c. Mplus input file for mediation model using onestep

2d. Mplus input file for mediation model using bias-adjusted three step (BCH)

2e. Mplus input file for mediation model using modal class assignment

2f. Mplus input file for mediation model using non-inclusive PCD

2g. Mplus input file for inclusive latent class model

2h. Mplus input file for mediation model using inclusive PCD

2i. Mplus input file for mediation model using updated PCD

3. R file to analyse the dataset of the estimates from the simulated data

4a. R file to analyse ALSPAC dataset for the applied example

4b. Mplus input file for unconditional latent class model

4c. Mplus input file for mediation model using onestep

4d. Mplus input file for mediation model using bias-adjusted three step (BCH)

4e. Mplus input file for mediation model using modal class assignment

4f. Mplus input file for mediation model using non-inclusive PCD

4g. Mplus input file for inclusive latent class model

4h. Mplus input file for mediation model using inclusive PCD

4i. Mplus input file for mediation model using updated PCD
