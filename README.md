# Balancing Selection Scan
## Overview

This is the repository for a future publication looking into scanning for balancing selection in an organism that has been sampled in a manner that allows questions to be asked about parapatric populations over a continuum of divergence. 

## NCD Scan

Most of the details for the NCD scan 

## Neutral Simulations

We used the software [SLiM v3.3](https://messerlab.org/slim/) to conduct our neutral simulations. Here, we aim to simulate neutral evolution in a population pair with migration. Mutation rate was taken from [Feulner et al. 2015](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004966) and recombination rate from [Roesti et al. 2013](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12322), and migration rate between populations was calculated using the software [MIGRATE](https://peterbeerli.com/migrate-html5/index.html). For more details, see the annotated SLiM configuration file ```SLiM_Config.slim```. 

To then perform the scans using the same metrics in the above genome scans, the R script ```R_Script``` parses the SLiM output table, and additionally taken input of a sequence of the same length as the simulated chromosome, and applies the mutations. The model of nucleotide substitution was JC69, but you can input any matrix by changing the values in the substitution ```if``` statements. 
