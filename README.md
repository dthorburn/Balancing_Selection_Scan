# Balancing Selection Scan
## Overview

This is the repository for a future publication looking into scanning for balancing selection in threespine stickelbacks that have been sampled in a structured manner that allows questions to be asked about parapatric populations over a continuum of divergence. 

## NCD Scan

Most of the details for the NCD scan [repository](https://github.com/dthorburn/NCD_Genome_Scan) I have. 

## Neutral Simulations

We used the software [SLiM v3.3](https://messerlab.org/slim/) to conduct our neutral simulations. Here, we aim to simulate neutral evolution in a population pair with migration. Mutation rate was taken from [Feulner et al. 2015](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004966) and recombination rate from [Roesti et al. 2013](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12322), and migration rate between populations was calculated from the genome-wide FST values. For more details, see the annotated SLiM configuration file ```SLiM_Config.slim```. 

To then perform the scans using the same metrics in the above genome scans, the R script ```R_Script``` parses the SLiM output table, and additionally taken input of a sequence of the same length as the simulated chromosome, and applies the mutations. The model of nucleotide substitution was JC69, but you can input any matrix by changing the values in the substitution ```if``` statements. 

## Linkage Disequilibrium Analysis

For the LD analysis, we used [Plink v1.9](https://www.cog-genomics.org/plink2/). As anyone who tries to do a genome-wide LD scan knows, you need to balance the number of pariwise calculations, to create a file size that is useable. My original attempts created single population and chromosome LD files that were several terabytes. Hence, for computational limitations, we limited the number of pariwise calculations to be calculated only if the two sites were within 1,000 vairants of one another. There was a 3Mb threshold too, but it was never met. 

To parse the .ld files, I wrote a function that will calculate summary statistics for each bin based on pairwise distance between variants. I used 500bp, but it's easily changed in the ```Calculate_Plink_LD_Decary.R``` script.

## GO Enrichment Analysis

The enrichment analysis was conducted in R, using the [TopGO](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf) package from the bioconductor project. The ```TopGO_Enrichment_Analysis.R``` script demonstrates how the autosomal universe was created, and how the orthologs were converted. Zebrafish orthologs were chosen at random if it there were multiple orthologous genes available. 

## Demography - MSMC

Demography was estiamted using the software [MSMC](https://www.nature.com/articles/ng.3015). For this, genomes were phased using [ShapeIt](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html), following the instruction manual for the read-aware approach. A second variant call was performed since information on all sites is required for this phasing approach. If requested, I will upload each step, but since I followed available tutorials online, I am not adding any more clarity. Here, the ```MSMC_GGplot.R``` file simply shows how to plot the output. 

## Diversity Statistics

The genome-wide diversity statistics were calculated using the [PopGenome](https://cran.r-project.org/web/packages/PopGenome/vignettes/An_introduction_to_the_PopGenome_package.pdf) package from the bioconductor project. The ```PopGenome_Scan.R``` script demonstrates how to read in a vcf, convert the scan to scanning windows, and perform the genome scan.
