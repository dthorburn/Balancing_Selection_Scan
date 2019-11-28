#########################################################################################################
#										Example PopGenome Analysis	 									
#########################################################################################################
library(PopGenome)
library(dplyr)
library(ggplot2)
library(tibble)

## Parameters for the scan
window_width = 7500
step_size = 1500

## This was all part of a series of nested loops, hence the paste function. Note: PopGenome can read in .vcf.gz files
input <- paste("~/Chapter1/VCFs/Attempt3/Population/Population/All_population_chr1_final.vcf.gz", sep = "")

## Reading in the vcf.gz file. "topos" is the end coordinate of interest. Here, I used the length of the chromosomes. "tid" is the chromosome name, chrVII for you. 
## "numcols" is just the amount of data to read in, I used a random high value, more than the number of SNPs in my VCFs. 
## "samplenames" is where you can identify which of the samples to read in if you just want one population, for example. 
## "chrom_lengths" is a matrix with the length of each chromosome in the genome. 
## "pop_names" is a list with each item being a vector of names that corresponds to the columns heads in the vcf file
pop.calls <- readVCF(input1, frompos = 1, topos = chrom_lengths[chrom,2], tid = r.chrom, numcols = 5000000, samplenames = pop_names[[pops]])

## Transforming the chromsome into a set of overlapping windows
sliding_window <- sliding.window.transform(pop.calls, width = window_width, jump = step_size, type = 2)

## Running Neutrality and Diversity tests on windows; Tajima, FuLi D+F, Watterson 
test1 <- neutrality.stats(sliding_window)
test2 <- diversity.stats(sliding_window, pi = TRUE, keep.site.info = FALSE)

## The following lines are my brute force way of getting what I wanted from the package. There are probably more elegant
## ways of doing this, but I wasn't so good at coding then. 
results1 <- as.data.frame(get.neutrality(test1, theta = TRUE)[1])
results2 <- tibble::rownames_to_column(results1, var = "Window_Start_Finish")
results3 <- as.data.frame(get.diversity(test2)[[1]])
results4 <- tibble::rownames_to_column(results3, var = "Window_Start_Finish")

## Creating a new dataframe with all the metrics I want to keep - again, pretty inelegant, but it did the trick. 
all.results <- data.frame( 	Chr = rep(as.factor(chrom),nrow(results2)),
							Tajima_D = as.numeric(results2$Tajima.D),
							N.Sites = as.numeric(results2$n.segregating.sites),
							Watterson = as.numeric(results2$theta_Watterson),
							Achaz_Watterson = as.numeric(results2$theta_Achaz.Watterson),
							FuLi_D = as.numeric(results2$Fu.Li.D),
							FuLi_F = as.numeric(results2$Fu.Li.F),
							Nuc_Diversity = as.numeric(results4$nuc.diversity.within),
							Pi = as.numeric(results4$Pi),
							Window_Start_Finish = results2$Window_Start_Finish)

wrte.csv(file = "CA_L_chr1_PopGenome_Scan_Results.csv", all.results)
