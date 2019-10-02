#########################################################################################################
#						              GO Analysis for Paralell Versus Non-Parallel Evolution											
######################################################################################################### 
## Date: 26/09/2019
## Involved: Miles
## Task: Create a plot from the LD analysis in PLINK

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

## Creating the vector of bins - here it's 500bp to 30Mb
bin_vector <- seq(from = 0, to =  30000000, by = 500)

## Sending the path to the *.ld file to the R environment
args <- commandArgs(trailingOnly = TRUE)
LD_File <- args[1]

## Reading in the data - use data.table's fread for large files, it's amazingly efficient.
data <- fread(LD_File, header = TRUE, sep = " ")

## Calculating the physical distance between pairwise sites
data[,"Raw_Difference":=data$BP_A - data$BP_B]
data[,"Difference":=abs(data$Raw_Difference)]

## Creating a column that tells me which 500bp bin they are in - the name was from my first attempt when it was 10Kb windows
data[,"bin_10k":=findInterval(data$Difference, bin_vector)]

## Function to calculate the summary statistics for each bin 
## Putting the name of the raw data object in the function is not ideal, but it'll do for now. 
bin_stats <- function(bin_xk){
	## subsetting the data
	temp_data <- subset(data, bin_10k == bin_xk)
	## Creating an output table with everything I could want
	output_df <- data.table(Difference = mean(temp_data$Difference), 
							Bin = bin_xk, 
							N_Sites = temp_data$BP_A %>% unique %>% length,
							Q1_R2 = summary(temp_data$R2)[2] %>% as.numeric, 
							Q3_R2 = summary(temp_data$R2)[5] %>% as.numeric, 
							Median_R2 = summary(temp_data$R2)[3] %>% as.numeric,
							Mean_R2 = mean(temp_data$R2))
	return(output_df)
}

## Creates a data.table object with the results. 
to_output <- lapply(unique(data$bin_10k), bin_stats) %>% rbindlist(.)
