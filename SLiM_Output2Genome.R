#########################################################################################################
#									SLiM Output Manipulation		 									
#########################################################################################################
library(stringr)
library(dplyr)
library(ape)
library(seqinr)
library(PopGenome)
library(ggplot2)
library(vcfR)
library(data.table)
library(lme4)
library(lmerTest)
library(ggsignif)
library(multcomp)

######################################################################## Task 1: Apply simulated mutations to nucleotide sequences
## Read in simulation data
#sim_data <- read.delim("SLiM_1789597437506_simplified_migration2.txt", skip = 6, header = FALSE, sep = " ")
sim_data <- read.csv("SLiM_1811605093866_simplified_migration_24M_Muts.csv", header = FALSE)

## Removing the extra data at the end of the text file, and naming columns
sim_data_muts <- subset(sim_data, V7 == "p0" | V7 == "p1")
names(sim_data_muts) <- c("Unique_ID1", "Unique_ID2", "Mutation_Type", "Position", "Selection_Coef", "Dominance_Coef", "Subpopulation",
							"Generation", "Prevelance")
#sim_data_muts$Unique_ID1  <- sim_data_muts$Unique_ID1 %>% as.numeric
## Reading in the genome data section. Here, I combined homologous chromosomes within each individual to make a single consensus sequence at the end of this.
sim_data_genomes <- read.csv("SLiM_1811605093866_simplified_migration_24M_Genomes.csv", header = TRUE)

## Read in what I am using as the ancestral sequence
ancestral_seq <- read.fasta("Gasterosteus_aculeatus_gasAcu1_chr1.fasta")
#seqs <- read.FASTA("Gasterosteus_aculeatus_gasAcu1_chr1.fasta", type = "DNA")
## Assigning each of the names in the simulated consensus genomes the ancestral sequence at the same length as the simulated contig.  
#for(name in names(sim_data_genomes)){
	print(name)
	sequence <- ancestral_seq$chrI[1:23999976]
	output <- paste(">", name, "\n", toupper(paste(sequence, collapse = "")), "\n", sep = "")
	cat(output, file = "Ancestral_sequence_SLiM.fasta", append = TRUE)1
}
ancestral_sequences <- read.fasta("Ancestral_sequence_SLiM.fasta")
ancestral_sequences2 <- read.FASTA("Ancestral_sequence_SLiM.fasta", type = "DNA")

str(ancestral_sequences2)
## Calculating the nucleotide substitution matrix when a genome ID and ancestral nucleotide is input - here using the F81 matrix.
create_matrix <- function(temporary_sequences){
	#temp_seqs <- as.DNAbin(temporary_sequences)
	temp_freqs <- temporary_sequences %>% base.freq
	temp_pi_tag <- c(temp_freqs[4], temp_freqs[2], temp_freqs[1], temp_freqs[3])
	temp_sub_mat <- sub_F81(temp_pi_tag)
	rownames(temp_sub_mat$Q) <- c("T", "C", "A", "G"); colnames(temp_sub_mat$Q) <- c("T", "C", "A", "G")
	return(temp_sub_mat$Q)
}

apply_mutation <- function(dist_matrix, ancestral_state){
	if(ancestral_state == "A"){
		sample(x = c("T", "C", "G"), size = 1, replace = FALSE, prob = as.numeric(dist_matrix["A",c(1,2,4)]))
	} else if(ancestral_state == "T"){
		sample(x = c("C", "A", "G"), size = 1, replace = FALSE, prob = as.numeric(dist_matrix["T",2:4]))
	} else if(ancestral_state == "G"){
		sample(x = c("T", "C", "A"), size = 1, replace = FALSE, prob = as.numeric(dist_matrix["G",1:3]))
	} else if(ancestral_state == "C"){
		sample(x = c("T", "A", "G"), size = 1, replace = FALSE, prob = as.numeric(dist_matrix["C",c(1,3,4)]))
	}	 
}
counter = 0
for(mutation in 1:nrow(sim_data_muts)){
	##Just keeping this line in case I need to debug something later
	#print(paste("Processing mutation in generation", sim_data_muts[mutation,"Generation"]))
	temp_mutation <- sim_data_muts[mutation,]
	## Returning which of the genome columns have the mutation ID value. 
	unique(colnames(sim_data_genomes))[which(sim_data_genomes == temp_mutation$Unique_ID1, arr.ind = TRUE)[,2]] -> ind_names
	## Ensuring the prevelance and the genomes columns adds up to the same value
	if(length(ind_names) == temp_mutation$Prevelance){
		## Because I couldn't figure out another way, I am using the original ancestral sequence to find the ancestral state, on the assumption that no site will have a mutation twice under neutral evolution 
		## and over such a low number of generations. Not great, but the best I can do for now.
		anc_state <- ancestral_seq$chrI[temp_mutation$Position] %>% toupper
		## Getting around N's being in the original contig by randomly choosing a nucleotide. 
		if(anc_state == "N"){
			anc_state <- sample(x = c("A", "C", "T", "G"), replace = FALSE, size = 1)
			print(paste("Altering position", temp_mutation$Position, "from N to", anc_state))
		}
		## If there were more than ~4000 mutation, I would update the base frequencies every generation, but given the low numbers, the changes to a ~28Mb contig are negligible.
		ancestral_sequences2[ind_names] %>% create_matrix %>% apply_mutation(., anc_state) -> new_nucleotide
		## Applying the mutation - ugly, but it'll do.
		if(sum(grepl("p0.0",  ind_names)) == 1){ancestral_sequences$p0.0[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p0.2",  ind_names)) == 1){ancestral_sequences$p0.2[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p0.4",  ind_names)) == 1){ancestral_sequences$p0.4[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p0.6",  ind_names)) == 1){ancestral_sequences$p0.6[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p0.8",  ind_names)) == 1){ancestral_sequences$p0.8[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p0.10", ind_names)) == 1){ancestral_sequences$p0.10[temp_mutation$Position] <- new_nucleotide}
		if(sum(grepl("p1.0",  ind_names)) == 1){ancestral_sequences$p1.0[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p1.2",  ind_names)) == 1){ancestral_sequences$p1.2[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p1.4",  ind_names)) == 1){ancestral_sequences$p1.4[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p1.6",  ind_names)) == 1){ancestral_sequences$p1.6[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p1.8",  ind_names)) == 1){ancestral_sequences$p1.8[temp_mutation$Position]  <- new_nucleotide}
		if(sum(grepl("p1.10", ind_names)) == 1){ancestral_sequences$p1.10[temp_mutation$Position] <- new_nucleotide}
		print(paste("Mutation in position", temp_mutation$Position, "of", paste(ind_names, collapse = ", "),"from", anc_state, "to", new_nucleotide))
	}	else	{
		counter <- counter + 1
		print("ERROR: Incorrect number of individuals in this mutation:")
		print(temp_mutation)
	}
}
counter 

## Note, this was performed for each population seperately, so that I can scan the population individually too. 
write.fasta(file = "SLiM_Output_Pair.fasta", ancestral_sequences[1:12], names = names(ancestral_sequences)[1:12])

########################################### Task 2: The genome scan
pop1_names <- c("p0.0", "p0.2", "p0.4", "p0.6", "p0.8", "p0.10")

## A directory with the fasta created above for just the P1 population. 
neutral_genomes <- readData("./Alignment")

## Only needed for the population pair section
pop.calls2 <- readVCF("./Alignment", frompos = 1, topos = 23999976, numcols = 3000000, samplenames = pop1_names)

pop.calls2 <- set.populations(neutral_genomes, list(pop1_names, pop2_names))
sliding_window <- sliding.window.transform(pop.calls2, width = 7500, jump = 1500, type = 2)

sliding_window <- sliding.window.transform(neutral_genomes, width = 7500, jump = 1500, type = 2)
test1 <- neutrality.stats(sliding_window)
test2 <- diversity.stats(sliding_window, pi = TRUE, keep.site.info = FALSE)
results1 <- as.data.frame(get.neutrality(test1, theta = TRUE)[1]);	results2 <- tibble::rownames_to_column(results1, var = "Window_Start_Finish")
results3 <- as.data.frame(get.diversity(test2)[[1]]);	results4 <- tibble::rownames_to_column(results3, var = "Window_Start_Finish")

all.results <- data.frame( 	Chr = "Sim_Chr1",
									Tajima_D = as.numeric(results2$Tajima.D),
									N.Sites = as.numeric(results2$n.segregating.sites),
									Watterson = as.numeric(results2$theta_Watterson)/7500,
									Achaz_Watterson = as.numeric(results2$theta_Achaz.Watterson),
									FuLi_D = as.numeric(results2$Fu.Li.D),
									FuLi_F = as.numeric(results2$Fu.Li.F),
									Nuc_Diversity = as.numeric(results4$nuc.diversity.within),
									Pi = as.numeric(results4$Pi)/7500,
									Window_Start_Finish = results2$Window_Start_Finish,
									Window_Width = rep(7500, nrow(results2)),
									Step_Size = rep(1500, nrow(results2)),
									Start = gsub(pattern = " .*", "", results2$Window_Start_Finish), 
									End = gsub(pattern = ".*- ", "", results2$Window_Start_Finish) %>% gsub(pattern = " :", replacement = "", .))


write.csv(file = "P1_both_Neutral_Scan_Results.csv", all.results)



############################################# Part 3: NCD
pop1_names <- c("p0.0", "p0.2", "p0.4", "p0.6", "p0.8", "p0.10")
pop2_names <- c("p1.0", "p1.2", "p1.4", "p1.6", "p1.8", "p1.10")
## The vcf was created by converting the fasta with the software PGDSpider. 
temp_chrom <- read.vcfR("SLiM_Output_Pair.vcf")
## The GT_separator is what separates gt info, I think this changes to | in phased VCFs. 
Freq_Fun <- function(input, value){
	str_count(string = input, pattern = value)/str_count(string = input, pattern = "[0-9]")
}
NCD_Gen <- function(vcf_obj, population_ID){
	output <- data.table(CHR = vcf_obj@fix[,"CHROM"] %>% gsub("Loci_", "", .) %>% as.integer(), 
							POS = vcf_obj@fix[,"POS"] %>% as.integer(), 
							ID = paste(vcf_obj@fix[,"CHROM"], vcf_obj@fix[,"POS"], sep = "|"),
							REF = vcf_obj@fix[,"REF"], ALT = vcf_obj@fix[,"ALT"],
							AF1 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("0") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric(),
							AF2 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("1") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric(),
							AF3 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("2") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric(),
							AF4 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("3") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric(),
							AF5 = gsub(pattern = ":.*", replacement =  "", vcf_obj@gt[,as.character(Population_list[,population_ID])]) %>% as.data.frame() %>%
									apply(MARGIN = 1, FUN = paste, collapse = "|") %>% Freq_Fun("[4-9]") %>% gsub(pattern = "^0$", replacement = NA) %>% as.numeric())
	## Finally, once the AFs are calculated, this creates the MAF column by taking the min values of the AF columns. 
	output %>% mutate(MAF = pmin(AF1, AF2, AF3, AF4, AF5, na.rm = TRUE)) %>% as.data.table()
}
Population_list <- data.frame(  P0=c("p0.0", "p0.2", "p0.4", "p0.6", "p0.8", "p0.10"), 
								P1=c("p1.0", "p1.2", "p1.4", "p1.6", "p1.8", "p1.10"))
P0_NCD <- NCD_Gen(temp_chrom, "P0")
P1_NCD <- NCD_Gen(temp_chrom, "P1")

P0_NCD1 <- NCD1(P0_NCD)
P1_NCD1 <- NCD1(P1_NCD)

IS_Filt_P0 <- subset(P0_NCD1, N_SNPs_cor >= 10)
IS_Filt_P1 <- subset(P1_NCD1, N_SNPs_cor >= 10)

write.csv(file = "P0_Simulated_NCD_Scan_Results.csv", IS_Filt_P0)
write.csv(file = "P1_Simulated_NCD_Scan_Results.csv", IS_Filt_P1)


