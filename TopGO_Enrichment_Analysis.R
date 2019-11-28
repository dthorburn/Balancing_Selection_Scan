#########################################################################################################
#										TopGO Enrichment Analaysis											
#########################################################################################################

library("BiocManager")
library("GO.db")
library("topGO")
library("biomaRt")
library("dplyr")
library("AnnotationDbi")
library("org.Dr.eg.db")
library("data.table")

set.seed(101011)

######################## Creating the base universe
base_universe <- read.delim(file = "All_autosomal_gasAcu1_Universe.txt", header = TRUE, sep = "\t"); names(base_universe) <- "Gene_ID"; base_universe["Significant"] <- 0

## Reading in the candidate genes
gene_list <- read.csv(paste("Population_Gene_List_Simplified.csv", sep = ""))

## Applying the candidate gene list to the universe
for(gene in 1:nrow(gene_list)){
	row <- which(grepl(pattern = as.character(gene_list$Gene_ID[gene]), x = base_universe$Gene_ID))
	base_universe$Significant[row] <- 1
}

write.table(file = "./Universes/Population_Autosomal_Universe.txt", x = base_universe, col.names = FALSE, row.names = FALSE)

######################## Running topGO
## Reading in orthologs
orthos <- read.delim(file = "Orthologs_Ga_Dr.txt", header = TRUE, sep = ","); names(orthos)[1] <- "Gene_ID"

## Reading universe back in 
temp_data <- read.delim("./Universes/Population_Autosomal_Universe.txt", header = FALSE)
names(temp_data) <- c("Gene_ID", "Significant")

## Setting both as data tables. Then, sampling the zebrafish orthologs at random if more than 1 is present.
setDT(orthos); setDT(temp_data)
temp_data2 <- orthos[temp_data, on = .(Gene_ID),  {ri <- sample(.N, 1L)
					.(Zebrafish.gene.stable.ID = Zebrafish.gene.stable.ID[ri], Significant = Significant)}, by = .EACHI]
## Taking just the zebrafish genes
all_genes_orthos <- subset(temp_data2, grepl("ENSDAR", Zebrafish.gene.stable.ID))
## Counting the genes converted
cat(paste0("Number of G.aculeatus candidate genes: ", sum(temp_data$Significant), "\n"))
cat(paste0("Total number of G.aculeatus genes: ", nrow(temp_data), "\n"))
## Taking the dataframe I want for the analysis - the Zebrafish gene set. 
Zeb_GeneSet <- as.factor(all_genes_orthos$Significant)
## Assigning the gene name as the column header as in tutorial
names(Zeb_GeneSet) <- all_genes_orthos$Zebrafish.gene.stable.ID
cat(paste0("Number of candidate orthologs found: ", sum(grepl("1", Zeb_GeneSet)), "\n"))
cat(paste0("Total number of orthologs found: ", length(Zeb_GeneSet), "\n"))

## Running the GO analysis for biological function
go_data <- new("topGOdata",
   	ontology = "BP",
   	allGenes = Zeb_GeneSet,
	nodeSize = 5,
	annotationFun = annFUN.org,
	mapping = "org.Dr.eg",
	ID = "ensembl")

## Running the fishers exact test on the results of the GO analysis
go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")

go_table <- GenTable(go_data, weightFisher = go_test,
    orderBy = "weightFisher", ranksOf = "weightFisher",
	topNodes = sum(score(go_test) < 1.00))
go_table2 <- subset(go_table, Significant >= 3 & Annotated >=3)

## Creating a new column in the dataframe and 
go_table2["FDR_p_value"] <- p.adjust(as.numeric(gsub(go_table2$weightFisher, pattern = "<", replacement ="")), method = "fdr")

write.csv(x = go_table2, file = paste("./Results/Population_Zeb_GO_Results.csv", sep = ""), row.names = FALSE)
