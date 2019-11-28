#########################################################################################################
#													Plotting MSMC Output														
#########################################################################################################
library(ggplot2)
library(data.table)
library(dplyr)


MSMC <- fread("All_MSMC_Results_31102019.csv")
MSMC[,"Ecotype" := gsub("^.._", "", Population)]

## Specify mutation rate
mu <- 1.4812e-8
## Specify generation time in years
gen <- 1

Colours <- c("#990099", "#990099", "#0000FF", "#0000FF", "#0099FF", "#0099FF", "#33FF66", "#33FF66", "#CC0000", "#CC0000")
ggplot(MSMC, aes(x = left_time_boundary/mu*gen, y = ((1/lambda_00)/(2*mu))/1000, linetype = Ecotype, colour = Population)) +
	geom_step() +
	ylim(0,100000)+
	theme_classic(base_size = 20) +
	scale_linetype_manual(values = c("solid", "longdash")) +
	scale_colour_manual(limits = c("CA_R", "CA_L", "G1_R", "G1_L", "G2_R", "G2_L", "NO_R", "NO_L", "US_R", "US_L"), values = Colours) +
	scale_x_continuous(limits = c(0,10000), expand = c(0,0), breaks = seq(0,10000,2000), labels = seq(0,10,2)) +
	scale_y_continuous(limits = c(0,80), expand = c(0,0)) +
	labs(x = expression(paste("Years ago (×10"^3*")")), y = expression(paste("Effective Population Size (×10"^3*")"))) +
	theme(legend.position = "none")
	ggsave(file = "./MSMC_Populations_NoLegend.png", width = 400, height = 200, units = "mm")
