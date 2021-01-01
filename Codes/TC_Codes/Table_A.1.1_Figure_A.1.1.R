############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_TC_analysis.R that contains all necessary functions.

## The outputs of TC_Model_NUMBER.R and TC_Model_NUMBER_cv.R for all TC values where NUMBER should be replaced by any one of the five candidate model numbers from Appendix A.1.1 in the supplementary materials - I, II, IIIA, IIIB and IV

############################################################################################
########################		Output     ##################################################
############################################################################################

## The outputs are Table A.1.1 and Figure A.1.1 of Appendix A.1.3 in the supplementary materials

###################################################################################################
###################################################################################################
###################################################################################################

library(ggplot2)
library(sp)
library(grid)
library(gridBase)
library(lattice)
library(proj4)
library(RColorBrewer) # for color

source("TC_Codes/Functions_for_TC_analysis.R")

############################################################################################
########################		Table A.1.1     ############################################
############################################################################################

model_name = c("I", "II", "IIIA", "IIIB", "IV")

LL_BPIC = matrix(0, 6, 5)

CV =  matrix(0, 9, 5)

moran_I_all = list()

for(model_no in 1:5)
{

model = model_name[model_no]

for(TC in 1:3)
	{

	load(paste0("Inputs_and_Outputs/Model_", model, "_TC_", TC , "_outputs.Rdata"))
	
	LL_BPIC[c(TC, 3 + TC), model_no] = fill_up_missing_data$LL_BPIC
	
	moran_I_all[[(model_no - 1)*3 + TC ]] = fill_up_missing_data$moran_I[ , 2]

	CV_output = read.csv(paste0("Inputs_and_Outputs/Model_", model, "_TC_", TC, "_CV_outputs.csv"), header = T)
	
	CV[ c(TC, 3 + TC, 6 + TC), model_no] =  CV_output[1:3,1]

	}

}

a = rep(c("LL (x 10^5)", "BPIC (x 10^5)", "Absolute Bias", "Uncertainty", "Empirical Coverage"),  each = 3)

b = rep( c("TC1", "TC2", "TC3"), 5)

LL_BPIC = round(LL_BPIC/10^5, 3)

outputs = rbind(LL_BPIC, CV)

Table_A.1.1 = data.frame(cbind(a, b, outputs))

colnames(Table_A.1.1) = c("Measure", "Feature", "Model I", "Model II", "Model III-A", "Model III_B", "Model IV")

write.csv(Table_A.1.1, file = "Inputs_and_Outputs/Table_A.1.1.csv", row.names = FALSE) 

############################################################################################
########################		Figure A.1.1     ############################################
############################################################################################

TC_temp = rep( c("TC1", "TC2", "TC3"), each = length(moran_I_all[[1]]))

TC = rep(TC_temp, 5)

Model = rep(c("Model I", "Model II", "Model III-A", "Model III-B", "Model IV"), each = 3*length(moran_I_all[[1]]))

Moran_I = c( do.call(cbind, moran_I_all) )

moran_I_full = data.frame(list(Model = Model, TC = TC, Moran_I = Moran_I))

ggsave("Inputs_and_Outputs/Figure_A.1.1.pdf", height = 6, width = 10)

m <- ggplot(moran_I_full, aes(x = Moran_I) ) 
m + geom_histogram(bins = 80) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
facet_grid(TC ~ Model, scales = "free") +
labs(x = "Moran's I", y = "")+
expand_limits(x=c(-.3,1))+
scale_x_continuous( breaks = seq(-.3, 1, by=.3))

dev.off()