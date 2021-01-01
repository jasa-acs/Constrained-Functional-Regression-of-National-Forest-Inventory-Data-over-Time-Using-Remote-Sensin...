############################################################################################
########################		Input     ##################################################
############################################################################################

## The output of TC_Model_BEST.R for all TCs where BEST should be replaced by the best TC model number (IV in our case) from the five candidate model numbers from Appendix A.1 in the supplementary materials - I, II, IIIA, IIIB and IV.

############################################################################################
########################		Output     ##################################################
############################################################################################

## The output is "best_TC_model_outputs.Rdata". This data contains a list reconstructed_TC consisting of four components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the reconstructed TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 

###################################################################################################
###################################################################################################
###################################################################################################

## TC data part

	BEST = readline(" Best model number from I, II, IIIA, IIIB and IV ")

	load(paste0("Inputs_and_Outputs/Model_", BEST , "_TC_1_outputs.Rdata"))
	
	locations = fill_up_missing_data$locations ## location of the cells

	S = nrow(fill_up_missing_data[[2]])

	T = ncol(fill_up_missing_data[[2]])

	TC_all = array (0, dim = c(S, T, 3))

	for(TC in 1:3)
		{

		load(paste0("Inputs_and_Outputs/Model_", BEST , "_TC_", TC , "_outputs.Rdata"))
	
		TC_all[ , , TC] = eval(parse(text = paste0("fill_up_missing_data$Transformed_TC", TC )))
	
		}

	TC1 = TC_all[ , , 1]

	TC2 = TC_all[ , , 2]

	TC3 = TC_all[ , , 3]

#########################################################################################################
#########################################################################################################

## saving locations and complete TC features

	reconstructed_TC = list(locations = locations, reconstructed_TC1 = TC1,  reconstructed_TC2 = TC2, reconstructed_TC3 = TC3)

	save(reconstructed_TC, file = paste0("Inputs_and_Outputs/best_TC_model_outputs.Rdata"))