############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_FIA_data_analysis.R that contains all necessary functions.

## The output of TC_best_model_output.R that is "best_TC_model_outputs.Rdata"

## best_TC_model_outputs.Rdata contains a list reconstructed_TC consisting of four components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the complete TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 

## hypothetical_basal_area_data.csv.
## hypothetical_basal_area_data.csv that is a matrix of five columns where first two entries in each row represents the coordinates of the field plot of live tree basal area data, next two entries represents the time point with year and month number (a number from 1 to 12) of the correspondig year, respectively and the last entry shows the corresponding the live tree basal area measurement.

############################################################################################
########################		Output     #################################################
############################################################################################

## The output is "best_TC_model_outputs_and_basal_area_data.Rdata". This data contains a list reconstructed_TC_BALIVE consisting of five components:

## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the reconstructed TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where rows and colums correspond to pixels and time points (months), respectively. 
## (iii) The fifth component is the live tree basal area in the form of a 10000x120 matrix where rows and colums correspond to pixels and time points (months), respectively. The coccurance of -99 within this matrix indicates unobserved data at those space-time combinations.

###################################################################################################
###################################################################################################
###################################################################################################

## TC data part

	source("FIA_Codes/Functions_for_FIA_data_analysis.R")

	load(paste0("Inputs_and_Outputs/best_TC_model_outputs.Rdata"))	## loading the data 

	locations = reconstructed_TC$locations
	
	S = nrow(locations)

	T = ncol(reconstructed_TC$Transformed_TC1)

#########################################################################################################
#########################################################################################################

## Basal area data part

## If the user wants to use own true basal area data, then
## (i) follow the instruction about the structure of the basal area data mentioned in 1(b) in the README.pdf file and 
## (ii) replace the file name as basal_area_data.csv in the following line.

	Basal_area_data = read.csv("Inputs_and_Outputs/hypothetical_basal_area_data.csv")

	BALIVE = matrix(-99, S, T)

	x_merge = 16		## The number of adjacent cells in TC features to be aggregated in X coordinate
	
	y_merge = 16		## The number of adjacent cells in TC features to be aggregated in Y coordinate
	
	res_loc = Basal_area_data[ ,c(1,2)]
 
    res_data = Basal_area_data[ ,5]
  
    sel_res_loc = which((res_loc[,1] >= min(locations[,1]) - 30*(x_merge/2)) & (res_loc[,1] <= max(locations[,1]) + 30*(x_merge/2)) & (res_loc[,2] >= min(locations[,2]) - 30*(y_merge/2)) & (res_loc[,2] <= max(locations[,2]) + 30*(y_merge/2)))

    sel_coarse_loc = c(apply(res_loc[sel_res_loc,], 1, position_covariate, locations))
  
	time = (Basal_area_data[,3]-2003)*12 + Basal_area_data[,4]
    
	loc_time = cbind(sel_coarse_loc, time[sel_res_loc])

	for (i in 1:dim(loc_time)[1]) BALIVE[sel_coarse_loc[i], time[sel_res_loc[i]]] = res_data[sel_res_loc[i]]

#########################################################################################################
#########################################################################################################

## saving locations, complete TC features and live tree basal area 

	reconstructed_TC_BALIVE = reconstructed_TC
	
	reconstructed_TC_BALIVE$BALIVE = BALIVE

	save(reconstructed_TC_BALIVE, file = paste0("Inputs_and_Outputs/best_TC_model_outputs_and_basal_area_data.Rdata"))