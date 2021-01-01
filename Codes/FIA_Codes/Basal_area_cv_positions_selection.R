############################################################################################
########################		Input     ##################################################
############################################################################################

## The output of TC_best_output_merge_basal_area.R that is "best_TC_model_outputs_and_basal_area_data.Rdata".

## best_TC_model_outputs_and_basal_area_data.Rdata contains a list reconstructed_TC_BALIVE consisting of five components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the complete TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 
## (iii) The fifth component is the live tree basal area in the form of a 10000x120 matrix where rows and colums correspond to pixels and time points (months), respectively. The coccurance of -99 within this matrix indicates unobserved data at those space-time combinations.

############################################################################################
########################		Output     #################################################
############################################################################################

## random_positions_for_CV_for_y1_model.Rdata is a list test_pos_list of 36 components where each component contains the postions of both zero, non-zero live tree basal area data those will be considered as test data in the correspondig holdout cross-validation method.

## random_positions_for_CV_for_y2_model.Rdata is a list test_pos_list of 36 components where each component contains only the postions of non-zero live tree basal area data those will be considered as test data in the correspondig holdout cross-validation method.

###############################################################################################################
################	selecting positions for test data for y^1 model		#######################################
###############################################################################################################

	set.seed(77777)

	load("Inputs_and_Outputs/best_TC_model_outputs_and_basal_area_data.Rdata")	## loading the data 

	#########################################################################################################
	
	locations_all = reconstructed_TC_BALIVE$locations
	
	S = nrow(locations_all) 

	T = 10
	
	count_full = matrix(-99,  S,  T)

	pos_0 = which(reconstructed_TC_BALIVE$BALIVE == 0, arr.ind = TRUE)

	pos_0 = pos_0[which(pos_0[ , 2] >  0 & pos_0[ , 2] < 121 ), ]

	pos_1 = which(reconstructed_TC_BALIVE$BALIVE > 0, arr.ind = TRUE)

	pos_1 = pos_1[which(pos_1[ , 2] >  0 & pos_1[ , 2] < 121 ), ]

	for(s in 1:(dim(pos_1)[1])) count_full[pos_1[s, 1], ceiling(pos_1[s, 2]/12)] = 1

	for(s in 1:(dim(pos_0)[1])) count_full[pos_0[s, 1], ceiling(pos_0[s, 2]/12)] = 0
	
	pos_avail =  which(count_full != -99)

	######################################################################################################
	
	n_CV = 36	## number of cross-validations

	n_total = length(pos_avail)

	l_test = 15	## length of test data set
	
	pos_1 = sample(1:n_total, floor(n_total/l_test)*l_test, replace = FALSE )

	pos_2 = setdiff(1:n_total, pos_1)

	pos_3 = sample(1:n_total, (n_CV*l_test- n_total), replace = FALSE )

	sel_pos = c(pos_1, pos_2, pos_3)
	
	test_pos_list = split(sel_pos, rep(1:n_CV, each = l_test))
	
	save(test_pos_list, file = "Inputs_and_Outputs/random_positions_for_CV_for_y1_model.Rdata")


###############################################################################################################
################	selecting positions for test data for y^2 model		#######################################
###############################################################################################################

	set.seed(7777)
	
	y_full = matrix(-99,  S,  T)

	pos_1 = which(reconstructed_TC_BALIVE$BALIVE > 0, arr.ind = TRUE)

	pos_1 = pos_1[which(pos_1[ , 2] >  0 & pos_1[ , 2] < 121 ), ]

	for(s in 1:(dim(pos_1)[1])) y_full [pos_1[s,1] , ceiling(pos_1[s,2]/12)] = log(reconstructed_TC_BALIVE$BALIVE[pos_1])[s]

	pos_avail =  which(y_full != -99)

	######################################################################################################
	
	n_CV = 36	## number of cross-validations

	n_total = length(pos_avail)

	l_test = 10	## length of test data set
	
	pos_1 = sample(1:n_total, floor(n_total/l_test)*l_test, replace = FALSE )

	pos_2 = setdiff(1:n_total, pos_1)

	pos_3 = sample(1:n_total, (n_CV*l_test- n_total), replace = FALSE )

	sel_pos = c(pos_1, pos_2, pos_3)
	
	test_pos_list = split(sel_pos, rep(1:n_CV, each = l_test))
	
	save(test_pos_list, file = "Inputs_and_Outputs/random_positions_for_CV_for_y2_model.Rdata")