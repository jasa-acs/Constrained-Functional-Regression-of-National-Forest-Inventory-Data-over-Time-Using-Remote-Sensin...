############################################################################################
############################################################################################

## selecting the positions of the test data for cross-validation

############################################################################################
########################		Input     ##################################################
############################################################################################

## the output of Aggregation.R that is "coarse_selected_data_with_missing_values.Rdata"
## "coarse_selected_data_with_missing_values.Rdata" contains a list coarse_creation that consists of two components:

## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.

## (ii) The second component contains the spatiotemporal data in form of a 10000x120x3 array where first and second dimensions correspond to pixels and time points (months), respectively. Along the third dimension, it contains values for TC1 (Brightness), TC2 (Greenness) and TC3 (Wetness). The occurrences of Zero (0) within this array indicate missing data at those space-time combinations.

############################################################################################
########################		Output     #################################################
############################################################################################

## the output is "random_positions_for_CV.Rdata" that contains a list test_pos_list_all consisting of 36 components where each component contains the positions of test data for the corresponding holdout method

############################################################################################
############################################################################################

set.seed(7777)

l_test = 25000		## number of samples that are considered as tested data

n_CV = 36	## Total number of times the holdout method is replicated 

load("Inputs_and_Outputs/coarse_selected_data_with_missing_values.Rdata") ## loading output of Aggregation.R that contains a list coarse_creation

data_TC1 = coarse_creation$coarse_data[ , , 1]

pos_0_all = which(data_TC1 == 0)	## vector of all missing positions

n_total = length(data_TC1)

non_missing_pos = setdiff(1:n_total, pos_0_all)

sel_pos = sample(non_missing_pos, l_test*n_CV, replace = FALSE)

test_pos_list_all = split(sel_pos, rep(1:n_CV, each = l_test))

## saving the list of the positions those will be considered as test data in cross-validation

save(test_pos_list_all, file = "Inputs_and_Outputs/random_positions_for_CV.Rdata")

