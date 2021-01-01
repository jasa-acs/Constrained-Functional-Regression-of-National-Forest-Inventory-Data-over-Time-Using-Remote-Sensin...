############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_TC_analysis.R that contains all necessary functions.

## The data file "TC_features_with_missing_values.Rdata" that contains a list pixel_wise_TC_data consisting of two components:

## (i) The first component is a 2560000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.

## (ii) The second component contains the spatiotemporal data in form of a 2560000x120x3 array where first and second dimensions correspond to pixels and time points (months), respectively. Along the third dimension, it contains values for TC1 (Brightness), TC2 (Greenness) and TC3 (Wetness). The occurrences of Zero (0) within this array indicate missing data at those space-time combinations.


############################################################################################
########################		Output     ##################################################
############################################################################################

## The output is "coarse_selected_data_with_missing_values.Rdata" that contains a list coarse_creation consisting of two components:

## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers after aggregating adjacent 16x16 pixels from the study area described in Section 2 of the manuscript.

## (ii) The second component contains the spatiotemporal data in form of a 10000x120x3 array where first and second dimensions correspond to aggregated pixels and time points (months), respectively. The aggregation is done by taking the average of the pixels with available data. Along the third dimension, it contains values for TC1 (Brightness), TC2 (Greenness) and TC3 (Wetness). The occurrences of Zero (0) within this array indicate missing data at those space-time combinations.


################################################################################################################
#################################		Merge Location and TC imagery		####################################
################################################################################################################

## loading necessary packages

library(abind)

source("TC_Codes/Functions_for_TC_analysis.R")

load("Inputs_and_Outputs/TC_features_with_missing_values.Rdata")

c = 1600			## Number of cells in X coordinate in each TC feature before aggregation 
d = 1600			## Number of cells in Y coordinate in each TC feature before aggregation 
x_merge = 16		## The number of adjacent cells in TC features to be aggregated in X coordinate
y_merge = 16		## The number of adjacent cells in TC features to be aggregated in Y coordinate


## merging the TC features 

coarse_creation = coarse_data_creation(data = pixel_wise_TC_data, c = c, d = d, x_merge = x_merge, y_merge = y_merge)

save(coarse_creation, file = "Inputs_and_Outputs/coarse_selected_data_with_missing_values.Rdata")
