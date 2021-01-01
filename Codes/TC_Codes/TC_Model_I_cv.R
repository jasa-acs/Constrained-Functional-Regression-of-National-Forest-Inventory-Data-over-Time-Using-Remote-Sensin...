############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_TC_analysis.R that contains all necessary functions.

## the output of Aggregation.R that is "coarse_selected_data_with_missing_values.Rdata"
## "coarse_selected_data_with_missing_values.Rdata" contains a list coarse_creation that consists of two components:

## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.

## (ii) The second component contains the spatiotemporal data in form of a 10000x120x3 array where first and second dimensions correspond to pixels and time points (months), respectively. Along the third dimension, it contains values for TC1 (Brightness), TC2 (Greenness) and TC3 (Wetness). The occurrences of Zero (0) within this array indicate missing data at those space-time combinations.




## the output of TC_cv_positions_selection.R that is random_positions_for_CV.Rdata
## random_positions_for_CV.Rdata contains a list test_pos_list_all that consists of 36 combinations where each component contains the positions of test data for the corresponding holdout method



## in the code 
## data_all: TC variable
## covariate: U matrix at any cell

############################################################################################
########################		Output     ##################################################
############################################################################################

## the output is "CV_outputs_Model_I_TC_tc.csv" where tc should be replaced by appropriate TC number from 1, 2 and 3 that is a 3x1 matrix containing mean of absolute bias, uncertainty and empirical coverage, repspectively.


############################################################################################
##########################     Cross validation part       #################################
############################################################################################

set.seed(7777)

## Loading necessary packages

library(mvtnorm)
library(coda)	## for HPD interval

source("TC_Codes/Functions_for_TC_analysis.R")

l_test = 25000	# number of samples that are considered as tested data

n_CV =  36	## total number of replication of holdout method 

load("Inputs_and_Outputs/coarse_selected_data_with_missing_values.Rdata") ## loading output of Aggregation.R that contains a list coarse_creation

load("Inputs_and_Outputs/random_positions_for_CV.Rdata")	# loading the output of TC_cv_positions_selection.R

TC = as.integer(readline(" input 1 for TC1, 2 for TC2, 3 for TC3 "))	## choose the TC covariate to model

## data representing the corresponding TC measurement chosen in the previous line where rows and columns represents grid-cells and months respectively

data_all = coarse_creation$coarse_data[ , , TC]

data_temp = data_all

## In data 0 represents missing value 

##########################################################################################

## To translate and rescale each TC features to the same range as that of TC1 

data_TC1 = coarse_creation$coarse_data[ , , 1]	## data_TC1 is the TC1 measurement

pos_0_all = which(data_all == 0)	## vector of the positions of all missing data 

## following available data are translated and rescaled to the same range as that of TC1, and then log transformation is performed 

data_all[-pos_0_all] = log( data_all[-pos_0_all]* diff(range(data_TC1[-pos_0_all]))/diff(range(data_all[-pos_0_all])) - 
			   min( data_all[-pos_0_all]* diff(range(data_TC1[-pos_0_all]))/diff(range(data_all[-pos_0_all]))) + min(data_TC1[-pos_0_all]))

location = coarse_creation$coarse_location	## location of the cells 

S = dim(location)[1] ## Total number of grid-cells in the data

T = dim(data_all)[2]  ## Total number of months in the data

n_total = S*T ## Total number of data points

covariate = t(sapply(1:T, function(t)  c(1, cos((1:2)*pi*t/6), sin((1:2)*pi*t/6))))	## Fourier series with first two harmonics at any cell 

p_beta = dim(covariate)[2]

cross_validation = list()

for (cv in 1:n_CV)

{

set.seed(7777)

data = data_all

test_pos = test_pos_list_all[[cv]]

test_data = data[test_pos]	## test data

data[test_pos] = 0	## treating test data as missing data

pos_0_cell = apply(data, 1, function(x) which(x == 0))	## list of size S where each list contains all missing months at the corresponding cell

pos_0_cell_length = sapply(pos_0_cell, length)	## vector of size S containing the number of missing months at the corresponding cell

pos_0_month = apply(data, 2, function(x) which(x == 0))	## list of size T where each list contains all missing cells at the corresponding month

pos_0_month_length = sapply(pos_0_month, length)	## vector of size T containing the number of missing cells at the corresponding month

pos_0 = which(data == 0)	## vector of the positions of all missing data 

n_non_missing = n_total - length(pos_0)	## Total number of non-missing data points

X_train = covariate[rep(1:T, (S - pos_0_month_length)), ]	## Construction of full covariate matrix for available TC feature.  

y_train = data[-pos_0]	## vector of available TC feature

XTX = crossprod(X_train)

XTy = colSums(X_train*y_train)

#####################################################################################

## set priors and initial values

## beta ~ MVN(beta0,c0*I)

beta0 = rep(0, p_beta)
c0 = 1000

## prior for sigma2 is IG(a10,b10)

a10 = 2.000001
b10 = 1.000001

#####################################################################################

## Initial values from least-square method

beta = colSums((inv_and_logdet.sym(XTX)[[2]]) * XTy)

e = y_train - colSums(t(X_train) * beta)

sigma2 = var(c(e))

#####################################################################################

## how many times you want to run the mcmc

nit = 10000      ## number of initial runs that you want to discard
nmc = 50000      ## number of runs that you want to do
nthin = 50    ## you want t thin the sample at what interval to avoid correlation ?

## during the mcmc, you need to store the simulated test data, so define storage variable

pred_test_store = matrix(0, l_test, nmc/nthin)	## store simulated test data

## start the mcmc loop

for (iter in 1:(nit+nmc))
	{
	
	#################################################################################
	
	## posterior distribution of beta is Multivariate Normal

		beta_post_dispersion = inv_and_logdet.sym((XTX/sigma2)+ diag(p_beta)/c0)[[2]]   

		beta_post_mean =  colSums(beta_post_dispersion * ((XTy/sigma2)+(beta0/c0)))         
    
	## then simulate sample of beta
    
		beta = c(rmvnorm (1, beta_post_mean, beta_post_dispersion))

	#################################################################################

		X_beta = colSums(t(covariate) * beta)
	
		residual = t(t(data) - X_beta)

	#################################################################################
		
	## posterior distribution of sigma2 is Inverse Gamma
	
		sigma2_post_shape = a10 + n_non_missing/2
	
		sigma2_post_rate = b10 + sum(residual[-pos_0]^2)/2
    
	## then simulate sample of sigma2
    
		sigma2 = 1/rgamma(1, shape = sigma2_post_shape, rate = sigma2_post_rate)	

	#################################################################################
		
	## Now store the predicted test data 
    
	if (iter > nit & (iter - nit)%%nthin==0)
		{

			pred_test_store[ , (iter - nit)/nthin] = - residual[test_pos] + rnorm(l_test, 0, sqrt(sigma2))			
		}
}

cross_validation[[cv]] = cbind(test_data, pred_test_store)

}

TC_test_pred = do.call(rbind, cross_validation)

## transforming the TC data into original scale

TC_test_pred_original = (exp(TC_test_pred) + min( data_temp[-pos_0_all]* diff(range(data_TC1[-pos_0_all]))/diff(range(data_temp[-pos_0_all]))) - min(data_TC1[-pos_0_all])) * diff(range(data_temp[-pos_0_all]))/diff(range(data_TC1[-pos_0_all]))

rm(TC_test_pred)

## Calculating necessary measures for cross validation

TC_test = TC_test_pred_original[ , 1]

TC_pred = apply(TC_test_pred_original[ , -1], 1, median)

TC_CI_HPD = apply(TC_test_pred_original[ , -1], 1, function(x) HPDinterval(mcmc(x), prob = 0.90)[1,] )

abs_bias = mean(abs(TC_test - TC_pred))

uncertainty_HPD = mean( diff(TC_CI_HPD) )

empirical_coverage_HPD = mean(as.numeric(TC_test >= TC_CI_HPD[1, ] & TC_test <= TC_CI_HPD[2, ]))

CV_outputs = round( matrix(c(abs_bias, uncertainty_HPD, empirical_coverage_HPD), 3,1),  3 )

write.csv(CV_outputs, file = paste0("Inputs_and_Outputs/Model_I_TC_", TC, "_CV_outputs.csv"), row.names = FALSE)