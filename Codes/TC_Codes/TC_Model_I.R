############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_TC_analysis.R that contains all necessary functions.

## the output of Aggregation.R that is "coarse_selected_data_with_missing_values.Rdata"
## "coarse_selected_data_with_missing_values.Rdata" contains a list coarse_creation that consists of two components:

## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.

## (ii) The second component contains the spatiotemporal data in form of a 10000x120x3 array where first and second dimensions correspond to pixels and time points (months), respectively. Along the third dimension, it contains values for TC1 (Brightness), TC2 (Greenness) and TC3 (Wetness). The occurrences of Zero (0) within this array indicate missing data at those space-time combinations.


## in the code: 
## data: TC variable
## covariate: U matrix at any cell

############################################################################################
########################		Output     #################################################
############################################################################################

## the output is "TC_tc_Model_I_outputs.Rdata" where tc should be replaced by appropriate TC number from 1, 2 and 3 that is a list fill_up_missing_data containing four components:

## (i) the first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.

## (ii) the second component contains the complete spatiotemporal map of corresponding TC data, after filling in the missing observations, in form of a 10000x120 matrix where first and second dimensions correspond to pixels and time points (months), respectively. 

## (iii) the third component is a vector with predictive uncertainty

## (iv) the fourth component is a vector of size 2 with log likelihood and BPIC, respectively. 

## (v) the fifth component is a 92x3 matrix whose first column represents the positions of the months at when more than 50% data are available. The second and third columns consist of corresponding month's Moran's I and p-value of the right tail test, respectively.

################################################################################################################
#####################     Data Analysis     ####################################################################
################################################################################################################

set.seed(7777)

## Loading necessary packages

library(mvtnorm)
library(fields)	## for rdist in neighborhood matrix construction 
library(ape)
library(coda)

source("TC_Codes/Functions_for_TC_analysis.R")

load("Inputs_and_Outputs/coarse_selected_data_with_missing_values.Rdata") ## loading output of Aggregation.R that contains a list coarse_creation

TC = as.integer(readline(" input 1 for TC1, 2 for TC2, 3 for TC3 "))	## choose the TC covariate to model

## data representing the corresponding TC measurement chosen in the previous line where rows and columns represents grid-cells and months respectively

data = coarse_creation$coarse_data[ , , TC]

## In data 0 represents missing value 

##########################################################################################

## To translate and rescale each TC features to the same range as that of TC1 

data_TC1 = coarse_creation$coarse_data[ , , 1] ## data_TC1 is the TC1 measurement

pos_0_cell = apply(data, 1, function(x) which(x == 0))  ## list of size S where each list contains all missing months at the corresponding cell

pos_0_cell_length = sapply( pos_0_cell, length)  ## vector of size S containing the number of missing months at the corresponding cell

pos_0_month = apply(data, 2, function(x) which(x == 0))  ## list of size T where each list contains all missing cells at the corresponding month

pos_0_month_length = sapply(pos_0_month, length)  ## vector of size T containing the number of missing cells at the corresponding month

pos_0 = which(data == 0)  ## vector of the positions of all missing data 

## following available data are translated and rescaled to the same range as that of TC1, and then log transformation is performed 

data[-pos_0] = log( data[-pos_0]* diff(range(data_TC1[-pos_0]))/diff(range(data[-pos_0])) - 
			   min( data[-pos_0]* diff(range(data_TC1[-pos_0]))/diff(range(data[-pos_0])) ) + min(data_TC1[-pos_0]) )

############################################################################################

location = coarse_creation$coarse_location ## location of the cells 

S = dim(location)[1] ## Total number of grid-cells in the data

T = dim(data)[2]  ## Total number of months in the data

n_total = S*T ## Total number of data points

n_non_missing = n_total - length(pos_0) ## Total number of non-missing data points

covariate = t( sapply(1:T, function(t)  c(1, cos((1:2)*pi*t/6), sin((1:2)*pi*t/6))) ) ## Fourier series with first two harmonics at any cell 

p_beta = dim(covariate)[2]
 
X_train = covariate[rep(1:T, (S - pos_0_month_length)),] ## Construction of full covariate matrix for available TC feature.  

y_train = data[-pos_0] ## vector of available TC feature

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

e = y_train - colSums(t(X_train) * beta) ## residual

sigma2 = var(c(e))

######################################################################################

## how many times the MCMC is run 

nit = 10000     ## number of initial runs that are discarded
nmc = 50000    	## number of runs that author wants to do
nthin = 50   	##  thinning interval to avoid correlation in the sample

## during the MCMC, author needs to store the simulated parameters, so define storage variables

beta_store = matrix(0, p_beta, nmc/nthin)

sigma2_store = c()

nlog_likelihood_store = c() ## negative log-likelihood storage

pred_store = matrix(0, length(pos_0), nmc/nthin) ## predicting the missing values

## start the MCMC loop

print(date())

for (iter in 1:(nit+nmc))
	{
	#################################################################################
	
	## posterior distribution of beta is Multivariate Normal

		beta_post_dispersion = inv_and_logdet.sym((XTX/sigma2)+ diag(p_beta)/c0)[[2]]   

		beta_post_mean =  colSums(beta_post_dispersion * ((XTy/sigma2) + (beta0/c0)))         
    
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
    
    if (iter > nit & (iter - nit)%%nthin==0)
		{
		
			beta_store[ , (iter - nit)/nthin] =  beta
			
			sigma2_store[ (iter - nit)/nthin ] = sigma2
			
			pred_store[ , (iter - nit)/nthin] = - residual[pos_0] + rnorm(length(pos_0), 0, sqrt(sigma2))
			
			nlog_likelihood_store[(iter - nit)/nthin] = n_non_missing*log(sigma2)/2 + sum((residual[-pos_0])^2)/(2*sigma2)
			
		}

	if(iter == 1000) print(date())	
	}

print(date())

	##########################################################################################

	## Complete the TC features by filling in the missing observations
	
if(TC == 1) 
    {
        data[pos_0] = apply(pred_store, 1, median)
        Transformed_TC1 = data
		
		## Predictive uncertainty (90% HPD interval width) for TC1  
		
		pred_HPD_interval = apply(pred_store, 1, function(x) HPDinterval(mcmc(x), prob = 0.90)[1,] )
		uncertainty_TC1 = c(diff(pred_HPD_interval))
	
	fill_up_missing_data = list(locations = coarse_creation$coarse_location, Transformed_TC1 = Transformed_TC1, uncertainty_TC1 = uncertainty_TC1)
    }

if(TC == 2) 
    {
        data[pos_0] = apply(pred_store, 1, median)
        Transformed_TC2 = data

		## Predictive uncertainty (90% HPD interval width) for TC2  
		
		pred_HPD_interval = apply(pred_store, 1, function(x) HPDinterval(mcmc(x), prob = 0.90)[1,] )
		uncertainty_TC2 = c(diff(pred_HPD_interval))

	fill_up_missing_data = list(locations = coarse_creation$coarse_location, Transformed_TC2 = Transformed_TC2, uncertainty_TC2 = uncertainty_TC2)
    }

if(TC == 3) 
    {
        data[pos_0] = apply(pred_store, 1, median)
        Transformed_TC3 = data

		## Predictive uncertainty (90% HPD interval width) for TC3  
		
		pred_HPD_interval = apply(pred_store, 1, function(x) HPDinterval(mcmc(x), prob = 0.90)[1,] )
		uncertainty_TC3 = c(diff(pred_HPD_interval))
		 
	fill_up_missing_data = list(locations = coarse_creation$coarse_location, Transformed_TC3 = Transformed_TC3, uncertainty_TC3 = uncertainty_TC3)
    }

	##########################################################################################
	
	# Log likelihood and BPIC calculation
	
	LL = median(-nlog_likelihood_store) ## Log likelihood
	
	D_bar = 2*mean(nlog_likelihood_store)

	mean_beta = apply(beta_store, 1, mean)
		
	mean_X_beta = colSums(t(covariate) * mean_beta)
	
	mean_res = sweep(data, 2, mean_X_beta)

	mean_sigma2 = mean(sigma2_store)
	
	D_theta_bar_mean = 2*(n_non_missing*log(mean_sigma2)/2 + sum((mean_res[-pos_0])^2)/(2*mean_sigma2))
	
	BPIC = 3*D_bar - 2*D_theta_bar_mean
	
	#############################################################################################################
	
	## Moran's I  test to calculate spatial association among residual

	## Construction of neighborhood Matrix

	x_loc = sort(unique(location[,1]))

	x_dist= x_loc[2] - x_loc[1]

	y_loc = sort(unique(location[,2]))

	y_dist= y_loc[2] - y_loc[1]

	W = sapply(1:S,f_W)	## list of size S that contains the neighbors of the corresponding cell

	W_mat = matrix(0, S, S)	## Neighborhood Matrix

	for(i in 1:S)
		{
			W_mat[i, W[[i]]] = 1
		}

	vec = which((pos_0_month_length/S) < 0.5)  ## months which have more than 50% available data. Moran's I are calculated only for these months 
	
	residual_store = array(0, dim = c(S, T, nmc/nthin))
	
	for (i in 1:(nmc/nthin) ) residual_store[ , , i] = sweep(data , 2 , colSums(t(covariate) * beta_store[ , i]) ) 
	
	data_na = matrix(0, S ,T)

	data_na[pos_0] = NA 

	med_residual = apply(residual_store, c(1,2), median)
	
	med_residual = med_residual + data_na ## rename the missing data as NA instead of 0

	moran_test = sapply(vec, function(i){unlist(Moran.I(med_residual[,i], W_mat, scaled = FALSE, na.rm = TRUE, alternative = "greater"))})

	#############################################################################################################################################
	
	## saving all necessary outputs to create first two rows in Table 1
	
	fill_up_missing_data$LL_BPIC = round(c(LL, BPIC), 3)
	
	fill_up_missing_data$moran_I = cbind(vec, round(t(moran_test)[ , c(1, 4)],3) )
	
	save(fill_up_missing_data, file = paste0("Inputs_and_Outputs/Model_I_TC_", TC ,"_outputs.Rdata"))