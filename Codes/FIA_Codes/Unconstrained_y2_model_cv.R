############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_FIA_data_analysis.R that contains all necessary functions.

## The output of TC_best_output_merge_basal_area.R that is "best_TC_model_outputs_and_basal_area_data.Rdata".

## best_TC_model_outputs_and_basal_area_data.Rdata contains a list reconstructed_TC_BALIVE consisting of five components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the complete TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 
## (iii) The fifth component is the live tree basal area in the form of a 10000x120 matrix where rows and colums correspond to pixels and time points (months), respectively. The coccurance of -99 within this matrix indicates unobserved data at those space-time combinations.


## one output of Basal_area_cv_positions_selction.R that is random_positions_for_CV_for_y2_model.Rdata
## random_positions_for_CV_for_y2_model.Rdata contains a list test_pos_list of 36 components where each component contains only the postions of non-zero live tree basal area data those will be considered as test data in the correspondig holdout cross-validation method.


############################################################################################
########################		Output     #################################################
############################################################################################

## the output is "CV_outputs_unconstrained_y2_model.csv" that is a 3x1 matrix containing mean of absolute bias, uncertainty and empirical coverage, respectively.

################################################################################################################
#####################     Cross Validation     #################################################################
################################################################################################################

set.seed(7777)

library(mvtnorm)
library(coda)	## for HPD interval

source("FIA_Codes/Functions_for_FIA_data_analysis.R")

load("Inputs_and_Outputs/best_TC_model_outputs_and_basal_area_data.Rdata")	## loading the data 

load("Inputs_and_Outputs/random_positions_for_CV_for_y2_model.Rdata")

	locations_all = reconstructed_TC_BALIVE$locations	## location of the cells 

	S = nrow(locations_all)	## Total number of grid-cells in response

	T = 10	## Total number of time points (years here) in response

	n_TC = 3	## Total number TC features

	n_total = S*T	## Total number of data points in response

	###############################################################################################################
	
	## construction of functional X matrix

	temp_col = ((1:T) - 1)*12 + 1
	
	col_no = c(sapply(0:11, function(i) temp_col + i))

	data = array(-999, c(S, 12*T, n_TC))
	
	for(TC in 1:n_TC) data[ , , TC] = reconstructed_TC_BALIVE[[TC+1]]

	data = data[ , col_no, ]

	dim(data) = c(S, T, 3*12)

	data = aperm(data, c(3, 2, 1))
		
	pred_mat_full = data
	
	dim(pred_mat_full) = c(3*12, T*S)
	
	## standardization of X matrix 

	center = apply(pred_mat_full, 1, mean)

	scale = apply(pred_mat_full, 1, sd)

	for (i in 1: length(center)) pred_mat_full[i , ] = (pred_mat_full[i , ] - center[i])/scale[i]

	pred_mat_full = rbind(1, pred_mat_full)
	
	p = nrow(pred_mat_full)
	
	#########################################################################################################

	## response 
	## here -99 represents data is not available
	
	y_full = matrix(-99,  S,  T)

	pos_1 = which(reconstructed_TC_BALIVE$BALIVE > 0, arr.ind = TRUE)

	pos_1 = pos_1[which(pos_1[,2]> 0 & pos_1[,2]<121),]

	for(s in 1:(dim(pos_1)[1])) y_full [pos_1[s,1] , ceiling(pos_1[s,2]/12)] = log(reconstructed_TC_BALIVE$BALIVE[pos_1])[s]

	y_full = t(y_full)
		
	pos_avail =  which(y_full != -99)
	
	######################################################################################################

	pred_mat_all = pred_mat_full[ , pos_avail]
		
	y_all = y_full[pos_avail]

	n_CV = 36		## total number of replication of holdout method 
		
	outputs = matrix(0, 1, 4)

	n_total = length(y_all)

	l_test = 10	## number of samples that are considered as tested data

	for(cv in 1:n_CV)
		
	{
	
set.seed(7777)

	pos = test_pos_list[[cv]]

	y = y_all[-c(pos)]
	
	y_test = y_all[c(pos)]	## test data
	
	n = length(y)

	pred_mat = t(pred_mat_all[ ,-c(pos)])
	
	pred_mat_test = t(pred_mat_all[ ,c(pos)])

	###########################################################################################

	## set priors
	
	## prior of sigma2 is IG(a0,b0)

	a0 = 2.000001

	b0 = 1.000001

	###########################################################################################

	## initial values

	beta = c(mean(y), rep(0, (p-1)) )
				
	sigma2 = var(y) 

	# Horseshoe related stuff:

	tau2_beta    = 1
	
	zai      = 1
    
	lambda2 = rep(1, (p - 1))
    
	nu      = rep(1, (p - 1))
    
	B = c(100, tau2_beta*lambda2)

	###########################################################################################

	mu = c( crossprod( t(pred_mat), beta) )

	######################################################################################################

	nit = 10000   ## number of initial runs that you want to discard
	nmc = 50000     ## number of runs that you want to do
	nthin = 10   ## you want t thin the sample at what interval to avoid correlation ?

	#################################################################################################
	
	## during the mcmc, you need to store necessary parameters, so define storage variable
			
	beta_store = matrix(0, p, nmc/nthin)
		
	sigma2_store = c()

	#################################################################################################

for (iter in 1:(nit+nmc))
	{
		
	#################################################################################################
	
	## update beta that follows multivariate normal distribution
	
		v1_inv = crossprod( pred_mat )/ sigma2
				
		v2_inv = diag(1/B)

		m1_by_v1 = crossprod( pred_mat, y ) / sigma2
				
		beta_post_var = inv_and_logdet.sym(v1_inv + v2_inv)[[2]]

		beta_post_mean = crossprod(beta_post_var, m1_by_v1)
		
		beta = c( rmvnorm(1, beta_post_mean, beta_post_var) )

	#################################################################################################
	
	# Horseshoe related parameters update

	##  while simulating beta parameters, their prior mean is 0 and prior dispersion is diag(B)

	## after you simulated the beta, simulate horseshoe related parameters
        
		tau2_beta = 1/rgamma( 1,shape = p/2, rate = 1/zai + 0.5*sum((beta[-1]^2.0)/lambda2) )
        
		zai = 1/rgamma(1, shape = 1, rate = 1 + 1/tau2_beta)

		lambda2 = 1/rexp((p - 1), rate = 1/nu + 0.5*(beta[-1]^2.0)/tau2_beta)
        
		nu = 1/rexp( (p - 1), rate = 1 + 1/lambda2)
    
		B[-1] = tau2_beta*lambda2

	#################################################################################################
	
	## update mean

		mu = c( crossprod( t( pred_mat ), beta) )

	#################################################################################################
	
	## update sigma2 that follows Inverse Gamma distribution

		S2 = (y - mu)^2
		   
		sigma2_post_shape = a0 + n/2 
		
		sigma2_post_rate = b0 + 0.5*sum(S2)

		sigma2 = 1/rgamma(1, shape = sigma2_post_shape, rate = sigma2_post_rate)

	##################################################################################################

	if (iter > nit & (iter - nit)%%nthin==0)
		{
		
			
			beta_store[ , (iter - nit)/nthin] = beta
			
			sigma2_store[(iter - nit)/nthin] = sigma2 

		}

	}

	mu_test = pred_mat_test%*%beta_store
	
	y_pred_MCMC_temp = mu_test + matrix( rep(sigma2_store, each = l_test)*rnorm(l_test*nmc/nthin) , l_test, nmc/nthin)

	y_pred_MCMC = exp(y_pred_MCMC_temp)

	y_pred = apply(y_pred_MCMC, 1, median)

	CI_HPD = apply(y_pred_MCMC, 1, function(x) HPDinterval(mcmc(x), prob = 0.90)[1,] )
	
	outputs = rbind(outputs, cbind(exp(y_test), y_pred, CI_HPD[1, ], CI_HPD[2, ]))

}

## the given unit of basal area is ft^2/acre. To transform it into m^2/hectare it is needed to mulitply the basal area by 0.229568411	

	outputs = outputs[-1, ]*0.229568411	## transforming  ft^2/acre  to m^2/hectare 

	colnames(outputs) = c("y_test", "y_pred",  "HPD LI", "HPD UI")
	
	abs_bias = mean(abs(outputs[ , 1] - outputs[ , 2]))

	uncertainty = mean(outputs[ , 4] - outputs[ , 3])

	empirical_coverage = sum(as.integer(outputs[ , 1] >= outputs[ , 3] & outputs[ , 1] <= outputs[ , 4]))/(l_test*n_CV)

	CV_outputs = round( matrix(c(abs_bias, uncertainty, empirical_coverage), 3, 1), 3)
	
	write.csv(CV_outputs, file = "Inputs_and_Outputs/CV_outputs_unconstrained_y2_model.csv", row.names = F)
