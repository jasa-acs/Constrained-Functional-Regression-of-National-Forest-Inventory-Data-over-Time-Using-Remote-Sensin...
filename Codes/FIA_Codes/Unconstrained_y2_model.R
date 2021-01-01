############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_FIA_data_analysis.R that contains all necessary functions.

## The output of TC_best_output_merge_basal_area.R that is "best_TC_model_outputs_and_basal_area_data.Rdata".

## best_TC_model_outputs_and_basal_area_data.Rdata contains a list reconstructed_TC_BALIVE consisting of five components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the complete TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 
## (iii) The fifth component is the live tree basal area in the form of a 10000x120 matrix where rows and colums correspond to pixels and time points (months), respectively. The coccurance of -99 within this matrix indicates unobserved data at those space-time combinations.

############################################################################################
########################		Output     #################################################
############################################################################################

## the output is "unconstrained_y2_model_output.Rdata" that contains a list named outputs consisting of two components:
## (i) the first component is the posterior median likelihoods at each available live tree basal area measurement.
## (ii) the second component is a vector of size 2 with log likelihood and BPIC, respectively. 

################################################################################################################
#####################     Data Analysis     ####################################################################
################################################################################################################

set.seed(7777)

## Loading necessary packages

library(mvtnorm)
library(matrixStats) # for colMedians

source("FIA_Codes/Functions_for_FIA_data_analysis.R")

load("Inputs_and_Outputs/best_TC_model_outputs_and_basal_area_data.Rdata")	## loading the data 

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

	pred_mat_full = t(rbind(1, pred_mat_full))
	
	p = ncol(pred_mat_full)
	
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

	## covariates and response for available response
	
	pred_mat = pred_mat_full[pos_avail, ]
	
	y = y_full[pos_avail]

	n = length(y)	## total number of available response
	
	n_UA = S*T - n	## total number of unavailable response
		
	#####################################################################################

	## set priors

	## prior of sigma2 is IG(a0,b0)

	a0 = 2.000001

	b0 = 1.000001

	#####################################################################################

	## Initial values

	beta = c(mean(y), rep(0, (p-1)) )
				
	sigma2 = var(y) 

	# Horseshoe related parameters:

	tau2_beta    = 1
	
	zai      = 1
    
	lambda2 = rep(1, (p - 1))
    
	nu      = rep(1, (p - 1))
    
	B = c(100, tau2_beta*lambda2)

	###########################################################################################

	mu = c( crossprod( t(pred_mat), beta) )

	#####################################################################################

	nit = 10000   ## number of initial runs that you want to discard
	nmc = 50000     ## number of runs that you want to do
	nthin = 10   ## you want t thin the sample at what interval to avoid correlation ?

	## during the mcmc, you need to store the simulated samples of necessary parameters, so define storage variables
			
	beta_store = matrix(0, p, nmc/nthin)
		
	sigma2_store = c()
	
	nlog_likelihood_store = c()	## negative log likelihood

	unconstrained_likelihood_store = matrix(0, nmc/nthin, n)	## likelihood for each available response
	
print(date())

for (iter in 1:(nit+nmc))
	{

	#####################################################################################

	## update beta that follows multivariate normal distribution

		v1_inv = crossprod( pred_mat )/ sigma2
				
		v2_inv = diag(1/B)

		m1_by_v1 = crossprod( pred_mat, y ) / sigma2
				
		beta_post_var = inv_and_logdet.sym(v1_inv + v2_inv)[[2]]

		beta_post_mean = crossprod(beta_post_var, m1_by_v1)
		
		beta = c( rmvnorm(1, beta_post_mean, beta_post_var) )

	#####################################################################################

	## Horseshoe related parameters update

	## while simulating beta parameters, their prior mean is 0 and prior dispersion is diag(B)

	## after you simulated the beta, simulate horseshoe related parameters
        
		tau2_beta = 1/rgamma( 1, shape = p/2, rate = 1/zai + 0.5*sum((beta[-1]^2.0)/lambda2) )
        
		zai = 1/rgamma(1, shape = 1, rate = 1 + 1/tau2_beta)

		lambda2 = 1/rexp((p - 1), rate = 1/nu + 0.5*(beta[-1]^2.0)/tau2_beta)
        
		nu = 1/rexp( (p - 1), rate = 1 + 1/lambda2)
    
		B[-1] = tau2_beta*lambda2

	#####################################################################################
		
	## update mean 

		mu = c( crossprod(t( pred_mat), beta) )

	#####################################################################################

	## update sigma2 that follows Inverse Gamma distribution

		S2 = (y - mu)^2
		   
		sigma2_post_shape = a0 + n/2 
		
		sigma2_post_rate = b0 + 0.5*sum(S2)

		sigma2 = 1/rgamma(1, shape = sigma2_post_shape, rate = sigma2_post_rate)

	if (iter > nit & (iter - nit)%%nthin==0)
		{
					
			beta_store[ , (iter - nit)/nthin] = beta
			
			sigma2_store[(iter - nit)/nthin] = sigma2 

			nlog_likelihood_store[(iter - nit)/nthin] = 0.5 * n * log(sigma2) + 0.5 * sum(S2)/sigma2

			unconstrained_likelihood_store[(iter - nit)/nthin, ] = exp(-(0.5 * log(sigma2) + 0.5 * S2/sigma2))
			
		}

		if(iter == 1000) print(date())
		
	}

print(date())

	#########################################################################################################
	
	# Log likelihood and BPIC calculation
	
	LL = median(-nlog_likelihood_store) ## Log likelihood

	mean_beta = apply(beta_store, 1, mean)
	
	mean_sigma2 = mean(sigma2_store)
	
	mean_mu = c( crossprod(t( pred_mat), mean_beta) )

	mean_S2 = (y - mean_mu)^2
		
	mean_nlog_likelihood = 0.5 * n * log(mean_sigma2) + 0.5 * sum(mean_S2)/mean_sigma2

	D_bar = 2*mean(nlog_likelihood_store)
	
	D_theta_bar_mean = 2*mean_nlog_likelihood
	
	BPIC = 3*D_bar - 2*D_theta_bar_mean

	########################################################################################################
	
	unconstrained_likelihood = colMedians(unconstrained_likelihood_store)
	
	#########################################################################################################
	#########################################################################################################

	## saving necessary outputs

	outputs = list(unconstrained_likelihood = round(unconstrained_likelihood,3), LL_BPIC = round(c(LL, BPIC), 3) )  
	
	save(outputs, file = "Inputs_and_Outputs/unconstrained_y2_model_output.Rdata")