############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_FIA_data_analysis.R that contains all necessary functions.

## The output of TC_best_output_merge_basal_area.R that is "best_TC_model_outputs_and_basal_area_data.Rdata".

## best_TC_model_outputs_and_basal_area_data.Rdata contains a list reconstructed_TC_BALIVE consisting of five components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the complete TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 
## (iii) The fifth component is the live tree basal area in the form of a 10000x120 matrix where rows and colums correspond to pixels and time points (months), respectively. The coccurance of -99 within this matrix indicates unobserved data at those space-time combinations.


## one output of Basal_area_cv_positions_selction.R that is random_positions_for_CV_for_y1_model.Rdata
## random_positions_for_CV_for_y1_model.Rdata contains a list test_pos_list of 36 components where each component contains the postions of both zero, non-zero live tree basal area data those will be considered as test data in the correspondig holdout cross-validation method.


############################################################################################
########################		Output     #################################################
############################################################################################

## it has two outputs:

## (i) CV_output_for_y1_model.Rdata that contains the MCMC outputs of probability of forest at each test data in the form of a 540x(nmc/nthin) matrix named prob_CV_storage where rows and columns correspond to the test data and MCMC outputs of corresponding test data after burn in and thinning, respectively.

## (ii) ROC_area_for_y1_model.Rdata that contains the ROC area 

################################################################################################################
#####################     Cross Validation     #################################################################
################################################################################################################

set.seed(7777)

## Loading necessary packages

library(mvtnorm)
library(matrixStats) # for rowMedians
library(verification) # to draw ROC curve

source("FIA_Codes/Functions_for_FIA_data_analysis.R")

load("Inputs_and_Outputs/best_TC_model_outputs_and_basal_area_data.Rdata")	## loading the data 

load("Inputs_and_Outputs/random_positions_for_CV_for_y1_model.Rdata")
	
	locations_all = reconstructed_TC_BALIVE$locations	## location of the cells
	
	S = nrow(locations_all) 	## Total number of grid-cells in response

	T = 10	## Total number of time points (years here) in response

	n_TC = 3		## Total number TC features

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
	## construction of count_full matrix
	
	## Count_full is a binary matrix where 0 and 1 represents absence and presence of forest respectively.
	
	## -99 is for unavailable data 

	count_full = matrix(-99,  S,  T)

	pos_0 = which(reconstructed_TC_BALIVE$BALIVE == 0, arr.ind = TRUE)

	pos_0 = pos_0[which(pos_0[ , 2] >  0 & pos_0[ , 2] < 121), ]

	pos_1 = which(reconstructed_TC_BALIVE$BALIVE > 0, arr.ind = TRUE)

	pos_1 = pos_1[which(pos_1[ , 2] >  0 & pos_1[ , 2] < 121), ]

	for(s in 1:(dim(pos_1)[1])) count_full[pos_1[s, 1], ceiling(pos_1[s, 2]/12)] = 1

	for(s in 1:(dim(pos_0)[1])) count_full[pos_0[s, 1], ceiling(pos_0[s, 2]/12)] = 0
	
	count_full = t(count_full)
	
	pos_avail =  which(count_full != -99)

	############################################################################################################
	
	pred_mat_all = pred_mat_full[ , pos_avail]
		
	count_all = count_full[pos_avail]

	n_CV = 36		## total number of replication of holdout method 

	first_year = ((1:S) - 1)*T + 1

	n_total = length(count_all)

	l_test = 15	## number of samples that are considered as tested data
	
	count_prob_CV = matrix(0,nrow = 1, ncol = 2)
	
	prob_CV = list()
	
	for(cv in 1:n_CV)
		
	{
	
	set.seed(7777)

	pos = test_pos_list[[cv]]

	count = count_all[-c(pos)]
	
	count_test = count_all[c(pos)]

	n = length(count)
		
	pos_0_small = which(count == 0)

	pos_1_small = which(count == 1)

	upper_bound = rep(Inf, n)

	upper_bound[pos_0_small] = 0

	lower_bound = rep(-Inf, n)

	lower_bound[pos_1_small] = 0

	pred_mat = t(pred_mat_all[ ,-c(pos)])
	
	pred_mat_test = t(pred_mat_all[ ,c(pos)])
		
	############################################################################################################

	## initial values

	beta = c( qnorm(mean(count)), rep(0, (p - 1) ) )
			
	sigma2 = 1

	# Horseshoe related stuff:

	tau2_beta    = 1
	
	zai      = 1
    
	lambda2 = rep(1, (p - 1) )
	
	nu = rep(1, (p - 1) )
    
	B = c(100, tau2_beta*lambda2)
	
	########################################################################################################
	
	mu = c( crossprod(t( pred_mat ), beta) )

	########################################################################################################
	
	nit = 10000   ## number of initial runs that you want to discard
	nmc = 50000    ## number of runs that you want to do
	nthin = 10   ## you want t thin the sample at what interval to avoid correlation ?

	## during the mcmc, you need to store the simulated samples of necessary parameters, so define storage variables
		
	beta_store = matrix(0, p, nmc/nthin)
	
for (iter in 1:(nit+nmc))
	{		
	########################################################################################################
	
	## update y (latent variable for probability model) from truncated normal distribution
	
		y = rtrunc_norm (number= n, a = lower_bound, b = upper_bound, mu = mu, sigma2 = 1)
	
	########################################################################################################
	
	## update beta that follows multivariate normal distribution

		v1_inv = crossprod( pred_mat )/ sigma2
				
		v2_inv = diag(1/B)

		m1_by_v1 = crossprod( pred_mat, y ) / sigma2
				
		beta_post_var = inv_and_logdet.sym(v1_inv + v2_inv)[[2]]

		beta_post_mean = crossprod(beta_post_var, m1_by_v1)

		beta = c( rmvnorm(1, beta_post_mean, beta_post_var) )
	
	########################################################################################################
	
	## Horseshoe related parameters update

	##  while simulating beta parameters, their prior mean is 0 and prior dispersion is diag(B)

	## after you simulated the beta, simulate horseshoe related parameters
        
		tau2_beta = 1/rgamma( 1,shape = p/2, rate = 1/zai + 0.5*sum((beta[-1]^2.0)/lambda2) )
        
		zai = 1/rgamma(1, shape = 1, rate = 1 + 1/tau2_beta)

		lambda2 = 1/rexp((p - 1), rate = 1/nu + 0.5*(beta[-1]^2.0)/tau2_beta)
        
		nu = 1/rexp( (p - 1), rate = 1 + 1/lambda2)
    
		B[-1] = tau2_beta*lambda2

	########################################################################################################
	
	## update mean

		mu = c( crossprod(t( pred_mat ), beta) )

		if (iter > nit & (iter - nit)%%nthin==0)

		{
			
			beta_store[ , (iter - nit)/nthin] = beta

		}

	}


	mu_test = pred_mat_test%*%beta_store
		
	prob_test_store = pnorm(mu_test)		

	output = cbind(count_test, rowMedians(prob_test_store))

	count_prob_CV = rbind(count_prob_CV, output)

	prob_CV[[cv]] = prob_test_store
	
	}


	prob_CV_storage = do.call(rbind, prob_CV)
	
	save(prob_CV_storage, file = "Inputs_and_Outputs/CV_output_for_y1_model.Rdata")	## this is needed to calculate basal area cross validation outputs

	count_prob_CV = count_prob_CV[-1,]	## count and probability storage
	
	ROC_area = round(roc.area(count_prob_CV[,1], count_prob_CV[,2])$A, digits = 3)	## ROC area 

	print(paste0("ROC area = " , ROC_area))
	
	save(ROC_area, file = "Inputs_and_Outputs/ROC_area_for_y1_model.Rdata")	## this is needed to calculate basal area cross validation outputs
