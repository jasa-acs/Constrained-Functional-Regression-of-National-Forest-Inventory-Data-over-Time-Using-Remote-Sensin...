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

## the output is "constrained_y2_model_output.Rdata" that contains a list named outputs consisting of three components:
## (i) the first component is the posterior median likelihoods at each available live tree basal area measurement.
## (ii) the second component is a vector of size 2 with log likelihood and BPIC, respectively. 
## (iii) the third component is the MCMC outputs of predicted live tree basal area measurement of the form of a TxSx(nmc/nthin) array where first and second dimensions correspond to time points (years) and pixels, respectively. Along the third dimension it contains the corresponding MCMC outputs after burn in and thinning.

################################################################################################################
#####################     Data Analysis     ####################################################################
################################################################################################################

set.seed(7777)

## Loading necessary packages

library(mvtnorm)
library(matrixStats) # for colMedians
library(abind)

source("FIA_Codes/Functions_for_FIA_data_analysis.R")

load("Inputs_and_Outputs/best_TC_model_outputs_and_basal_area_data.Rdata")	## loading the data 

	locations_all = reconstructed_TC_BALIVE$locations	## location of the cells 

	S = nrow(locations_all)	## Total number of grid-cells in response

	T = 10	## Total number of time points (years here) in response

	n_TC = 3	## Total number TC features

	n_total = S*T	## Total number of data points in response

	###############################################################################################################

	## Construction of functional X matrix
	
	temp_col = ((1:T) - 1)*12 + 1
	
	col_no = c(sapply(0:11, function(i) temp_col + i))

	data = array(-999, c(S, 12*T, n_TC))
	
	for(TC in 1:n_TC) data[ , , TC] = reconstructed_TC_BALIVE[[TC+1]]

	data = data[ , col_no, ]

	dim(data) = c(S, T, 3* 12)

	data = aperm(data, c(3, 2, 1))
		
	pred_mat_full = data
	
	dim(pred_mat_full) = c(3*12, T*S)

	## standardization of X matrix 
	
	center = apply(pred_mat_full, 1, mean)

	scale = apply(pred_mat_full, 1, sd)

	for (i in 1: length(center)) pred_mat_full[i , ] = (pred_mat_full[i , ] - center[i])/scale[i]

	pred_mat_full = rbind(1, pred_mat_full)
	
	p = nrow(pred_mat_full)

	pred_mat_full_3D = pred_mat_full
		
	dim(pred_mat_full_3D) = c( p, T, S)
	
	##############################################################################################

	y_full = matrix(-99,  S,  T)	## y_full is same as y^(2) in the manuscript

	pos_1 = which(reconstructed_TC_BALIVE$BALIVE > 0, arr.ind = TRUE)

	pos_1 = pos_1[which(pos_1[,2]> 0 & pos_1[,2]<121),]

	for(s in 1:(dim(pos_1)[1])) y_full [pos_1[s, 1] , ceiling(pos_1[s, 2]/12)] = log(reconstructed_TC_BALIVE$BALIVE[pos_1])[s]

	available_locations_all = list()

	for (t in 1:T) available_locations_all[[t]] =   sort(which(y_full[,t] != -99)) 

	available_locations_unique = sort(unique(unlist(available_locations_all)) ) ## this stores all the cell locations, where any nonzero data was available during 10 years

	S_A =  length(available_locations_unique)

	unavailable_locations_unique = setdiff(1:S, available_locations_unique)

	S_UA =  S - S_A

	y_full = t(y_full)
	
	pos_avail_full =  which(y_full != -99)	## position of all available non-zero basal area data 

	############################################################################################################################
	
	## Some calculations to update z from all possible combination of z using multinomial distribution

	z_all = as.matrix( unname( cbind( 1, expand.grid(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1)) ) ) )
	
	n_z = nrow(z_all)
	
	Bt_all = t( apply(z_all, 1, Bt_matrix_create) )
	
	pred_mat_3D = pred_mat_full_3D[ , , available_locations_unique]
	
	## Construction of covariate matrix for Omega
		
	W_pred_mat = array(-99999, c(n_z, p, (T-1), S_A))
	
	for(i in 1:n_z) W_pred_mat[i, , , ] = pred_mat_3D[ , Bt_all[i, -10], ] - pred_mat_3D[ , -1, ]
	
	W_pred_mat = aperm(W_pred_mat, c(2, 3, 4, 1))

	dim(W_pred_mat) = c( p,  (T - 1)*S_A*n_z )  

	######################################################################
	
	y_temp = y_full[ , available_locations_unique]

	pos_avail = (y_temp != -99)
	
	pos_unavail = (y_temp == -99 )
	
	y = y_temp[pos_avail]	
		
	n = length(y)

	##########################################################################
	
	## Construction of covariate matrix for all possible combinations of Z
	
	A_st = array(-99999, c(n_z, p, T, S_A)) 

	for(i in 1:n_z) A_st[i, , , ] = pred_mat_3D[ , Bt_all[i, ], ]

	A_st = aperm(A_st, c(2, 3, 4, 1))

	dim(A_st) = c( p,  T*S_A*n_z )

	pos_avail_3D = rep(pos_avail, n_z)
	
	dim(pos_avail_3D) = c(T, S_A, n_z)

	## Construction of y^(2) to model all possible combination of Z
		
	y_temp_3D = array(0, c(T, S_A, n_z))

	for(i in 1:n_z) y_temp_3D[ , , i] = y_temp

	aa = lapply(1:S_A, function(i) which(y_temp[ , i] != -99))
	
	group_y = rep(1:S_A, sapply(aa, length))
	
	group_z = rep(1:S_A, each = (T - 1))
	
	cell = 1:S_A

	n_omega = S_A*(T-1)
	
	####################################################################################################################
	
	## X matrix for available and unavailable y^(2) 
	
	pred_mat = pred_mat_3D
	
	dim(pred_mat) = c(p, T*S_A)
	
	pred_mat_UA = pred_mat_full_3D[ , , unavailable_locations_unique]
	
	############################################################################################################################

	## prior distributions and initial values
	
	beta = c( mean(y_full[pos_avail_full]), rep(0, (p-1) ) )
				
	sigma2 = var(y_full[pos_avail_full]) 
	
	# Horseshoe related stuff:

	tau2_beta    = 1
	zai      = 1
    
	lambda2 = rep(1, (p - 1))
	nu = rep(1, (p - 1))
    
	B = c(100 , tau2_beta*lambda2)

	## prior of sigma2 is IG(a0,b0)

	a0 = 2.000001
	b0 = 1.000001

	nit = 10000   ## number of initial runs that you want to discard
	nmc = 50000     ## number of runs that you want to do
	nthin = 10   ## you want t thin the sample at what interval to avoid correlation ?

	## during the mcmc, you need to store the simulated samples of mu and sigma, so define storage variables
	
	y_store = array(0, dim = c(T, S_A, nmc/nthin) )
	
	z_store = array(0, dim = c(T,  S_A, nmc/nthin) )
	
	omega_store = array(0, dim = c(n_omega, nmc/nthin) )
	
	beta_store = matrix(0, p, nmc/nthin)
		
	sigma2_store = c()
			
	nlog_likelihood_store = c()

	constrained_likelihood_store = matrix(0,  nmc/nthin, n)
				
	z_pred_store = array(-99, dim = c(T, S_UA, nmc/nthin) )
	
	y_pred_store = array(0, dim = c(T, S_UA, nmc/nthin) )

print(date())

for (iter in 1:(nit+nmc))
	{
	
	###################################################################################

	## updating z from multinomial distribution
	
	## data part
	
		A_st_beta = colSums(A_st*beta)

		dim(A_st_beta) = c( T, S_A, n_z )

		temp_1 =  - 0.5* (y_temp_3D[pos_avail_3D] - A_st_beta[pos_avail_3D])^2/sigma2

		dim(temp_1) = c(n, n_z)

		data_part = rowsum(temp_1, group_y)

	## prior part 
	
		temp_2 = colSums(W_pred_mat*beta)

		dim(temp_2) = c( (T - 1), S_A, n_z )

		temp_2 = aperm(temp_2, c(3, 1, 2))

		temp_a = c(z_all[ , -1])*c(pnorm(temp_2, log.p = TRUE))

		temp_b = c(1 - z_all[ , -1])*c(pnorm(temp_2, lower.tail = FALSE, log.p = TRUE))

		dim(temp_a) = c(n_z, (T - 1)*S_A)

		dim(temp_b) = c(n_z, (T - 1)*S_A)

		temp_d = cbind(temp_a, temp_b)

		z_part = rowsum(t(temp_d), c(group_z, group_z))

		prob_mat = exp(z_part + data_part)
				
		pos_z = rmult_new_multiple( prob_mat/rowSums(prob_mat) )
		
		z_matrix = t( z_all[pos_z, ] )

		Bt_matrix = t( Bt_all[pos_z, ] )

###################################################################################
		
	## Updating omega from truncated normal distribution

		upper_bound_W = matrix(Inf, (T - 1), S_A)

		lower_bound_W = matrix(-Inf, (T - 1), S_A)
		
		pos_z_1 = (z_matrix[-1, ] == 1)
		
		pos_z_0 = (z_matrix[-1, ] == 0)

		upper_bound_W[pos_z_0] = 0
		
		lower_bound_W[pos_z_1] = 0

		temp_pos_W =  rep( ( (T-1)*S_A*(pos_z - 1) + (T - 1)*(cell - 1) + 1 ), each = (T - 1) ) + (0:(T - 2)) 

		temp_W_pred_mat = W_pred_mat[, temp_pos_W]
		
		mu_omega = colSums(temp_W_pred_mat*beta) 
		
		omega = rtrunc_norm(number = n_omega, a = c(lower_bound_W), b = c(upper_bound_W), mu = mu_omega, sigma2 = 1)
		
	##############################################################################################################################

	## Update Beta from multivariate normal distribution
	
		temp_pos_beta = rep( ( T*S_A*(pos_z - 1) + T*(cell - 1) + 1 ), each = T ) + (0:(T - 1))
		
		temp_X =  A_st[ , temp_pos_beta]
	
		X_mat = t( temp_X[ , c(pos_avail)] )
			
		v1_inv = crossprod( X_mat )/ sigma2
		
		v2_inv = crossprod( t(temp_W_pred_mat) )
		
		v3_inv = diag(1/B)

		m1_by_v1 = crossprod( X_mat, y) / sigma2
		
		m2_by_v2 = crossprod( t(temp_W_pred_mat), omega )	## since sigma2 = 1
		
		beta_post_var = inv_and_logdet.sym(v1_inv + v2_inv + v3_inv)[[2]]

		beta_post_mean = crossprod(beta_post_var, m1_by_v1 + m2_by_v2)
		
		beta = c( rmvnorm(1, beta_post_mean, beta_post_var) )

	############################################################################################################	

	## Horseshoe related parameters update

	##  while simulating beta parameters, their prior mean is 0 and prior dispersion is diag(temp_b)

	## after you simulated the beta, simulate horseshoe related parameres
        
		tau2_beta = 1/rgamma( 1,shape = p/2, rate = 1/zai + 0.5*sum((beta[-1]^2.0)/lambda2) )
        
		zai = 1/rgamma(1, shape = 1, rate = 1 + 1/tau2_beta)
        
		lambda2 = 1/rexp((p - 1), rate = 1/nu + 0.5*(beta[-1]^2.0)/tau2_beta)
        
		nu = 1/rexp( (p - 1 ), rate = 1 + 1/lambda2)
    
		B[-1] = tau2_beta*lambda2
		
	############################################################################################################
	
	## mean of the data

		mu = colSums(t(X_mat)*beta)

	############################################################################################################

	## update sigma2  from Inverse Gamma distribution

		S2 = (y - mu)^2
   
		sigma2_post_shape = a0 + n/2
		
		sigma2_post_rate = b0 + 0.5*sum(S2)

		sigma2 = 1/rgamma(1, shape = sigma2_post_shape, rate = sigma2_post_rate)

	if (iter > nit & (iter - nit)%%nthin==0)
		{
					
			y_temp[pos_unavail] = colSums(temp_X[ , c(pos_unavail)]*beta) + sqrt(sigma2)* rnorm(T*S_A - n)
			
			y_store[ , , (iter - nit)/nthin] = y_temp             
			
			z_store[ , , (iter - nit)/nthin] = z_matrix
			
			omega_store[ , (iter - nit)/nthin] = omega
			
			beta_store[ , (iter - nit)/nthin] = beta
			
			sigma2_store[(iter - nit)/nthin] = sigma2 
			
			nlog_likelihood_store[(iter - nit)/nthin] = 0.5 * n * log(sigma2) + 0.5 * sum(S2)/sigma2
			
			constrained_likelihood_store[(iter - nit)/nthin, ] = exp( - (0.5 * log(sigma2) + 0.5 * S2/sigma2) )
						
		}
		
		if(iter == 1000) print(date())

	}

print(date())

	###################################################################################################################
	
	## predication of missing values

	mu_pred_store = array(-99999, c(T, S_UA, nmc/nthin))
	
	mu_pred_store[1, , ] = t( pred_mat_UA[ , 1,  ]) %*% beta_store
	
	y_pred_store[ 1 , , ] = mu_pred_store[1, , ] + t( sqrt(sigma2_store)* matrix(rnorm(S_UA*nmc/nthin), nmc/nthin) )

	z_pred_store[1, , ] = 1
	
	for (t in  2:T)	
		{
			mu_temp =  t( pred_mat_UA[ , t,  ] )%*% beta_store
			
			z_pred_store[t, , ] = ceiling( c(pnorm( t(t(mu_pred_store[(t-1), , ] - mu_temp )) ) ) - matrix(runif(S_UA*nmc/nthin), S_UA) ) 

			mu_pred_store[t, , ] = z_pred_store[t, , ]*mu_temp + (1 - z_pred_store[t, , ])*mu_pred_store[(t-1), , ]
			
			y_pred_store[t, , ] = mu_pred_store[t, , ] + t( sqrt(sigma2_store)* matrix(rnorm(S_UA*nmc/nthin), nmc/nthin) )
		}

	##########################################################################################
	
	# Log likelihood and BPIC calculation
	
	LL = median(-nlog_likelihood_store) ## Log likelihood

	mean_beta = c( apply(beta_store, 1, mean) )
	
	mean_sigma2 = mean(sigma2_store)
	
	mean_z = apply(z_store, c(1, 2), mean)
		
	mean_mu = matrix(-9999, T, S_A)
	
	mean_mu[1, ] = colSums(pred_mat_3D[ , 1, ]*mean_beta)
	
	for (t in 2:T) mean_mu[t, ] = mean_z[t, ]* colSums(pred_mat_3D[ , t, ]*mean_beta) + (1 - mean_z[t, ])* mean_mu[(t - 1), ]
	
	mean_S2 = (y - mean_mu[pos_avail])^2
		
	mean_nlog_likelihood = 0.5 * n * log(mean_sigma2) + 0.5 * sum(mean_S2)/mean_sigma2

	D_bar = 2*mean(nlog_likelihood_store)
	
	D_theta_bar_mean = 2*mean_nlog_likelihood

	BPIC = 3*D_bar - 2*D_theta_bar_mean

	###########################################################################################################################

	basal_area_full_temp = abind(y_store, y_pred_store, along = 2)

	basal_area_MCMC = array(0, dim = dim(basal_area_full_temp))

	basal_area_MCMC[ ,  c(available_locations_unique, unavailable_locations_unique), ]   = basal_area_full_temp

## the given unit of basal area is ft^2/acre. To transform it into m^2/hectare it is needed to mulitply the basal area by 0.229568411	

	basal_area_MCMC_exp = exp(basal_area_MCMC)*0.229568411	## transforming  ft^2/acre  to m^2/hectare 

	rm(basal_area_full_temp, basal_area_MCMC, y_pred_store)

	#########################################################################################################

	constrained_likelihood = colMedians(constrained_likelihood_store)
	
	#########################################################################################################
	#########################################################################################################

	## saving necessary outputs

	outputs = list(constrained_likelihood = round(constrained_likelihood,3), LL_BPIC = round(c(LL, BPIC), 3),  basal_area_MCMC_exp = basal_area_MCMC_exp) 
	
	save(outputs, file = "Inputs_and_Outputs/constrained_y2_model_output.Rdata")