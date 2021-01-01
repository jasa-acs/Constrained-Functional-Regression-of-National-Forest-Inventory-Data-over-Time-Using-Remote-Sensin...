############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_FIA_data_analysis.R that contains all necessary functions.

## The output of TC_best_output_merge_basal_area.R that is "best_TC_model_outputs_and_basal_area_data.Rdata" where BEST should be replaced by the best TC model number (IV in our case) from the five candidate model numbers from Appendix A.1 in the supplementary materials - I, II, IIIA, IIIB and IV.

## best_TC_model_outputs_and_basal_area_data.Rdata contains a list reconstructed_TC_BALIVE consisting of five components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the complete TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 
## (iii) The fifth component is the live tree basal area in the form of a 10000x120 matrix where rows and colums correspond to pixels and time points (months), respectively. The coccurance of -99 within this matrix indicates unobserved data at those space-time combinations.


## one output of Basal_area_cv_positions_selction.R that is random_positions_for_CV_for_y1_model.Rdata
## random_positions_for_CV_for_y1_model.Rdata contains a list test_pos_list of 36 components where each component contains the postions of both zero, non-zero live tree basal area data those will be considered as test data in the correspondig holdout cross-validation method.


## one output of Model_for_y1_cv.R that is CV_output_for_y1_model.Rdata
## (i) CV_output_for_y1_model.Rdata that contains the MCMC outputs of probability of forest at each test data in the form of a 540x(nmc/nthin) matrix named prob_CV_storage where rows and columns correspond to the test data and MCMC outputs of corresponding test data after burn in and thinning, respectively.

############################################################################################
########################		Output     #################################################
############################################################################################

## the output is "CV_output_basal_area_model.csv" that is a 4x2 matrix where the second column represents the corresponding comparison statistic (mean of absolute bias, uncertainty and empirical coverage) as mentioned in the first column.

################################################################################################################
#####################     Cross Validation     #################################################################
################################################################################################################

set.seed(77777)

library(mvtnorm)
library(matrixStats) # for colMedians, rowCumsums
library(coda)	## for HPD interval

source("FIA_Codes/Functions_for_FIA_data_analysis.R")

load("Inputs_and_Outputs/best_TC_model_outputs_and_basal_area_data.Rdata")	## loading the data 

load("Inputs_and_Outputs/random_positions_for_CV_for_y1_model.Rdata")

	locations_all = reconstructed_TC_BALIVE$locations	## location of the cells 

	S = nrow(locations_all)	## Total number of grid-cells in response

	T = 10	## Total number of time points (years here) in response

	n_TC = 3	## Total number TC features

	n_total = S*T	## Total number of data points in response

	l_test = 15	## number of samples that are considered as tested data

	n_CV = 36		## total number of replication of holdout method 

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

	for (i in 1:length(center)) pred_mat_full[i , ] = (pred_mat_full[i , ] - center[i])/scale[i]

	pred_mat_full = rbind(1, pred_mat_full)
	
	p = nrow(pred_mat_full)

	pred_mat_full_3D = pred_mat_full
		
	dim(pred_mat_full_3D) = c( p, T, S)

	##############################################################################################

	## response 

	## here -99 represents data is not available
	
	y_full = matrix(-99,  S,  T)

	pos_0 = which(reconstructed_TC_BALIVE$BALIVE == 0, arr.ind = TRUE)

	pos_0 = pos_0[which(pos_0[ , 2] >  0 & pos_0[ , 2] < 121), ]

	pos_1 = which(reconstructed_TC_BALIVE$BALIVE > 0, arr.ind = TRUE)

	pos_1 = pos_1[which(pos_1[,2]> 0 & pos_1[,2]<121),]

	for(s in 1:(dim(pos_1)[1])) y_full [pos_1[s,1] , ceiling(pos_1[s,2]/12)] = log(reconstructed_TC_BALIVE$BALIVE[pos_1])[s]

	for(s in 1:(dim(pos_0)[1])) y_full[pos_0[s, 1], ceiling(pos_0[s, 2]/12)] = 0

	available_locations_all = list()

	for (t in 1:T) available_locations_all[[t]] =   sort(which(y_full[ , t] != -99)) 

	available_locations_unique_temp = sort(unique(unlist(available_locations_all)) ) ## this stores all the cell locations, where any nonzero data was available during 10 years

	y_full = t(y_full)

	y_all = y_full[ , available_locations_unique_temp]

	pos_avail_temp = which(y_all != -99)

cross_validation = list()

for (cv in 1:n_CV)
		
	{
	
	set.seed(77777)

	pos = test_pos_list[[cv]]
	
	y_temp_all = y_all

	y_test = y_temp_all[pos_avail_temp[pos]]

	y_temp_all[pos_avail_temp[pos]] = -99

	pos_0_temp = which(y_temp_all == 0)

	y_temp_all[pos_0_temp] = -99

	y = y_temp_all[y_temp_all != -99]
		
	n = length(y)

	y_temp_all[pos_avail_temp[pos]] = y_test

	temp_available_locations = list()

	for (t in 1:T) temp_available_locations[[t]] = sort(which(y_temp_all[ t , ] != -99)) 
	
	available_locations = sort(unique(unlist(temp_available_locations)) )

	available_locations_unique = available_locations_unique_temp[available_locations]

	S_A =  length(available_locations_unique)
	
	############################################################################################################################
	
	## Some calculations to update z from all possible combination of z using multinomial distribution

	z_all = as.matrix( unname( cbind( 1, expand.grid(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1)) ) ) )
	
	n_z = nrow(z_all)
	
	Bt_all = t( apply(z_all, 1, Bt_matrix_create) )
	
	pred_mat_3D = pred_mat_full_3D[ , , available_locations_unique]
	
	A_st = array(-99999, c(p, T, S_A, n_z)) 

	for(i in 1:n_z) A_st[ , , , i] = pred_mat_3D[ , Bt_all[i, ], ]

	dim(A_st) = c( p,  T*S_A*n_z )

	y_temp = y_temp_all[ , available_locations]
		
	pos_test = rep(-999, l_test)

	pos_0_test = which(y_test == 0)

if(length(pos_0_test)>0 ) { pos_test[pos_0_test] = which(y_temp == 0); pos_test[- pos_0_test] = match(y_test[-pos_0_test], y_temp) } else {pos_test = match(y_test, y_temp) }

	y_temp[pos_test] = -99
	
	pos_avail = (y_temp != -99)
	
	pos_avail_3D = rep(pos_avail, n_z)
	
	dim(pos_avail_3D) = c(T, S_A, n_z)

	y_temp_3D = array(0, c(T, S_A, n_z))

	for(i in 1:n_z) y_temp_3D[ , , i] = y_temp
	
	## Construction of X matrix for omega
	
	temp_pos = apply(pos_avail, 2, sum)
	
	group_y = rep(1:S_A, temp_pos)

	temp_pos_0 = which(temp_pos == 0)
	
	n_temp_pos_0 = length(temp_pos_0)
		
	W_pred_mat = array(-99999, c(p, (T-1), (S_A - n_temp_pos_0), n_z))

if(n_temp_pos_0 == 0){ for(i in 1:n_z) W_pred_mat[ , , , i] = pred_mat_3D[ , Bt_all[i, -10],  ] - pred_mat_3D[ , -1,  ]} else {for(i in 1:n_z) W_pred_mat[ , , , i] = pred_mat_3D[ , Bt_all[i, -10], -c(temp_pos_0) ] - pred_mat_3D[ , -1, -c(temp_pos_0) ] }
	
	dim(W_pred_mat) = c( p,  (T - 1)*(S_A - n_temp_pos_0)*n_z )  

	group_z = rep(1:(S_A - n_temp_pos_0), each = (T - 1))
	
	cell_X = setdiff(1:S_A, temp_pos_0)

	cell_W = 1:(S_A - n_temp_pos_0)

n_omega = if( n_temp_pos_0 == 0) { length(y_temp[-1, ])} else { length(y_temp[-1, -c(temp_pos_0) ]) }
	
	
if(n_temp_pos_0 > 0)
	{

	aa = y_temp_all[ , available_locations]

	aa[-pos_test] = -99

	pos_test_1 = which(aa[ , -c(temp_pos_0)] != -99)

	pos_test_2 = which(aa[ , c(temp_pos_0)] != -99)

	xx = yy = aa

	xx[ , c(temp_pos_0)] = -99

	yy[ , -c(temp_pos_0)] = -99

	pos_match_1 = match(which(xx != -99), pos_test)

	pos_match_2 = match(which(yy != -99), pos_test)

	pos_match = match(1:l_test, c(pos_match_1, pos_match_2))

	}
	####################################################################################################################
	
	## z and omega is updated from the second time point and for first year they are 1 (one)
	
	pred_mat = pred_mat_3D
	
	dim(pred_mat) = c(p, T*S_A)
	
	############################################################################################################################

	## initial values
	
	beta = c( mean(y), rep(0, (p-1) ) )
				
	sigma2 = var(y) 
	
	# Horseshoe related stuff:

	tau2_beta = 1
	
	zai = 1
    
	lambda2 = rep(1, (p - 1))
	
	nu = rep(1, (p - 1))
    
	B = c(100 , tau2_beta*lambda2)

	## prior of sigma2 is IG(a0,b0)

	a0 = 2.000001
	b0 = 1.000001

#########################################################################################################

	nit = 10000   	## number of initial runs that you want to discard
	nmc = 50000     ## number of runs that you want to do
	nthin = 10   	## you want t thin the sample at what interval to avoid correlation ?

	## during the mcmc, you need to store the simulated test data, so define storage variables
	
	y_pred_store = matrix(0, l_test, nmc/nthin)
			
for (iter in 1:(nit+nmc))
	{
	
	#########################################################################################################
	
	## updating z from multinomial distribution
	
	## data part
	
		A_st_beta = colSums(A_st*beta)

		dim(A_st_beta) = c( T, S_A, n_z )

		temp_1 =  - 0.5* (y_temp_3D[pos_avail_3D] - A_st_beta[pos_avail_3D])^2/sigma2

		dim(temp_1) = c(n, n_z)

		data_part = rowsum(temp_1, group_y)

	## prior part 
	
		temp_2 = colSums(W_pred_mat*beta)

		dim(temp_2) = c( (T - 1), (S_A - n_temp_pos_0), n_z )

		temp_2 = aperm(temp_2, c(3, 1, 2))

		temp_a = c(z_all[ , -1])*c(pnorm(temp_2, log.p = TRUE))

		temp_b = c(1 - z_all[ , -1])*c(pnorm(temp_2, lower.tail = FALSE, log.p = TRUE))

		dim(temp_a) = c(n_z, (T - 1)*(S_A - n_temp_pos_0))

		dim(temp_b) = c(n_z, (T - 1)*(S_A - n_temp_pos_0) )

		temp_d = cbind(temp_a, temp_b)

		z_part = rowsum(t(temp_d), c(group_z, group_z))

		prob_mat = exp(z_part + data_part)
				
		pos_z = rmult_new_multiple( prob_mat/rowSums(prob_mat) )
		
		z_matrix = t( z_all[pos_z, ] )

		Bt_matrix = t( Bt_all[pos_z, ] )
	
	#########################################################################################################
	
	## Updating omega from truncated normal distribution

		upper_bound_W = matrix(Inf, (T - 1), (S_A - n_temp_pos_0) )

		lower_bound_W = matrix(-Inf, (T - 1), (S_A - n_temp_pos_0))
		
		pos_z_1 = (z_matrix[-1, ] == 1)
		
		pos_z_0 = (z_matrix[-1, ] == 0)

		upper_bound_W[pos_z_0] = 0
		
		lower_bound_W[pos_z_1] = 0

		temp_pos_W = ( rep( ( (T-1)*(S_A - n_temp_pos_0)*(pos_z - 1) + (T - 1)*(cell_W - 1) + 1 ), each = (T - 1) ) + (0:(T - 2)) )

		temp_W_pred_mat = W_pred_mat[, temp_pos_W]
		
		mu_omega = colSums(temp_W_pred_mat*beta) 
		
		omega = rtrunc_norm(number= n_omega, a = c(lower_bound_W), b = c(upper_bound_W), mu = mu_omega, sigma2 = 1)
		
	##############################################################################################################################

	## Beta update from multivariate normal distribution
	
		temp_pos_beta = rep( ( T*S_A*(pos_z - 1) + T*(cell_X - 1) + 1 ), each = T ) + (0:(T - 1))
		
		temp_X =  A_st[ , temp_pos_beta]

	X_mat = if(n_temp_pos_0 == 0) {t( temp_X[ , c(pos_avail)] )} else {t( temp_X[ , c(pos_avail[ , - c(temp_pos_0) ])] )}
				
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

	###########################################################################################################
	
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
						
if(n_temp_pos_0 == 0)	{
			
			y_pred_store[ , (iter - nit)/nthin] = colSums(temp_X[ , c(pos_test)]*beta) + sqrt(sigma2)* rnorm(l_test)

			}            

if(n_temp_pos_0 > 0)	{
			
			temp_y_pred_1 = colSums(temp_X[ , c(pos_test_1)]*beta) + sqrt(sigma2)* rnorm(length(c(pos_test_1)))             
			
			temp_mu_pred = matrix(-99999, T, n_temp_pos_0)
	
			temp_mu_pred[1, ] = t( pred_mat_3D[ , 1, temp_pos_0 ]) %*% beta
	
			temp_y_pred = matrix(-99999, T, n_temp_pos_0)
	
			temp_y_pred[ 1 , ] = c(temp_mu_pred[1, ]) + sqrt(sigma2)* rnorm(n_temp_pos_0)

			temp_z_pred = matrix(-99999, T, n_temp_pos_0)
			
			temp_z_pred[1, ] = 1
	
	for (t in  2:T)	
		{
			mu_temp =  t( pred_mat_3D[ , t, temp_pos_0 ] )%*% beta
			
			temp_z_pred[t, ] = ceiling( pnorm( temp_mu_pred[(t-1), ] - mu_temp ) - runif(n_temp_pos_0) ) 

			temp_mu_pred[t, ] = temp_z_pred[t, ]*mu_temp + (1 - temp_z_pred[t, ])*temp_mu_pred[(t-1), ]
			
			temp_y_pred[t, ] = temp_mu_pred[t, ] + sqrt(sigma2)* rnorm(n_temp_pos_0)
		}

			y_pred_store[ , (iter - nit)/nthin] = c( temp_y_pred_1, temp_y_pred[pos_test_2] )[pos_match]            

			}			
		}
				
	}

	###################################################################################################################

	cross_validation[[cv]] = cbind(y_test, y_pred_store)

}


	test_pred = do.call(rbind, cross_validation)

	pos_0_test = which(test_pred[, 1] == 0)

## the given unit of basal area is ft^2/acre. To transform it into m^2/hectare it is needed to mulitply the basal area by 0.22956841. But in log scale multiplication becomes addition.
	
	test_pred = test_pred + log(0.229568411)	## transforming  ft^2/acre  to m^2/hectare 
	
	test_pred[pos_0_test, 1] = 0
	
	y_test_all = exp(test_pred[ , 1])
	
	y_test_all[pos_0_test] = 0

	load("Inputs_and_Outputs/CV_output_for_y1_model.Rdata")
	
	prob_matrix = ceiling(prob_CV_storage - matrix(runif(length(prob_CV_storage)), nrow(prob_CV_storage), ncol(prob_CV_storage)))

	y_test_MCMC = exp(test_pred[, -1]) * prob_matrix

	y_test_pred = apply(y_test_MCMC, 1, median)

	y_CI_HPD = apply(y_test_MCMC, 1, function(x) HPDinterval(mcmc(x), prob = 0.90)[1,] )

	outputs = cbind( y_test_all, y_test_pred, y_CI_HPD[1, ], y_CI_HPD[2, ])

	colnames(outputs) = c("y_test", "y_pred",  "HPD LI", "HPD UI")

	abs_bias = mean(abs(outputs[ , 1] - outputs[ , 2]))

	uncertainty = mean(outputs[ , 4] - outputs[ , 3])

	empirical_coverage = sum(as.integer(outputs[ , 1] >= outputs[ , 3] & outputs[ , 1] <= outputs[ , 4]))/nrow(outputs)

	CV_outputs = round( c(abs_bias, uncertainty, empirical_coverage), 3)

	a = c("Absolute Bias", "Uncertainty", "Empirical Coverage")
	
	CV_outputs = data.frame(cbind(a, CV_outputs))

	colnames(CV_outputs) = c("Comparison Statistic", "Basal Area Model")
	
	write.csv(CV_outputs, file = "Inputs_and_Outputs/CV_output_basal_area_model.csv", row.names = F)
