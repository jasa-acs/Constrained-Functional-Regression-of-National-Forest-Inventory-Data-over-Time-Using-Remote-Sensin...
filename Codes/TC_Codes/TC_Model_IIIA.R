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

## the output is "TC_tc_Model_IIIA_outputs.Rdata" where tc should be replaced by appropriate TC number from 1, 2 and 3 that is a list fill_up_missing_data containing four components:

## (i) the first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.

## (ii) The second component contains the complete spatiotemporal map of corresponding TC data, after filling in the missing observations, in form of a 10000x120 matrix where first and second dimensions correspond to pixels and time points (months), respectively. 

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
library(ltsa)
library(coda)

source("TC_Codes/Functions_for_TC_analysis.R")

load("Inputs_and_Outputs/coarse_selected_data_with_missing_values.Rdata") ## loading output of Aggregation.R that contains a list coarse_creation

TC = as.integer(readline(" input 1 for TC1, 2 for TC2, 3 for TC3 "))	## choose the TC covariate to model

## data representing the corresponding TC measurement chosen in the previous line where rows and columns represents grid-cells and months respectively

data = coarse_creation$coarse_data[ , , TC]

## In data 0 represents missing value 

##########################################################################################

## To translate and rescale each TC features to the same range as that of TC1 

data_TC1 = coarse_creation$coarse_data[ , , 1]	## data_TC1 is the TC1 measurement

pos_0_cell = apply(data, 1, function(x) which(x == 0))	## list of size S where each list contains all missing months at the corresponding cell

pos_0_cell_length = sapply(pos_0_cell, length)	## vector of size S containing the number of missing months at the corresponding cell

pos_0_month = apply(data, 2, function(x) which(x == 0))	 ## list of size T where each list contains all missing cells at the corresponding month																				

pos_0_month_length = sapply(pos_0_month, length)	## vector of size T containing the number of missing cells at the corresponding month

pos_0 = which(data == 0) 	## vector of the positions of all missing data 

## following available data are translated and rescaled to the same range as that of TC1, and then log transformation is performed 

data[-pos_0] = log( data[-pos_0]* diff(range(data_TC1[-pos_0]))/diff(range(data[-pos_0])) - 
			   min( data[-pos_0]* diff(range(data_TC1[-pos_0]))/diff(range(data[-pos_0]))) + min(data_TC1[-pos_0]))

############################################################################################

location = coarse_creation$coarse_location	## location of the cells 

S = dim(location)[1]	## Total number of grid-cells in the data

T = dim(data)[2]	 ## Total number of months in the data

n_total = S*T	## Total number of data points

n_non_missing = n_total - length(pos_0)	## Total number of non-missing data points

month_missing = which(pos_0_month_length == S)	## Months with 100% missing data

month_all = which(pos_0_month_length == 0)	## Months with no missing data

## Binary vector where 0 represents 100% missing data at the corresponding months and 1 otherwise

z_psi = rep(1, T)

z_psi[month_missing] = 0

covariate = t(sapply(1:T, function(t)  c(1, cos((1:2)*pi*t/6), sin((1:2)*pi*t/6))))	## Fourier series with first two harmonics at any cell

p_beta = dim(covariate)[2]

X_train = covariate[rep(1:T, (S- pos_0_month_length)),]	## Construction of full covariate matrix for available TC feature

XTX = crossprod(X_train)

#############################################################################################
	
## Construction of neighborhood Matrix

x_loc = sort(unique(location[,1]))

x_dist= x_loc[2] - x_loc[1]

y_loc = sort(unique(location[,2]))

y_dist= y_loc[2] - y_loc[1]

W = sapply(1:S,f_W)

W_mat = matrix(0, S, S)

for(i in 1:S)
{
   W_mat[i, W[[i]]] = 1
}

W_plus = sapply(W, length)

################################################################################################

## set priors and initial values

## beta ~ MVN(beta0,c0*I)

beta0 = matrix(0, p_beta, 1)
c0 = 1000

## prior for sigma2 is IG(a10,b10)

a10 = 2.000001
b10 = 1.000001

## prior for tau2 is IG(a20,b20)

a20 = 2.000001
b20 = 1.000001

## prior for kappa2 is IG(a30,b30)

a30 = 2.000001
b30 = 1.000001

## prior for eta is S(0, nu2)I_{-1,+1}

nu2 = 1

## initial values 

phi = rep(0, S)

tau2 = 1

psi = rep(0, T)

psi_0 = 0

eta = 0

kappa2 = 0.01

#####################################################################################

## Initial values from least-square method

beta = crossprod(t(inv_and_logdet.sym(XTX)[[2]]), crossprod(X_train, data[-pos_0]))

e = data[-pos_0] - crossprod(t(X_train), beta)

sigma2 = var(c(e))

#############################################################################################################

## how many times you want to run the mcmc

nit = 10000    ## number of initial runs that you want to discard
nmc = 50000    ## number of runs that you want to do
nthin = 50	## you want t thin the sample at what interval to avoid correlation ?

## during the mcmc, you need to store the simulated parameters, so define storage variables

beta_store = matrix(0, nmc/nthin , p_beta)

sigma2_store = c()

phi_store = matrix(0, nmc/nthin, S)

tau2_store = matrix(0, nmc/nthin, 1)

psi_store = matrix(0, nmc/nthin, T)

psi_0_store = matrix(0, nmc/nthin, 1)

eta_store = matrix(0, nmc/nthin, 1)

kappa2_store = matrix(0, nmc/nthin, 1)

nlog_likelihood_store = c()	## negative log-likelihood storage

pred_store = matrix(0, nmc/nthin, length(pos_0))	## predicting the missing values

## start the mcmc loop

print(date())

for (iter in 1:(nit+nmc))
	{
	
	#############################################################################################################
	
	## posterior distribution of beta is Multivariate Normal

		y_train = t(t(data - phi) - psi)
	
		beta_post_dispersion = inv_and_logdet.sym((XTX/sigma2)+ diag(p_beta)/c0)[[2]]     ## do the calculation in your notebook and then write the formula here

		beta_post_mean =  crossprod(t(beta_post_dispersion),((crossprod(X_train, y_train[-pos_0])/sigma2)+(beta0/c0)))         ## complete the formula here
    
	## then simulate sample of beta
    
		beta = c(rmvnorm (1, beta_post_mean, beta_post_dispersion))

	#############################################################################################################
	
		X_beta = c(crossprod(t(covariate), beta))
	
		y_minus_X_beta = t(t(data) - X_beta)

	#############################################################################################################
	
	## posterior distribution of sigma2 is Inverse Gamma
	
		sum_S2 = sum((t(t(y_minus_X_beta - phi) - psi)[-pos_0])^2)
	
		sigma2_post_shape = a10 + n_non_missing/2
		
		sigma2_post_rate = b10 + sum_S2/2
    
	## then simulate sample of sigma2
    
		sigma2 = 1/rgamma(1, shape = sigma2_post_shape, rate = sigma2_post_rate)

	#############################################################################################################
	###############################################################################################################
	
	## simulation of all necessary parameters related to spatial random effect
		
	## posterior distribution of phi is normal		
 
		for (j in 1:S)
			{
				phi_var = 1/((T-pos_0_cell_length[j])/sigma2 + W_plus[j]/tau2)
				
				phi_mu = phi_var * ((T - pos_0_cell_length[j]) * mean((y_minus_X_beta[j,] - psi)[-pos_0_cell[[j]]])/sigma2 + sum(phi[W[[j]]])/tau2)
				
    ## then simulate sample of phi
	
				phi[j] = rnorm(1, phi_mu, sqrt(phi_var))
			}

	## standardizing phi	
	
		phi = phi - mean(phi)

	## posterior distribution of tau2 is Inverse Gamma
 
		tau2_post_shape = a20 + (S-1)/2  	
		
		tau2_post_rate = b20 + sum((phi[rep(1:S, W_plus)] - phi[unlist(W)])^2)/4

	## then simulate sample of tau2
	
		tau2 = 1/rgamma(1, shape = tau2_post_shape, rate = tau2_post_rate)

	###############################################################################################################
	###############################################################################################################
	
	## simulation of all necessary parameters related to temporal random effect

	## simulate sample of psi_0 whose posterior distribution is normal

		psi_0 = rnorm(1, eta* psi[1], sqrt(kappa2))

		psi_temp = c(psi_0, psi[1:(T-1)])
	
	## posterior distribution of kappa2 given others is Inverse Gamma

		kappa2_post_shape = a30 + (T+1)/2 
		
		kappa2_post_rate = b30 + sum((psi - eta*psi_temp)^2)/2 + (1 - eta^2)* psi_0^2/ 2 
	
	## then simulate sample of kappa2
	
		kappa2 = 1/rgamma(1, shape = kappa2_post_shape, rate = kappa2_post_rate)
	
	## simulate eta via slice sampling
	
		lim_eta = sqrt(1 - (runif(1, 0, sqrt(1-eta^2)))^2)
	
		eta_post_var = 1/ (sum((psi[1:(T-1)])^2)/kappa2 + 1/nu2)
		
		eta_post_mean = eta_post_var * sum(psi*psi_temp)/ kappa2
			
		eta = rtrunc_norm(number=1, a = -lim_eta, b = lim_eta, mu = eta_post_mean, sigma2 = eta_post_var) 

	## posterior distribution of psi is Multivariate Normal
	
		psi_dispersion =  create_AR_cov_mat(kappa2 = kappa2, eta = eta, total_months = T)
		
		psi_mean_temp = matrix(0, T, 1)
		
		for (t in month_all) psi_mean_temp[t,] = mean( y_minus_X_beta[,t] )
		
		for (t in setdiff(1:T, union(month_all, month_missing))) psi_mean_temp[t,] = mean( y_minus_X_beta[-pos_0_month[[t]],t] - phi[- pos_0_month[[t]]] )

		psi_post_dispersion = inv_and_logdet.sym(diag(z_psi*(S - pos_0_month_length)/sigma2) 
							+ inv_and_logdet.sym(psi_dispersion)[[2]] )[[2]]		

		psi_post_mean = crossprod(t(psi_post_dispersion), z_psi*(S - pos_0_month_length) * psi_mean_temp/sigma2)

	## then simulate sample of psi for observations

		psi = c(rmvnorm(1, psi_post_mean, psi_post_dispersion) )

	###############################################################################################################
	###############################################################################################################
		
	## Now store the samples in the store vector you already created
   
	if (iter > nit & (iter - nit)%%nthin==0)
    {
		beta_store[(iter - nit)/nthin, ] =  beta

		sigma2_store[(iter - nit)/nthin] = sigma2

		tau2_store[(iter - nit)/nthin, ] = tau2

		phi_store[(iter - nit)/nthin, ] = phi 

		psi_store[(iter - nit)/nthin, ] = psi

		psi_0_store[(iter - nit)/nthin, ] = psi_0

		eta_store[(iter - nit)/nthin, ] = eta

		kappa2_store[(iter - nit)/nthin, ] = kappa2

		res_temp = t(t(y_minus_X_beta - phi) - psi)
		
		nlog_likelihood_store[(iter - nit)/nthin] = n_non_missing*log(sigma2)/2 + sum((res_temp[-pos_0])^2)/(2*sigma2)		

		pred_store[(iter - nit)/nthin, ] =  t(t(- y_minus_X_beta  + phi) + psi) [pos_0] + rnorm(length(pos_0), 0, sqrt(sigma2))
		
	}
	
		if(iter == 1000) print(date())
	}

print(date())

	###################################################################################
	
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

	###################################################################################

	# Log likelihood and BPIC calculation
	
	LL = median(-nlog_likelihood_store) ## Log likelihood
	
	D_bar = 2*mean(nlog_likelihood_store)

	mean_beta = apply(beta_store, 2, mean)

	mean_phi = apply(phi_store, 2, mean)
	
	mean_psi = apply(psi_store, 2, mean)
    
	mean_sigma2 = mean(sigma2_store)
	
	mean_X_beta = c(crossprod(t(covariate), c(mean_beta)))
	
	mean_resid = sweep(data, 2, (mean_X_beta + mean_psi)) - mean_phi

	D_theta_bar_mean = 2*(n_non_missing*log(mean_sigma2)/2 + sum((mean_resid[-pos_0])^2)/(2*mean_sigma2))

	BPIC = 3*D_bar - 2*D_theta_bar_mean

	#########################################################################################

	## Moran's I test

	vec = which((pos_0_month_length/S) < 0.5)  ## months which have more than 50% available data. Moran's I are calculated only for these months 
	
	residual_store = array(0, dim = c(S, T, nmc/nthin))
	
	for (i in 1:(nmc/nthin) ) residual_store[ , , i] = sweep(data , 2 , (colSums(t(covariate) * beta_store[ i , ]) + psi_store[i , ] ) ) - phi_store[i , ]
	
	data_na = matrix(0, S ,T)

	data_na[pos_0] = NA 

	med_residual = apply(residual_store, c(1,2), median)
	
	med_residual = med_residual + data_na

	moran_test = sapply(vec, function(i){unlist(Moran.I(med_residual[,i], W_mat, scaled = FALSE, na.rm = TRUE, alternative = "greater"))})

	#############################################################################################################################################
	
	## saving all necessary outputs to create first two rows in Table 1
	
	fill_up_missing_data$LL_BPIC = round(c(LL, BPIC), 3)
	
	fill_up_missing_data$moran_I = cbind(vec, round(t(moran_test)[ , c(1, 4)],3) )
	
	save(fill_up_missing_data, file = paste0("Inputs_and_Outputs/Model_IIIA_TC_", TC ,"_outputs.Rdata"))