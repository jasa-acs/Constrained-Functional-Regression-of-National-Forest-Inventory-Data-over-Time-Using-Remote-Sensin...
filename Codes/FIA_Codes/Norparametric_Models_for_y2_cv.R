############################################################################################
########################		Input     ##################################################
############################################################################################

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

## the output is "CV_outputs_nonparametric_y2_models.csv" that is a 2X3 matrix where columns represends support vector machine (SVM), stochastic gradient boosting (SGB) as well as generalized additive model (GAM), and rows represends relative MSPE (Mean Square Prediction Error) as well as correlation between test and prediction data, respectively for the corresponding method.

################################################################################################################
#####################     Cross Validation     #################################################################
################################################################################################################

set.seed(77777)

library(e1071)		## For support vector machine
library(gbm)		## For stochastic gradient boosting 
library(gam)		## For generalized additive model

load("Inputs_and_Outputs/best_TC_model_outputs_and_basal_area_data.Rdata")	## loading the data 

load("Inputs_and_Outputs/random_positions_for_CV_for_y2_model.Rdata")

	locations_all = reconstructed_TC_BALIVE$locations	## location of the cells 

	S = nrow(locations_all)	## Total number of grid-cells in response

	T = 10	## Total number of time points (years here) in response

	n_TC = 3	## Total number TC features

	##########################################################################################################
	
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

	pred_mat_full = t(pred_mat_full)
	
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

	pred_mat_all = pred_mat_full[pos_avail, ]
		
	y_all = y_full[pos_avail]

	y_train_mean = c()
	
	n_CV = 36		## total number of replication of holdout method 

	l_test = 10	## number of samples that are considered as tested data

	## Store the simulated test data
		
	svm_outputs = matrix(0, 1, 2)

	sgb_outputs = matrix(0, 1, 2)

	gam_outputs = matrix(0, 1, 2)

	for(cv in 1:n_CV)
		
	{
	
	set.seed(77777)

	pos = test_pos_list[[cv]]

	y = y_all[-c(pos)]	## train data
	
	y_test = y_all[c(pos)]	## test data
	
	n = length(y)

	## mean of the train data after transforming into original scale by taking exponential
	## the given unit of basal area is ft^2/acre. To transform it into m^2/hectare it is needed to mulitply the basal area by 0.229568411	
	
	y_train_mean[cv] = mean(exp(y)*0.229568411)

	## covariates for train data

	pred_mat = pred_mat_all[-c(pos), ]

	## covariates for test data
	
	pred_mat_test = pred_mat_all[c(pos), ]

	data_train <- data.frame(y = y, pred_mat)
	
	data_test <- data.frame(pred_mat_test)

	###########################################################################################################
	###########################################################################################################

	## Fitting support vector machine (svm)

	set.seed(77777)

	svm_y2_model <- svm(y ~ . , data = data_train, type = "eps")

	svm_y_pred <- predict(svm_y2_model, newdata = data_test, type = "response")

	###########################################################################################
	
	## storing the outputs
	
	svm_outputs = rbind(svm_outputs, cbind( y_test, svm_y_pred ) )

	#####################################################################################
	#####################################################################################

	## Fitting stochastic gradient boosting (sgb)

	set.seed(77777)

	sgb_y2_model <- gbm(y ~ . , data = data_train,  distribution = "gaussian")

	## Estimating the optimal number of boosting iterations for a gbm object
	
	best.iter <- gbm.perf(sgb_y2_model, plot.it = FALSE)

	sgb_y_pred <- predict(sgb_y2_model, newdata = data_test, n.trees = best.iter, type = "response")

	###########################################################################################
	
	## storing the outputs
	
	sgb_outputs = rbind( sgb_outputs, cbind( y_test, sgb_y_pred ) )
	
	#######################################################################################################	#######################################################################################################

	## Fitting generalized additive model (gam)

	set.seed(77777)
	
	xnam <- paste0( "s(X", 1:ncol(pred_mat), ")" )
	
	fmla <- as.formula(paste0("y ~ ", paste0(xnam, collapse= "+")))
	
	gam_y2_model <- gam(fmla, data = data_train,  family = gaussian)

	gam_y_pred <- predict(gam_y2_model, newdata = data_test, type = "response")

	###########################################################################################
	
	## storing the outputs
	
	gam_outputs = rbind(gam_outputs, cbind( y_test, gam_y_pred ) )
	
}

	###########################################################################################################
	###########################################################################################################

	## support vector machine (svm) outputs

	## transforming basal area data into original scale by taking exponential	
	## the given unit of basal area is ft^2/acre. To transform it into m^2/hectare it is needed to mulitply the basal area by 0.229568411	

	svm_outputs = exp(svm_outputs[-1, ])*0.229568411	## transforming  ft^2/acre  to m^2/hectare 

	colnames(svm_outputs) = c("y_test", "y_pred")
	
	y_train_mean_full = rep(y_train_mean, each = l_test)
	
	svm_MSPE_pred = mean( (svm_outputs[ , 1] - svm_outputs[ , 2])^2 )
	
	svm_MSPE_mean = mean( (svm_outputs[ , 1] - y_train_mean_full)^2 )
	
	svm_relative_MSPE = svm_MSPE_pred/svm_MSPE_mean
	
	svm_correlation = cor(svm_outputs[ , 1], svm_outputs[ , 2])

	#####################################################################################
	#####################################################################################

	## stochastic gradient boosting (sgb) outputs

	## transforming basal area data into original scale by taking exponential	
	## the given unit of basal area is ft^2/acre. To transform it into m^2/hectare it is needed to mulitply the basal area by 0.229568411	

	sgb_outputs = exp(sgb_outputs[-1, ])*0.229568411	## transforming  ft^2/acre  to m^2/hectare 

	colnames(sgb_outputs) = c("y_test", "y_pred")
	
	sgb_abs_bias = mean(abs(sgb_outputs[ , 1] - sgb_outputs[ , 2]))

	sgb_MSPE_pred = mean( (sgb_outputs[ , 1] - sgb_outputs[ , 2])^2 )

	sgb_MSPE_mean = mean( (sgb_outputs[ , 1] - y_train_mean_full)^2 )
	
	sgb_relative_MSPE = sgb_MSPE_pred/sgb_MSPE_mean
	
	sgb_correlation = cor(sgb_outputs[ , 1], sgb_outputs[ , 2])

	#######################################################################################################	#######################################################################################################

	## generalized additive model (gam) outputs

	## transforming basal area data into original scale by taking exponential	
	## the given unit of basal area is ft^2/acre. To transform it into m^2/hectare it is needed to mulitply the basal area by 0.229568411	

	gam_outputs = exp(gam_outputs[-1, ])*0.229568411	## transforming  ft^2/acre  to m^2/hectare 

	colnames(gam_outputs) = c("y_test", "y_pred")
	
	gam_abs_bias = mean(abs(gam_outputs[ , 1] - gam_outputs[ , 2]))

	gam_MSPE_pred = mean( (gam_outputs[ , 1] - gam_outputs[ , 2])^2 )
	
	gam_MSPE_mean = mean( (gam_outputs[ , 1] - y_train_mean_full)^2 )
	
	gam_relative_MSPE <- gam_MSPE_pred/gam_MSPE_mean
	
	gam_correlation <- cor(gam_outputs[ , 1], gam_outputs[ , 2])
	
	#####################################################################################
	#####################################################################################

	outputs <- round(matrix( c(svm_relative_MSPE, svm_correlation, sgb_relative_MSPE, sgb_correlation, gam_relative_MSPE, gam_correlation), 2, 3 ), 3)
	
	colnames(outputs) <- c("SVM", "SGB", "GAM")

	rownames(outputs) <- c("Relative MSPE", "Correlation")
		
	write.csv(outputs, file = "Inputs_and_Outputs/CV_outputs_nonparametric_y2_models.csv")
