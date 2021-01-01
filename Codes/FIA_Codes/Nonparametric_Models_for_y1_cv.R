############################################################################################
########################		Input     ##################################################
############################################################################################

## The output of TC_best_output_merge_basal_area.R that is "best_TC_model_outputs_and_basal_area_data.Rdata" where BEST should be replaced by the best TC model number (IV in our case) from the five candidate model numbers from Appendix A.1 in the supplementary materials - I, II, IIIA, IIIB and IV.

## best_TC_model_outputs_and_basal_area_data.Rdata contains a list reconstructed_TC_BALIVE consisting of five components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the complete TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 
## (iii) The fifth component is the live tree basal area in the form of a 10000x120 matrix where rows and colums correspond to pixels and time points (months), respectively. The coccurance of -99 within this matrix indicates unobserved data at those space-time combinations.


## one output of Basal_area_cv_positions_selction.R that is random_positions_for_CV_for_y1_model.Rdata
## random_positions_for_CV_for_y1_model.Rdata contains a list test_pos_list of 36 components where each component contains the postions of both zero, non-zero live tree basal area data those will be considered as test data in the correspondig holdout cross-validation method.


############################################################################################
########################		Output     #################################################
############################################################################################

## the output is "CV_outputs_nonparametric_y1_models.Rdata" that is a list containing ROC areas for support vector machine (SVM), stochastic gradient boosting (SGB) and generalized additive model (GAM).
 
################################################################################################################
#####################     Cross Validation     #################################################################
################################################################################################################

set.seed(77777)

library(e1071)			## For support vector machine
library(gbm)			## For stochastic gradient boosting 
library(gam)			## For generalized additive model
library(verification) 	## to draw ROC curve

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

	pred_mat_full <- t(pred_mat_full)
		
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
	
	pred_mat_all = pred_mat_full[pos_avail, ]
		
	count_all = count_full[pos_avail]

	n_CV = 36		## total number of replication of holdout method 

	first_year = ((1:S) - 1)*T + 1

	n_total = length(count_all)

	l_test = 15	## number of samples that are considered as tested data
	
	svm_count_prob_CV = matrix(0,nrow = 1, ncol = 3)

	sgb_count_prob_CV = matrix(0,nrow = 1, ncol = 2)

	gam_count_prob_CV = matrix(0,nrow = 1, ncol = 2)
	
	
	for(cv in 1:n_CV)
		
	{
	
	set.seed(77777)

	pos = test_pos_list[[cv]]

	count = count_all[-c(pos)]
	
	count_test = count_all[c(pos)]

	pred_mat = pred_mat_all[-c(pos), ]
	
	pred_mat_test = pred_mat_all[c(pos), ]
	
	data_train <- data.frame(y = count, pred_mat)
	
	data_test <- data.frame(pred_mat_test)
	
	###########################################################################################################
	###########################################################################################################

	## Fitting support vector machine (svm)

	set.seed(77777)

	svm_y1_model <- svm(y ~ . , data = data_train, type = "C", probability=TRUE, )
	
	svm_count_pred <- predict(svm_y1_model, newdata = data_test, probability=TRUE, , type = "response")

	############################################################################################################
	
	svm_output = cbind( count_test, attr(svm_count_pred, "probabilities")[ , 1], svm_count_pred )

	svm_output[ , 3] = svm_output[ , 3] - 1
	
	svm_count_prob_CV = rbind(svm_count_prob_CV, svm_output)

	#####################################################################################
	#####################################################################################

	## Fitting stochastic gradient boosting (sgb)

	set.seed(77777)

	sgb_y1_model <- gbm(y ~ . , data = data_train,  distribution = "bernoulli")

	## Estimating the optimal number of boosting iterations for a gbm object
	
	best.iter <- gbm.perf(sgb_y1_model, plot.it = FALSE)

	#print(best.iter)

	sgb_prob_test <- predict(sgb_y1_model, newdata = data_test, n.trees = best.iter, type = "response")

	############################################################################################################
	
	sgb_output = cbind(count_test, sgb_prob_test)

	sgb_count_prob_CV = rbind(sgb_count_prob_CV, sgb_output)

	#######################################################################################################	#######################################################################################################

	## Fitting generalized additive model (gam)

	set.seed(77777)

	xnam <- paste0( "s(X", 1:ncol(pred_mat), ")" )
	
	fmla <- as.formula(paste0("y ~ ", paste0(xnam, collapse= "+")))

	gam_y1_model <- gam(fmla, data = data_train,  family = binomial)
	
	gam_prob_test <- predict(gam_y1_model, newdata = data_test, type = "response")

	#########################################################################################################
	
	gam_output = cbind(count_test, gam_prob_test)

	gam_count_prob_CV = rbind(gam_count_prob_CV, gam_output)

	}
	
	#####################################################################################
	#####################################################################################

	## support vector machine (svm) outputs

	svm_count_prob_CV = svm_count_prob_CV[-1, ]	## count, probability and prediction storage

	svm_ROC_area = round(roc.area(svm_count_prob_CV[ , 1], svm_count_prob_CV[ , 2])$A, digits = 3)	## ROC area 

	## stochastic gradient boosting (sgb) outputs

	sgb_count_prob_CV = sgb_count_prob_CV[-1,]	## count and probability storage

	sgb_ROC_area = round(roc.area(sgb_count_prob_CV[,1], sgb_count_prob_CV[,2])$A, digits = 3)	## ROC area 

	## generalized additive model (gam) outputs

	gam_count_prob_CV = gam_count_prob_CV[-1,]	## count and probability storage
	
	gam_ROC_area = round(roc.area(gam_count_prob_CV[,1], gam_count_prob_CV[,2])$A, digits = 3)	## ROC area 

	## storing all nonparametric models outputs
	
	CV_outputs_nonparametric_y1_models <- list(svm_ROC_area = svm_ROC_area, sgb_ROC_area = sgb_ROC_area, gam_ROC_area = gam_ROC_area)

	save(CV_outputs_nonparametric_y1_models, file = "Inputs_and_Outputs/CV_outputs_nonparametric_y1_models.Rdata")
