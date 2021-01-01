############################################################################################
########################		Input     ##################################################
############################################################################################

## The output of Unconstrained_y2_model.R that is "unconstrained_y2_model_output.Rdata".
## "unconstrained_y2_model_output.Rdata" contains a list named outputs consisting of two components:
## (i) the first component is the posterior median likelihoods at each available live tree basal area measurement.
## (ii) the second component is a vector of size 2 with log likelihood and BPIC, respectively. 


## The output of Constrained_y2_model.R  that is "constrained_y2_model_output.Rdata".
##"constrained_y2_model_output.Rdata" contains a list named outputs consisting of three components:
## (i) the first component is the posterior median likelihoods at each available live tree basal area measurement.
## (ii) the second component is a vector of size 2 with log likelihood and BPIC, respectively. 
## (iii) the third component is the MCMC outputs of predicted live tree basal area measurement of the form of a TxSx(nmc/nthin) array where first and second dimensions correspond to time points (years) and pixels, respectively. Along the third dimension it contains the corresponding MCMC outputs after burn in and thinning.


## One output of Model_for_y1_cv.R that is "ROC_area_for_y1_model.Rdata".
## "ROC_area_for_y1_model.Rdata" contains the ROC area.

## The output of Unconstrained_y2_model_cv.R named "CV_outputs_unconstrained_y2_model.csv".
## "CV_outputs_unconstrained_y2_model.csv" is a 3x1 matrix containing mean of absolute bias, uncertainty and empirical coverage, respectively.

## The output of Constrained_y2_model_cv.R that is "CV_outputs_constrained_y2_model.Rdata".
## "CV_outputs_constrained_y2_model.Rdata" is list of 5 measures containing mean absolute bias, mean uncertainty, empirical coverage, relative MSPE (Mean Square Prediction Error) and correlation between test and prediction data.

## The output of Non_Parametric_Models_for_y1_cv.R that is "CV_outputs_nonparametric_y1_models.Rdata".
## "CV_outputs_nonparametric_y1_models.Rdata" is a list containing ROC areas for support vector machine (SVM), stochastic gradient boosting (SGB) and generalized additive model (GAM). 

## The output of Nor_Parametric_Models_for_y2_cv.R that is "CV_outputs_nonparametric_y2_models.csv".
## CV_outputs_nonparametric_y2_models.csv is a 2X3 matrix where columns represends support vector machine (SVM), stochastic gradient boosting (SGB) as well as generalized additive model (GAM), and rows represends relative MSPE (Mean Square Prediction Error) as well as correlation between test and prediction data, respectively for the corresponding method.

############################################################################################
########################		Output     #################################################
############################################################################################

## Table 1, 2 and Figure 5

############################################################################################
############################################################################################

## loading necessary packages

library(ggplot2)
library(reshape2)

############################################################################################
############################################################################################

## outputs of Unconstrained_y2_model.R

	load("Inputs_and_Outputs/unconstrained_y2_model_output.Rdata")	## loading outputs of Unconstrained_y2_model.R 

	unconstrained_likelihood = outputs$unconstrained_likelihood
	
	unconstrained_LL_BPIC = outputs$LL_BPIC
	
############################################################################################
############################################################################################

## outputs of Constrained_y2_model.R

	load("Inputs_and_Outputs/constrained_y2_model_output.Rdata")	## loading outputs of Unconstrained_y2_model.R 

	constrained_likelihood = outputs$constrained_likelihood
	
	constrained_LL_BPIC = outputs$LL_BPIC

############################################################################################
############################################################################################

## outputs of Unconstrained_y2_model_CV.R

	unconstrained_CV_output = read.csv("Inputs_and_Outputs/CV_outputs_unconstrained_y2_model.csv", header = T)[1:3,1]

############################################################################################
############################################################################################

## outputs of Constrained_y2_model_cv.R

	load("Inputs_and_Outputs/CV_outputs_constrained_y2_model.Rdata")	## loading outputs of Unconstrained_y2_model_cv.R 

	constrained_CV_output = matrix(c(CV_outputs_constrained_y2_model$abs_bias, CV_outputs_constrained_y2_model$uncertainty, CV_outputs_constrained_y2_model$empirical_coverage), ncol = 1)
	
############################################################################################
########################		Table 1     ###############################################
############################################################################################

	LL_BPIC = cbind(unconstrained_LL_BPIC, constrained_LL_BPIC)

	CV_outputs = cbind(unconstrained_CV_output, constrained_CV_output)

	Table_1 = round(rbind (LL_BPIC, CV_outputs), 3)

	a = c("LL", "BPIC", "Absolute Bias", "Uncertainty", "Empirical Coverage")
	
	Table_1 = data.frame(cbind(a, Table_1))

	colnames(Table_1) = c("Comparison Statistic", "Unconstrained Model", "Constrained Model")

	write.csv(Table_1, file = "Inputs_and_Outputs/Table_1.csv", row.names = FALSE) 
	
############################################################################################
########################		Figure 5     ###############################################
############################################################################################

	order_pos = order(unconstrained_likelihood)
	
	data_all = data.frame(list(sample = 1:length(constrained_likelihood), Constrained_Model = constrained_likelihood[order_pos], Unconstrained_Model = unconstrained_likelihood[order_pos]))

	names(data_all) <- c("sample", "Constrained Model", "Unconstrained Model")

	data_all_melt <- melt(data_all, id.vars = 'sample', variable.name = 'Model', value.name = 'likelihood')

	pdf("Inputs_and_Outputs/Figure_5.pdf", height = 5, width = 7)
			
	p1 <- ggplot(data_all_melt, aes(x = sample, y = likelihood, shape = Model, color = Model))+
	geom_point() + 
	scale_shape_manual(values = c('Constrained Model' = 17, 'Unconstrained Model' = 16))+
	labs(x = "", y = "Likelihood") +
	theme(legend.position = c(0.8, 0.1), legend.title=element_blank()) +
	scale_color_manual(values=c("red", "blue"))
		
	print(p1)
		
	dev.off()

############################################################################################
########################		Table 2     ###############################################
############################################################################################

## loading outputs of Model_for_y1_cv.R

	load("Inputs_and_Outputs/ROC_area_for_y1_model.Rdata")	## loading outputs of Unconstrained_y2_model_cv.R 

	y1_model_ROC_area = ROC_area
	
## outputs of Constrained_y2_model_cv.R

	constrained_y2_model_CV_output = c(CV_outputs_constrained_y2_model$relative_MSPE, CV_outputs_constrained_y2_model$correlation)

## loading outputs of Nonparametric_Models_for_y1_cv.R

	load("Inputs_and_Outputs/CV_outputs_nonparametric_y1_models.Rdata")
	
	nonparametric_y1_models_ROC_area = c(CV_outputs_nonparametric_y1_models$svm_ROC_area, CV_outputs_nonparametric_y1_models$sgb_ROC_area, CV_outputs_nonparametric_y1_models$gam_ROC_area)
	
## loading outputs of Norparametric_Models_for_y2_cv.R

	nonparametric_y2_models_CV_output = read.csv("Inputs_and_Outputs/CV_outputs_nonparametric_y2_models.csv", header = T)[ , -1]
	
## Table 3

	Table_2 = round(rbind(c(y1_model_ROC_area, nonparametric_y1_models_ROC_area), cbind(constrained_y2_model_CV_output, nonparametric_y2_models_CV_output)), 3)

	colnames(Table_2) <- c("Proposed Model", "SVM", "SGB", "GAM")
	
	rownames(Table_2) <- c("ROC area", "Relative MSPE", "Correlation")
		
	write.csv(Table_2, file = "Inputs_and_Outputs/Table_2.csv")
