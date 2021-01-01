############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_FIA_data_analysis.R that contains all necessary functions.


## The output of TC_best_output_merge_basal_area.R that is "best_TC_model_outputs_and_basal_area_data.Rdata".
## best_TC_model_outputs_and_basal_area_data.Rdata contains a list reconstructed_TC_BALIVE consisting of five components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the complete TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 
## (iii) The fifth component is the live tree basal area in the form of a 10000x120 matrix where rows and colums correspond to pixels and time points (months), respectively. The coccurance of -99 within this matrix indicates unobserved data at those space-time combinations.


## The output of Model_for_y1.R which is "Model_y1_output.Rdata".
## "Model_y1_output.Rdata" contains the MCMC outputs of probability of the forest of the form of a TxSx(nmc/nthin) array prob_construct_MCMC where first and second dimensions correspond to time points (years) and pixels, respectively. Along the third dimension it contains the corresponding MCMC outputs after burn in and thinning.
 

## The output of Constrained_y2_model.R that is "constrained_y2_model_output.Rdata".
## "constrained_y2_model_output.Rdata" contains a list named outputs consisting of three components:
## (i) the first component is the posterior median likelihoods at each available live tree basal area measurement.
## (ii) the second component is a vector of size 2 with log likelihood and BPIC, respectively. 
## (iii) the third component is the MCMC outputs of predicted live tree basal area measurement of the form of a TxSx(nmc/nthin) array where first and second dimensions correspond to time points (years) and pixels, respectively. Along the third dimension it contains the corresponding MCMC outputs after burn in and thinning.

############################################################################################
########################		Output     #################################################
############################################################################################

## Figure 6, Figure 7 and Figure 8

################################################################################################################
################## 		Figure 6, 7, 8 		 ###################################################################
################################################################################################################

set.seed(77777)

## Loading necessary packages

library(sp)
library(grid)
library(gridBase)
library(lattice)
library(proj4)
library(RColorBrewer) # for color
library(matrixStats) # for colMedians
library(coda)

source("FIA_Codes/Functions_for_FIA_data_analysis.R")

load("Inputs_and_Outputs/best_TC_model_outputs_and_basal_area_data.Rdata")	## loading the data 

	locations_all = reconstructed_TC_BALIVE$locations	## location of the cells 

	S = nrow(locations_all)	## Total number of grid-cells in response

	T = 10	## Total number of time points (years here) in response

	########################################################################################################

	## Outputs of y1 model
	
	load("Inputs_and_Outputs/Model_y1_output.Rdata")
		
	prob_forest = prob_construct_MCMC
		
	rm(prob_construct_MCMC)
	
	nmc_by_nthin = dim(prob_forest)[3]
	
	count_temp = matrix(-99,  S,  T)

	pos_0 = which(reconstructed_TC_BALIVE$BALIVE == 0, arr.ind = TRUE)

	pos_0 = pos_0[which(pos_0[ , 2] >  0 & pos_0[ , 2] < 121), ]

	pos_1 = which(reconstructed_TC_BALIVE$BALIVE > 0, arr.ind = TRUE)

	pos_1 = pos_1[which(pos_1[ , 2] >  0 & pos_1[ , 2] < 121), ]

	for(s in 1:(dim(pos_1)[1])) count_temp[pos_1[s, 1], ceiling(pos_1[s, 2]/12)] = 1

	for(s in 1:(dim(pos_0)[1])) count_temp[pos_0[s, 1], ceiling(pos_0[s, 2]/12)] = 0

	count_temp = t(count_temp)
	
	count_full = ceiling(prob_forest - array(runif(dim(prob_forest)[1] * dim(prob_forest)[2]* dim(prob_forest)[3]), dim = dim(prob_forest)))
		
	dim(count_full) = c(T*S, nmc_by_nthin)
		
	nonforest = which(count_temp == 0)
		
	forest = which(count_temp == 1)
		
	count_full[nonforest, ] = 0
		
	count_full[forest, ] = 1
		
	dim(count_full) = c(T, S, nmc_by_nthin)

	########################################################################################################

	## cells where the basal area is available at two time points
	## In count_temp, 0 represents 0 basal area and 1 represents positive basal area 
	
	available_locations_basal_area = list()

	for (t in 1:T) available_locations_basal_area[[t]] =   sort(which(count_temp[t, ] != -99)) 

	## this stores all the cell locations, where is available at two time points	
	available_locations_unique_basal_area = sort(unique(unlist(available_locations_basal_area)) ) 
	
	########################################################################################################

	## Outputs of constrained y2 model
	
	load("Inputs_and_Outputs/constrained_y2_model_output.Rdata")
	
	basal_area_MCMC_exp = outputs$basal_area_MCMC_exp
	
	basal_area_exp_times_prob_forest = basal_area_MCMC_exp * count_full
		
	basal_area_summary_exp = t( apply(basal_area_exp_times_prob_forest, c(1, 2), median) )
		
	########################################################################################################

	pdf("Inputs_and_Outputs/Figure_6_7_8.pdf", height = 10, width = 21)
	
	## Figure 6 (Posterior median of yearly average live tree basal area)
			
	L=grid.layout(1,1)

	pushViewport(viewport(layout=L))  #enter the device/layout

	p1 = spatial_plot(locations_all, (basal_area_summary_exp[,c(ncol(basal_area_summary_exp):6,1:5)]), constant = 1, log_scale=TRUE, grey_image=FALSE, titleplot=paste("",sep=""),layout_value=c(5,2,1),subplot=as.character(c(2012:2008,2003:2007)))

	pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))  #enter the first column

	print(p1,newpage=F)

	popViewport()
		
	grid.newpage()

	########################################################################################################

	# Figure 7(a) (10-year averaged 90% HPD interval width)

	basal_area_exp_uncertainty_90 = t( apply(basal_area_exp_times_prob_forest, c(1, 2), function(x) diff(HPDinterval(mcmc(x), prob = 0.90)[1,] ) ) )

	mean_uncertainty <- rep(-9999, S)

	## available_locations_unique_basal_area contains the gridcells that have two observed basal area measurements
	## So, for these cells mean uncertainty are calculated based on the posterior predictions for the remaining eight years
	
	mean_uncertainty[available_locations_unique_basal_area] <- rowSums(basal_area_exp_uncertainty_90[available_locations_unique_basal_area, ])/(T-2)

	## mean uncertainty for the cells where there is no observed data
	
	mean_uncertainty[-available_locations_unique_basal_area] <- rowMeans(basal_area_exp_uncertainty_90[-available_locations_unique_basal_area, ])
	
	L=grid.layout(1,1)

	pushViewport(viewport(layout=L))  #enter the device/layout

	p1 = spatial_plot_uncertainty(locations_all, mean_uncertainty, log_scale = FALSE, grey_image=FALSE, titleplot=paste("",sep=""),layout_value=c(1,1,1))

	pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))  #enter the first column

	print(p1,newpage=F)

	popViewport()
		
	grid.newpage()

	#####################################################################################################		

	## Figure 7(b) (10-year averaged ratio of 90% HPD interval width and median for yearly average live tree basal area
	
	const = quantile(basal_area_summary_exp[basal_area_summary_exp>0], prob = .001)
		
	pos_rej = which(basal_area_summary_exp <= const)
		
	uncertainty_by_basal = matrix(-9999, nrow(basal_area_summary_exp), ncol(basal_area_summary_exp))
		
	uncertainty_by_basal[-pos_rej] = basal_area_exp_uncertainty_90[-pos_rej]/basal_area_summary_exp[-pos_rej]
		
	uncertainty_by_basal[pos_rej] = basal_area_exp_uncertainty_90[pos_rej]/const - (basal_area_exp_uncertainty_90[pos_rej]/const^2) *(basal_area_summary_exp[pos_rej] - const)
		
	mean_uncertainty_by_basal <- rep(-9999, S)

	## available_locations_unique_basal_area contains the gridcells that have two observed basal area measurements
	## So, for these cells mean of relative uncertainty are calculated based on the posterior predictions for the remaining eight years
	
	mean_uncertainty_by_basal[available_locations_unique_basal_area] <- rowSums(uncertainty_by_basal[available_locations_unique_basal_area, ])/(T-2)

	## mean of relative uncertainty for the cells where there is no observed data
	
	mean_uncertainty_by_basal[-available_locations_unique_basal_area] <- rowMeans(uncertainty_by_basal[-available_locations_unique_basal_area, ])
	
	L=grid.layout(1,1)

	pushViewport(viewport(layout=L))  #enter the device/layout

	p1 = spatial_plot_uncertainty(locations_all, mean_uncertainty_by_basal, constant = 1, log_scale = TRUE, grey_image=FALSE, titleplot=paste("",sep=""),layout_value=c(1,1,1))

	pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))  #enter the first column

	print(p1,newpage=F)

	popViewport()
		
	grid.newpage()

	#####################################################################################################		

	## Figure 8(a) Change in estimated live tree basla area averaged over 200=8-2012 compared to its average from the period 2003-2007

	five_year_change_basal_area = apply(basal_area_summary_exp[ ,6:10], 1, mean) - apply(basal_area_summary_exp[ ,1:5], 1, mean)

	L=grid.layout(1,1)
	
	pushViewport(viewport(layout=L))  #enter the device/layout

	p1 = spatial_plot_five_year_change(locations_all, five_year_change_basal_area, log_scale=FALSE, grey_image=FALSE, titleplot=paste("",sep=""),layout_value=c(1,1,1))

	pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))  #enter the first column

	print(p1,newpage=F)

	popViewport()

	grid.newpage()

	###############################################################################################
	
	## Figure 8(b) Standard deviation of yearly average live tree basal area estimates
	
	sd_temporally = apply(basal_area_summary_exp, 1, sd)
	
	L=grid.layout(1,1)

	pushViewport(viewport(layout=L))  #enter the device/layout

	p1 = spatial_plot_uncertainty(locations_all, sd_temporally, log_scale = TRUE, grey_image=FALSE, titleplot=paste("",sep=""),layout_value=c(1,1,1))

	pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))  #enter the first column

	print(p1,newpage=F)

	popViewport()

	dev.off()
