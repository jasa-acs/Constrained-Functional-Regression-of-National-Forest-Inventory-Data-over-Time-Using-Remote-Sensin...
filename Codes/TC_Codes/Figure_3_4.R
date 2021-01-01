############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_TC_analysis.R that contains all necessary functions.

## the output of Aggregation.R that is "coarse_selected_data_with_missing_values.Rdata"
## "coarse_selected_data_with_missing_values.Rdata" contains a list coarse_creation that consists of two components:

## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.

## (ii) The second component contains the spatiotemporal data in form of a 10000x120x3 array where first and second dimensions correspond to pixels and time points (months), respectively. Along the third dimension, it contains values for TC1 (Brightness), TC2 (Greenness) and TC3 (Wetness). The occurrences of Zero (0) within this array indicate missing data at those space-time combinations.

## The output of TC_Model_BEST.R for all TCs where BEST should be replaced by the best TC model number (IV in our case) from the five candidate model numbers from Appendix A.1.1 in the supplementary materials - I, II, IIIA, IIIB and IV.

############################################################################################
########################		Output     ##################################################
############################################################################################

## The outputs are Figure 3 and Figure 4

###################################################################################################
###################################################################################################
###################################################################################################

library(ggplot2)
library(sp)
library(grid)
library(gridBase)
library(lattice)
library(proj4)
library(RColorBrewer) # for color

source("TC_Codes/Functions_for_TC_analysis.R")

########################################################################################
########################################################################################

## All TC features are in log scale and before doing the log transformation TC2 and TC3 are scaled to the same as that of TC1
## So before plotting the TC features, it is required to convert them into original scale

########################################################################################
########################################################################################

## Calculating necessary values to transformation the TC features into original scale

load("Inputs_and_Outputs/coarse_selected_data_with_missing_values.Rdata")

data_TC1 = coarse_creation$coarse_data[,,1]

data_TC2 = coarse_creation$coarse_data[,,2]

data_TC3 = coarse_creation$coarse_data[,,3]

S = dim(data_TC1)[1]

T = dim(data_TC1)[2]

pos_0 = which(data_TC1 == 0) 

range_TC1 = diff(range(data_TC1[-pos_0]))

range_TC2 = diff(range(data_TC2[-pos_0]))

range_TC3 = diff(range(data_TC3[-pos_0]))

min_TC2 =  min( data_TC2[-pos_0]* range_TC1/range_TC2)

min_TC3 =  min( data_TC3[-pos_0]* range_TC1/range_TC3)

##############################################################################################

pos_0_month_length = apply(data_TC1, 2, function(x) length(which(x == 0)) )

month_all = which(pos_0_month_length == 0)

pos_0_month_length = pos_0_month_length[-month_all]

group = rep(1:length(pos_0_month_length),  pos_0_month_length)

##############################################################################################
##############################################################################################

## loading the complete TC features after filling in the missing observations using the best model 

best_model = readline(" Best model number from I, II, IIIA, IIIB and IV ")

locations_all = coarse_creation$coarse_location ## location of the cells

TC_all = array (0, dim = c(S, T, 3))

uncertainty_HPD_list = list()

for(TC in 1:3)
	{

	load(paste0("Inputs_and_Outputs/Model_", best_model , "_TC_", TC , "_outputs.Rdata"))
	
	TC_all[ , , TC] = eval(parse(text = paste0("fill_up_missing_data$Transformed_TC", TC )))
	
	uncertainty_HPD_list[[TC]] = eval(parse(text = paste0("fill_up_missing_data$uncertainty_TC", TC )))
	
	}

TC1 = TC_all[ , , 1]

TC2 = TC_all[ , , 2]

TC3 = TC_all[ , , 3]

uncertainty_HPD_mat = do.call(cbind, uncertainty_HPD_list)

uncertainty_HPD = rowsum(uncertainty_HPD_mat, group = group)/pos_0_month_length

############################################################################################
########################		Figure 3     ###############################################
############################################################################################

## transforming TC features into their original scale 

TC1_original = exp(TC1)

TC2_original = (exp(TC2) - min(data_TC1[-pos_0]) + min_TC2)*range_TC2/range_TC1

TC3_original = (exp(TC3) - min(data_TC1[-pos_0]) + min_TC3 )*range_TC3/range_TC1

## estimating monthly 10-year average TC features

TC1_original_mean = matrix(0, dim(TC1)[1], 12)

TC2_original_mean = matrix(0, dim(TC2)[1], 12)

TC3_original_mean = matrix(0, dim(TC3)[1], 12)

for (t in 1:12)
{

TC1_original_mean[,t] = rowMeans(TC1_original[ , (t + (0:9)*12)])

TC2_original_mean[,t] = rowMeans(TC2_original[ , (t + (0:9)*12)])

TC3_original_mean[,t] = rowMeans(TC3_original[ , (t + (0:9)*12)])
}

## Following calculations helps in plotting TC features 

TC1_mean =  log(TC1_original_mean)

TC2_mean = log(TC2_original_mean*range_TC1/range_TC2 - min_TC2 + min(data_TC1[-pos_0]))

TC3_mean = log(TC3_original_mean*range_TC1/range_TC3 - min_TC3 + min(data_TC1[-pos_0]))

month_names = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

pdf("Inputs_and_Outputs/Figure_3.pdf", height = 15, width = 20)

#######################################################################################################
#######################################################################################################

## TC 1
 
   TC_number = 1
   
   p1 =  spatial_plot(locations_all,TC1_mean[ ,c(9:12,8:5,1:4)],log_scale=FALSE,grey_image=FALSE,titleplot="",layout_value=c(4,3,1),subplot=month_names[c(9:12,8:5,1:4)])

   L=grid.layout(1,1)

   pushViewport(viewport(layout=L))  #enter the device/layout

   pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))  #enter the first column

   print(p1,newpage=F)

   grid.newpage()


#######################################################################################################
#######################################################################################################

## TC 2
 
   TC_number = 2
   
   p1 =  spatial_plot(locations_all,TC2_original_mean[ ,c(9:12,8:5,1:4)],log_scale=FALSE,grey_image=FALSE,titleplot= "" ,layout_value=c(4,3,1),subplot=month_names[c(9:12,8:5,1:4)])
    
   L=grid.layout(1,1)
    
   pushViewport(viewport(layout=L))  #enter the device/layout
    
   pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))  #enter the first column
    
   print(p1,newpage=F)
    
   grid.newpage()

   
#######################################################################################################
#######################################################################################################

## TC 3
 
   TC_number = 3
	  
   p1 =  spatial_plot(locations_all,TC3_original_mean[ ,c(9:12,8:5,1:4)],log_scale=FALSE,grey_image=FALSE,titleplot= "",layout_value=c(4,3,1),subplot=month_names[c(9:12,8:5,1:4)])
    
   L=grid.layout(1,1)
    
   pushViewport(viewport(layout=L))  #enter the device/layout
    
   pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))  #enter the first column
    
   print(p1,newpage=F)
    
dev.off ( )

############################################################################################
########################		Figure 4     ###############################################
############################################################################################

missing_proportion_cut_off = c(0, .1, .3, .5, .7, .9, .99, 1)
pos_list = list()
missing_proportion_interval = c()
for (i in 1:(length(missing_proportion_cut_off) - 1)) 
	{
		missing_proportion_interval = c(missing_proportion_interval, paste0(missing_proportion_cut_off[i]*100, "-", missing_proportion_cut_off[i+1]*100))
		pos_list[[i]] = which((pos_0_month_length/S) > missing_proportion_cut_off[i]  & (pos_0_month_length/S) <= missing_proportion_cut_off[i+1])

	}

pos_list_length = sapply(pos_list, length)

group_1 = rep(1:length(pos_list),  pos_list_length)

mean_uncertainty_HPD = matrix(0, ncol = 3)
for (i in 1:length(pos_list)) mean_uncertainty_HPD = rbind(mean_uncertainty_HPD, colMeans(uncertainty_HPD[pos_list[[i]] , ]))
mean_uncertainty_HPD = mean_uncertainty_HPD[-1, ]

# create a dataset
missing_prop <- rep(missing_proportion_interval, each = 3)
TCS <- rep(c("TC1" , "TC2" , "TC3") , length(missing_proportion_interval))
uncertainty = c(t(mean_uncertainty_HPD))
data <- data.frame(missing_prop, TCS, uncertainty)
 
pdf("Inputs_and_Outputs/Figure_4.pdf", height = 4, width = 10)

figure_4 <- ggplot(data, aes(fill=TCS, y=uncertainty, x=missing_prop)) + 
			geom_bar(position="dodge", stat="identity", width = .5)+
			labs(x = "Monthly percentage of missing data", y = "Uncertainty") +
			scale_x_discrete(labels = c(expression(phantom(x) <=10), missing_proportion_interval[2:6], "100"))+
			scale_fill_manual("", values = c("yellow3", "springgreen3", "lightskyblue"))+
			theme(legend.position = c(.08, .8),
			legend.title = element_blank(),
			text = element_text(size=15))

print(figure_4)

dev.off()