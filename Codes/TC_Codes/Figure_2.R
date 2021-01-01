############################################################################################
########################		Input     ##################################################
############################################################################################

## Functions_for_TC_analysis.R that contains all necessary functions.

## the output of Aggregation.R that is "coarse_selected_data_with_missing_values.Rdata"
## "coarse_selected_data_with_missing_values.Rdata" contains a list coarse_creation that consists of two components:

## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.

## (ii) The second component contains the spatiotemporal data in form of a 10000x120x3 array where first and second dimensions correspond to pixels and time points (months), respectively. Along the third dimension, it contains values for TC1 (Brightness), TC2 (Greenness) and TC3 (Wetness). The occurrences of Zero (0) within this array indicate missing data at those space-time combinations.

############################################################################################
########################		Output     #################################################
############################################################################################

## Figure 2

############################################################################################################
############################################################################################################
############################################################################################################

set.seed(7777)

library(sp)
library(grid)
library(gridBase)
library(lattice)
library(proj4)
library(RColorBrewer) # for color
library(ggplot2)

source("TC_Codes/Functions_for_TC_analysis.R")

load("Inputs_and_Outputs/coarse_selected_data_with_missing_values.Rdata") ## loading output of Aggregation.R that contains a list coarse_creation

TC = 1

data = coarse_creation$coarse_data[,,TC]

S = dim(data)[1]

T = dim(data)[2]

pos_0_list = apply(data, 1, function(x) which(x == 0))

pos_0_cell_ratio = sapply(pos_0_list, length)/T	## proportion of missing data at each cell

pos_0_list_month = apply(data, 2, function(x) which(x == 0))

pos_0_month_ratio = sapply(pos_0_list_month, length)/S		## proportion of missing data at each month

location = coarse_creation$coarse_location

month_cal = function(t) t - floor((t-1)/12)*12

month_name = c("January", "February", "March", "April","May", "June", "July", "August", "September", "October", "November", "December")

Months = c(1:120)

pdf("Inputs_and_Outputs/Figure_2.pdf", height = 6, width = 6)

###############################################################################################################
###############################################################################################################
###############################################################################################################

## Figure 2 (a)

p1 <-ggplot(data= data.frame(Months, pos_0_month_ratio), aes(x=Months, y=pos_0_month_ratio)) +
	geom_bar(stat="identity")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin=unit(c(1,1,1,1),"cm")) +
	scale_x_continuous( labels = as.character(c(1,seq(12,120, by = 12))), breaks = c(1,seq(12,120, by = 12)) )+
	labs(x = "Months", y = "Proportion of missing data")

print(p1)

grid.newpage()

###############################################################################################################
###############################################################################################################
###############################################################################################################

## Figure 2 (b)

L=grid.layout(1,1)

pushViewport(viewport(layout=L))  #enter the device/layout

p1 = spatial_plot_Figure_2b(location, pos_0_cell_ratio, log_scale=FALSE, grey_image=TRUE, titleplot=paste("",sep=""), layout_value=c(1,1,1),subplot=none)

pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))  

print(p1,newpage=F)

popViewport()

dev.off()