## This code is optional, it needs to be run only if the user wants to generate a hypothetical basal area data different from what is already included in Inputs_and_Outputs folder by changing the seed specified below and/or using a different set of values for model parameters.

############################################################################################
########################		Input     ##################################################
############################################################################################

## The output of TC_best_model_output.R that is "best_TC_model_outputs.Rdata"

## best_TC_model_outputs.Rdata contains a list reconstructed_TC consisting of four components:
## (i) The first component is a 10000x2 location matrix with its rows representing coordinates of the pixel centers from the study area described in Section 2 of the manuscript.
## (ii) The second, third and fourth components contain the complete TC1, TC2 and TC3 data, respectively, (after filling in the missing values using the best performing model) in form of a 10000x120 matrix where  rows and colums correspond to pixels and time points (months), respectively. 

## parameters_used_to_simulate_basal_area_data.Rdata that contains a list named beta_sigma2 consisting of three components:
## (i) The first component is a vector containing regression coefficients to simulate y^(1) data described in Section 3.1.
## (ii) The second and third components are a vector containing regression coefficients and a scalar containing variance parameter, respectively, to simulate y^(2) data from the constrained model described in Section 3.2.

############################################################################################
########################		Output     #################################################
############################################################################################


## the output "hypothetical_basal_area_data.csv" is a matrix of five columns where first two entries in each row represents the simulated coordinates of the field plot of live tree basal area data, next two entries represents the simulated time points with year and month number (a number from 1 to 12) of the correspondig year, respectively and the last entry shows the corresponding simulated live tree basal area measurement.

 
################################################################################################################
#######################     basal area data simulation     #####################################################
################################################################################################################

# If the user wants to simulate a hypothetical basal area data different from what is already included in Inputs_and_Outputs folder, the seed must be changed in the following line

set.seed(7777)

load("Inputs_and_Outputs/best_TC_model_outputs.Rdata")	## loading the data 

locations_all = reconstructed_TC$locations	## location of the cells 

S = nrow(locations_all)	## Total number of grid-cells in response

T = 10	## Total number of time points (years here) in response

n_TC = 3	## Total number TC features

n_total = S*T	## Total number of data points in response

###############################################################################################################

## Construction of functional X matrix

temp_col = ((1:T) - 1)*12 + 1

col_no = c(sapply(0:11, function(i) temp_col + i))

TC_data = array(-999, c(S, 12*T, n_TC))

for(TC in 1:n_TC) TC_data[ , , TC] = reconstructed_TC[[TC+1]]

TC_data = TC_data[ , col_no, ]

dim(TC_data) = c(S, T, 3* 12)

TC_data = aperm(TC_data, c(3, 2, 1))

pred_mat_full = TC_data

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

## If one does not want to use the following parameter values and instead wants to specify own set of parameter values, those need to be mentioned here.

load("Inputs_and_Outputs/parameters_used_to_simulate_basal_area_data.Rdata") ## beta and sigma2 used to simulate basal area data

beta_y1 = beta_sigma2$beta_y1

beta_y2 = beta_sigma2$beta_y2

sigma2_y2 = beta_sigma2$sigma2_y2

##############################################################################################

## since the average number of observations per year in the original dataset is about 29, we randomly select number of observations per year from a +- 5 interval around that number.

## to account for a 5-year remeasurement period, we made sure that from all of these locations data were recorded twice at a 5 year interval, although , in actual dataset, this remeasurement period was seen to vary between 4 to 6 years.

S_sim_year = sample(24:34, 5, replace = TRUE)

cum_S_sim_year = c(0,cumsum(S_sim_year))

loc_pos_sim = sample(1:S, sum(S_sim_year), replace = FALSE)     ## simulating locations where data is available

loc_pos_sim_list = split(loc_pos_sim, rep(1:5, S_sim_year))

year_sim = rep(2003:2007, S_sim_year)
             
BALIVE_loc_year = cbind( c(loc_pos_sim, loc_pos_sim), c(year_sim, year_sim + 5))

n = nrow(BALIVE_loc_year)

##############################################################################################

BALIVE = matrix(-99, n, 5)

colnames(BALIVE) = c("X", "Y", "Year", "Month Number", "Basal area")

x_merge = 16		## The number of adjacent cells in TC features to be aggregated in X coordinate

y_merge = 16		## The number of adjacent cells in TC features to be aggregated in Y coordinate

x_window = 30*(x_merge/2)    ## each pixel tile of TC features is 30X30 meter^2 before aggregation

y_window = 30*(y_merge/2)    ## each pixel tile of TC features is 30X30 meter^2 before aggregation

temp_loc = locations_all[BALIVE_loc_year[ ,1], ]

## simulating locations for basal area data

BALIVE[ , 1] = runif(n, temp_loc[ , 1] - x_window, temp_loc[ , 1] + x_window)

BALIVE[ , 2] = runif(n, temp_loc[ , 2] - y_window, temp_loc[ , 2] + y_window)

## years for basal area data

BALIVE[ , 3] = BALIVE_loc_year[ , 2]

## simulating months within years for basal area data

BALIVE[ , 4] = sample(1:12, n, replace = TRUE)

##############################################################################################

count_data = matrix(-99,  T,  S)

#for (i in 1:5)   count_data[c(i, i+5), loc_pos_sim[(cum_S_sim_year[i]+1):cum_S_sim_year[i+1]] ] = -99999

for (i in 1:5)   count_data[c(i, i+5), loc_pos_sim_list[[i]] ] = -99999

pos_avail = which(count_data == -99999)

##############################################################################################

## Simulating y_1 data (zero and non-zero basal area data)

pred_mat_y1 = t(pred_mat_full[ , pos_avail])

mean_y1 = pred_mat_y1%*%beta_y1

y1_sim = mean_y1 + rnorm(n)

pos_0 = which(y1_sim <= 0)

pos_1 = which(y1_sim > 0) 

count_data[ pos_avail[pos_0] ] = 0

count_data[ pos_avail[pos_1] ] = 1

##############################################################################################

## the simulated basal area data will be stored in basal_area

basal_area = matrix(-99,  T,  S)

#### Storing zero basal areas 

basal_area[ pos_avail[pos_0] ] = 0

##############################################################################################

## Simulating y_2 data

available_locations_all = list()

for (t in 1:T) available_locations_all[[t]] =   sort(which(count_data[t, ] == 1)) 

available_locations_unique = sort(unique(unlist(available_locations_all)) ) ## this stores all the cell locations, where any nonzero data was available during 10 years

S_A =  length(available_locations_unique)

pred_mat_3D = pred_mat_full_3D[ , , available_locations_unique]

mean_y2 = matrix(-99, T, S_A)

mean_y2[1, ] = c( t(pred_mat_3D[ , 1,  ]) %*% beta_y2 )

z_mat = matrix(-99, T, S_A)

z_mat[1, ] = 1

for (t in  2:T)	
{
  mean_y2_temp =  c( t( pred_mat_3D[ , t,  ]) %*% beta_y2 )
  
  z_mat[t, ] = ceiling( c(pnorm( mean_y2[(t-1), ] - mean_y2_temp ) ) - runif(S_A) ) 
  
  mean_y2[t, ] = z_mat[t, ]*mean_y2_temp + (1 - z_mat[t, ])*mean_y2[(t-1), ]
  
}

pos_avail_y2 = which( count_data[ , available_locations_unique] == 1)

y2_sim = mean_y2[pos_avail_y2] + sqrt(sigma2_y2)*rnorm(length(pos_avail_y2))

##############################################################################################

#### Storing non-zero basal areas 

basal_area[ pos_avail[pos_1] ] = exp(y2_sim)

basal_area = t(basal_area)

##############################################################################################

loc_pos_sim_list_all = c(loc_pos_sim_list, loc_pos_sim_list) 

basal_area_avail_list = list()

for(i in 1:T) basal_area_avail_list[[i]] = basal_area[loc_pos_sim_list_all[[i]], i]

BALIVE[ , 5] = unname(unlist(basal_area_avail_list))

##############################################################################################

write.csv(BALIVE, "Inputs_and_Outputs/hypothetical_basal_area_data.csv", row.names = FALSE)

