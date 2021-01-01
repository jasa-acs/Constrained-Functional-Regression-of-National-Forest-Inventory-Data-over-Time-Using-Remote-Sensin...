################################################################################################################
################## 		Necessary Functions 		 ###########################################################
################################################################################################################

f_mean = function(x)  ifelse(max(abs(x))>0, mean(x[x !=0]),  0)	## to calculate the mean of non-zero values in a vector

## the coarse_data_creation function merges the 16x16 adjacent pixels
## input of this function:
## data: the list pixel_wise_TC_data of the data TC_features_with_missing_values.Rdata
## c: total number of pixels in X coordinate (Longitude) in the data
## d: total number of pixels in Y coordinate (Latitude) in the data
## x_merge: number of adjacent pixels to be merged in X axis 
## y_merge: number of adjacent pixels to be merged in Y axis 
 
coarse_data_creation = function(data, c, d, x_merge, y_merge)
{
  location = data$location
  sel_data = data$TC_features
  
  coarse_location = matrix(0, nrow = 1, ncol = 2)
  colnames(coarse_location) = c("X","Y")
  coarse_data = array(0,dim = c(1,120,3))

  for(j in 1:(d/y_merge))
    {
      y_location = which( (location[,2] <= max(location[,2]) - (j-1)*y_merge*30) & (location[,2] > max(location[,2]) - j*y_merge*30) )
      for(i in 1:(c/x_merge))
        {
          x_location = which( (location[,1] >= min(location[,1]) + (i-1)*x_merge*30) & (location[,1] < min(location[,1]) + i*x_merge*30) )
          position = intersect(x_location, y_location)
          coarse_location = rbind(coarse_location, c(mean(location[position,1]), mean(location[position,2])))
          coarse_data = abind(coarse_data, apply(sel_data[position,,], MARGIN = c(2,3), f_mean), along = 1)
        }
    }
    
	coarse_location = coarse_location[-1,]
	coarse_data = coarse_data[-1,,]
  
	coarse_output = list(coarse_location = coarse_location, coarse_data = coarse_data)
  	
	return(coarse_output)
}


## The input (A) of this function is a matrix and the output is a list of two components
## The first component is Cholesky decomposition of A and second component is the inverse of A
	
inv_and_logdet.sym = function(A)
{
  B = chol(A)
  result = list()
  result[[1]] = B
  result[[2]] = chol2inv(B)
  return(result)
}

## function to find the neighbors of any cell

f_W = function(b)	## the input of this function is any cell number from 1 to S and the output is the neighbors of that cell
{
out = which(rdist(matrix(location[b,],nrow=1), location) < 0.95*sqrt(x_dist^2 + y_dist^2) )

## Here x_dist and y_dist represents the distance in X axis and Y axis between two neighbors. If two cells share common edges, then they are called neighbor.
## Hence, the distance between the center of two neighbors should be less than sqrt(x_dist^2 + y_dist^2). So, we multiply sqrt(x_dist^2 + y_dist^2) by a number (0.95) close to 1.

out = setdiff(out,b)
return(out)
}

create_AR_cov_mat = function(kappa2, eta, total_months)	## Covariance matrix for AR(1) model
{
cov = c(kappa2)*tacvfARMA(phi = eta, maxLag = total_months-1)
covmat = toeplitz(cov)
return(covmat)
}

rtrunc_norm = function(number=1, a, b, mu, sigma2)	## random number generation from truncated normal distribution
  {
    F_a = pnorm (a, mean = mu, sd = sqrt(sigma2))
    F_b = pnorm (b, mean = mu, sd = sqrt(sigma2))
    x = qnorm((F_a + runif(number)*(F_b - F_a)), mean = mu, sd = sqrt(sigma2))
        return(x)
  }

phi_generate_Model_IIIB = function(input) ## generate phi and tau2

{
    phi_construct = phi
 
    mu_temp =  y_minus_X_beta/sigma2
	
    mu_temp[pos_0] = 0
 
    for (j in 1:S)
 
        {
            phi_var = 1/(1/sigma2 + W_plus[j]/tau2)

            phi_var[pos_0_cell[[j]]] =  tau2[pos_0_cell[[j]]]/W_plus[j]
 
            phi_construct[j,] = phi_var*(mu_temp[j,] + colSums(phi_construct[W[[j]], ])/tau2) + sqrt(phi_var) * rnorm(T)
						
        }

	phi_construct = sweep(phi_construct, 2, colMeans(phi_construct))

	tau2_post_shape = a20 + (S-1)*T/2  
	
	tau2_post_rate = b20 + sum((phi_construct[rep(1:S, W_plus),] - phi_construct[unlist(W),])^2)/4

	tau2_construct = 1/rgamma(1, shape = tau2_post_shape, rate = tau2_post_rate)
	
	output_spatial = rbind(phi_construct, rep(tau2_construct, T))
	
	return(output_spatial)	
}

phi_generate_Model_IV = function(input)	## generate phi and tau2 for Model IV 

{
    phi_construct = phi
 
    mu_temp =  sweep(y_minus_X_beta, 2, psi)/sigma2
	
    mu_temp[pos_0] = 0
 
    for (j in 1:S)
 
        {
            phi_var = 1/(1/sigma2 + W_plus[j]/tau2)

            phi_var[pos_0_cell[[j]]] =  tau2[pos_0_cell[[j]]]/W_plus[j]
 
            phi_construct[j,] = phi_var*(mu_temp[j,] + colSums(phi_construct[W[[j]], ])/tau2) + sqrt(phi_var) * rnorm(T)
						
        }

	phi_construct = sweep(phi_construct, 2, colMeans(phi_construct))

	tau2_post_shape = a20 + (S-1)*T/2  
	
	tau2_post_rate = b20 + sum((phi_construct[rep(1:S, W_plus),] - phi_construct[unlist(W),])^2)/4

	tau2_construct = 1/rgamma(1, shape = tau2_post_shape, rate = tau2_post_rate)
	
	output_spatial = rbind(phi_construct, rep(tau2_construct, T))
	
	return(output_spatial)	
}


## to plot Figure 2 (b)

spatial_plot_Figure_2b = function(locations,dataset,log_scale=F,grey_image=F,titleplot="Spatial Plot",layout_value = c(1,1,1),subplot_name="Plot")
{
    
    spatial_data= SpatialPixelsDataFrame(locations, data.frame(dataset), tolerance=0.01,proj4string=CRS("+proj=longlat +datum=WGS84"))
    
    n_point=6
    n_color = 1000
	
    if (grey_image==TRUE) g.col<-colorRampPalette((c(grey(.999),grey(.90),grey(.8),grey(.7), grey(.5),grey(.0001))))  #colors scheme - grayscale
    if (grey_image==FALSE) g.col<-colorRampPalette(rev(c("blue", "forestgreen","lightgreen", "yellow", "darkgoldenrod","black")))  #colors scheme - color

    if(log_scale==TRUE)
    {
        spatial_data@data[spatial_data@data<=0]=min(spatial_data@data[spatial_data@data>0])
        spatial_data@data=log(spatial_data@data)
        zlim=range(spatial_data@data,na.rm=T)
        at=seq(zlim[1],zlim[2],length.out=n_point)  #log-scaled color key values
	c=seq((min(spatial_data@data)),(max(spatial_data@data)),length.out=n_point)
        l=list(at=c,cex=1,labels=format(round((exp(c) -1),3),format="scientific"))
    }
	   
    if(log_scale==FALSE)
    {
        zlim=range(spatial_data@data,na.rm=T)
        at=seq(zlim[1],zlim[2],length.out=n_point)
        l=list(at=at,cex=1,labels=format(round(at,2),format="scientific"))
    }
    
    p1= spplot(spatial_data,col.regions =g.col(n_color),
    at=seq(min(spatial_data@data),max(spatial_data@data),length.out=n_color), 
    par.strip.text=list(cex=1,font=3),
    par.settings=list(strip.background=list(col="transparent")),aspect=1,
    scales=list(draw=F),
    colorkey=list(height=1,col=g.col(n_color),labels=l),
    sub="",xlab="",ylab="", 
    layout=layout_value,   
    names.attr=subplot_name, main=list(label=titleplot,cex=1))
    
    return(p1)
}

## this function is used to plot TC features

spatial_plot = function(locations,dataset,log_scale=F,grey_image=F,titleplot="Spatial Plot",layout_value = c(1,1,1),subplot_name="Plot")
{
    
    spatial_data= SpatialPixelsDataFrame(locations, data.frame(dataset), tolerance=0.01,proj4string=CRS("+proj=longlat +datum=WGS84"))
    
	n_point = 7
	n_color = 1000
	
	if (grey_image==TRUE) g.col<-colorRampPalette(rev(c(grey(.98),grey(.90),grey(.8),grey(.7),grey(.0001))))  #colors scheme - grayscale
	if (grey_image==FALSE)
		{	
			if (TC_number == 1) g.col<-colorRampPalette(hsv(1/6,1,seq(0.15,1, by = .001),1)) #colors scheme - color
			if (TC_number == 2) g.col<-colorRampPalette(hsv(1/3,.755,seq(0,1, by = .001),1)) #colors scheme - color
			if (TC_number == 3) g.col<-colorRampPalette(c(hsv(209/360,.7,seq(0.15,1, by = .05),1), rev(c("#bee7ff", "#83d0ff", "#5cc2ff", "#34b3ff"))))   #colors scheme - color
		}

    if(log_scale==TRUE)
    {
        spatial_data@data[spatial_data@data<=0]=min(spatial_data@data[spatial_data@data>0])
        spatial_data@data=log(spatial_data@data)
        zlim=range(spatial_data@data,na.rm=T)
        at=seq(zlim[1],zlim[2],length.out=n_point)  #log-scaled color key values
	c=seq((min(spatial_data@data)),(max(spatial_data@data)),length.out=n_point)
        l=list(at=c,cex=1,labels=format(round((exp(c) -1),3),format="scientific"))
    }
	   
    if(log_scale==FALSE)
	{
        zlim=range(spatial_data@data,na.rm=T)
 
	if (TC_number == 1) 
		{
			at=seq(zlim[1],zlim[2],length.out=n_point)
			l=list(at=at,cex=2,labels=format(round(exp(at),2),format="scientific"))
		}
	if (TC_number == 2 ) 
		{
			at= union(seq(zlim[1], 0 ,length.out = floor(n_point/2) + 1 ), seq(0, zlim[2],length.out = floor(n_point/2) + 1))
			l=list(at=at,cex=2,labels=format(round(at, 2)),format="scientific")
		}
		
	if ( TC_number == 3) 
		{
			at= union(c(zlim[1], 0 ), seq(0, zlim[2],length.out = 6))
			l=list(at=at,cex=2,labels=format(round(at, 2)),format="scientific")
		}
	}
    
    p1= spplot(spatial_data,col.regions =g.col(n_color), cex = 4, font = 4,
    at=seq(min(spatial_data@data),max(spatial_data@data),length.out=n_color), 
    par.strip.text=list(cex=2,font=4),
    par.settings=list(strip.background=list(col="transparent")),aspect=1,
    scales=list(draw=F),#alternating=1),
    colorkey=list(height=1,col=g.col(n_color),labels=l),
    sub="",xlab="",ylab="",
    layout=layout_value,   #uncomment this to fix rows and columns in display
    names.attr=subplot_name, main=list(label=titleplot,cex=4))
    
    return(p1)
}
