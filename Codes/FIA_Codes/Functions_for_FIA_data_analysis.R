################################################################################################################
################## 		Necessary Functions 		 ###########################################################
################################################################################################################

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

Bt_vector_create = function(z, t)
	{
		B =which( c( rev(1 - cummax( rev(z[2:t]) ) ), 1) * z[1:t] == 1)
		return(B)
	}

Bt_matrix_create = function(z)
	{
		B = c(1)
#		for(t in 1) B[t] = 1
		for(t in 2:T) B[t] = Bt_vector_create(z, t)
		return(B)
	}

rmult_new_multiple = function(prob_mat) {  rowSums(ceiling(runif(nrow(prob_mat)) - rowCumsums(prob_mat))) + 1  }

rtrunc_norm = function(number=1, a, b, mu, sigma2) 
  {
    F_a = pnorm (a, mean = mu, sd = sqrt(sigma2))
    F_b = pnorm (b, mean = mu, sd = sqrt(sigma2))
    x = qnorm((F_a + runif(number)*(F_b - F_a)), mean = mu, sd = sqrt(sigma2))
        return(x)
  }

position_covariate = function(a, coarse_location)
{
   pos_x =  which((a[1] >= coarse_location[,1] - 30*(x_merge/2)) & a[1] < coarse_location[,1] + 30*(x_merge/2))
   pos_y =  which((a[2] >= coarse_location[pos_x,2] - 30*(y_merge/2)) & a[2] < coarse_location[pos_x,2] + 30*(y_merge/2))

   pos = pos_x[pos_y]
   return(pos) 
}

## to plot Figure 7

spatial_plot = function(locations,dataset,log_scale=F,grey_image=F,titleplot="",layout_value = c(1,1,1),subplot_name="Plot", constant = 1)
{
	min_data = min(dataset, na.rm=T)
	dataset = dataset - min_data + constant
	spatial_data= SpatialPixelsDataFrame(locations, data.frame(dataset), tolerance=0.01,proj4string=CRS("+proj=longlat +datum=WGS84"))
    
	n_point=6
	n_color = 1000
	
    if (grey_image==TRUE) g.col<-colorRampPalette(rev(c(grey(.98),grey(.90),grey(.8),grey(.7),grey(.0001))))  #colors scheme - grayscale
    if (grey_image==FALSE) g.col<-colorRampPalette(rev(c("blue", "forestgreen", "lightgreen", "yellow", "darkgoldenrod","black")))  
 

 if(log_scale==TRUE)
    {
        spatial_data@data=log(spatial_data@data)
        zlim=range(spatial_data@data,na.rm=T)
        at=seq(zlim[1],zlim[2],length.out=n_point)  #log-scaled color key values
	c=seq((min(spatial_data@data)),(max(spatial_data@data)),length.out=n_point)
        l=list(at=c,cex=2,labels=format(round((exp(c) + min_data - constant), 2),format="scientific"))
    }
	   
    if(log_scale==FALSE)
    {
        zlim=range(spatial_data@data,na.rm=T)
        at=c(seq(zlim[1],zlim[2],length.out=n_point), (constant - min_data))
        l=list(at=at,cex=2,labels=format(round((at + min_data - constant),2 ),format="scientific"))
    }
    
    p1= spplot(spatial_data,col.regions =g.col(n_color),
    at=seq(min(spatial_data@data),max(spatial_data@data),length.out=n_color), 
    par.strip.text=list(cex=2,font=4),
    par.settings=list(strip.background=list(col="transparent")),aspect=1,
    scales=list(draw=F),
    colorkey=list(height=1,col=g.col(n_color),labels=l),
    sub="",xlab="",ylab="",		
    layout=layout_value,   
    names.attr=subplot_name, main=list(label=titleplot,cex=1))
    
    return(p1)
}

## to plot Figure 8 and 9(b)

spatial_plot_uncertainty = function(locations,dataset,log_scale=F,grey_image=F,titleplot="",layout_value = c(1,1,1),subplot_name="Plot", constant = 1)
{
	min_data = min(dataset, na.rm=T)
	dataset = dataset - min_data + constant
	spatial_data= SpatialPixelsDataFrame(locations, data.frame(dataset), tolerance=0.01,proj4string=CRS("+proj=longlat +datum=WGS84"))
    
	n_point=6
	n_color = 1000
	
    if (grey_image==TRUE) g.col<-colorRampPalette(rev(c(grey(.98),grey(.90),grey(.8),grey(.7),grey(.0001))))  #colors scheme - grayscale
    if (grey_image==FALSE) g.col<-colorRampPalette(rev(c("blueviolet", "orangered", "orange", "yellow", "darkgoldenrod","black")))  
	
	
    if(log_scale==TRUE)
    {
        spatial_data@data=log(spatial_data@data)
        zlim=range(spatial_data@data,na.rm=T)
        at=seq(zlim[1],zlim[2],length.out=n_point)  #log-scaled color key values
	c=seq((min(spatial_data@data)),(max(spatial_data@data)),length.out=n_point)
        l=list(at=c,cex=2,labels=format(round((exp(c) + min_data - constant), 2),format="scientific"))
    }
	   
    if(log_scale==FALSE)
    {
        zlim=range(spatial_data@data,na.rm=T)
        at=c(seq(zlim[1],zlim[2],length.out=n_point), (constant - min_data))
        l=list(at=at,cex=2,labels=format(round((at + min_data - constant),2 ),format="scientific"))
    }
    
    p1= spplot(spatial_data,col.regions =g.col(n_color),
    at=seq(min(spatial_data@data),max(spatial_data@data),length.out=n_color), 
    par.strip.text=list(cex=2,font=3),
    par.settings=list(strip.background=list(col="transparent")),aspect=1,
    scales=list(draw=F),
    colorkey=list(height=1,col=g.col(n_color),labels=l),
    sub="",xlab="",ylab="",		
    layout=layout_value,   
    names.attr=subplot_name, main=list(label=titleplot,cex=1))
    
    return(p1)
}

## to plot Figure 9(a)

spatial_plot_five_year_change = function(locations,dataset,log_scale=F,grey_image=F,titleplot="",layout_value = c(1,1,1),subplot_name="Plot")
{
	min_data = min(dataset, na.rm=T)
	dataset = dataset - min_data + 1
	spatial_data= SpatialPixelsDataFrame(locations, data.frame(dataset), tolerance=0.01,proj4string=CRS("+proj=longlat +datum=WGS84"))
    
	n_point=6
	n_color = 1000
	
    if (grey_image==TRUE) g.col<-colorRampPalette(rev(c(grey(.98),grey(.90),grey(.8),grey(.7),grey(.0001))))  #colors scheme - grayscale
    if (grey_image==FALSE) g.col_1 <- colorRampPalette(rev(c("yellow", "darkgoldenrod","black")))  #colors scheme - color
    if (grey_image==FALSE) g.col_2 <- colorRampPalette(rev(c("blue", "forestgreen","lightgreen")))  #colors scheme - color

    if(log_scale==TRUE)
    {
        spatial_data@data=log(spatial_data@data)
        zlim=range(spatial_data@data,na.rm=T)
        at=seq(zlim[1],zlim[2],length.out=n_point)  #log-scaled color key values
	c=seq((min(spatial_data@data)),(max(spatial_data@data)),length.out=n_point)
        l=list(at=c,cex=2,labels=format(round((exp(c) + min_data -1), 2),format="scientific"))
    }
	   
    if(log_scale==FALSE)
    {
        zlim=range(spatial_data@data,na.rm=T)
        at= unique( c(seq(zlim[1], (1 - min_data),length.out = (n_point/2) + 1), seq((1 - min_data), zlim[2], length.out = n_point/2 )) )
        l=list(at=at,cex=2,labels=format(round((at + min_data -1),2 ),format="scientific"))
    }
    
    p1= spplot(spatial_data,col.regions = c( g.col_1(n_color), g.col_2(n_color)) ,
    at= unique( c( seq(min(spatial_data@data), (1 - min_data), length.out = n_color ), seq( (1 - min_data), max(spatial_data@data), length.out = n_color ) )) ,
    par.strip.text=list(cex=2,font=3),
    par.settings=list(strip.background=list(col="transparent")),aspect=1,
    scales=list(draw=F),
    colorkey=list(height=1,col= c( g.col_1(n_color), g.col_2(n_color)),labels=l),
    sub="",xlab="",ylab="",
    layout=layout_value,  
    names.attr=subplot_name, main=list(label=titleplot,cex=1))
    
    return(p1)
}

