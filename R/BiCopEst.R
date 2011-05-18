BiCopEst <-
function(u1,u2,family, method="mle", se=FALSE, max.df=30, max.BB=list(BB1=c(5,6),BB7=c(5,6)))
{
  # Function that estimates the parameter(s) of the bivatiate copula
  #---------------------------------------------------------
  # INPUT:
  #   u1,u2      Data for which to estimate parameter
  #   family            The array definig the copulas in the pcc copula construction
  # OUTPUT:
  #   theta      Estimated Parameters
  #----------------------------------------------------------
  # Author: Carlos Almeida  <calmeida at ma.tum.de>
  # Update: Ulf Schepsmeier <schepsmeier at ma.tum.de>
  # Date: 2008-12-08
  # Update date: 2010-12-27
  # Version: 1.1
  #---------------------------------------------------------------

  # Sicherheitsabfragen
  if(is.null(u1)==TRUE || is.null(u2)==TRUE) stop("u1 and/or u2 are not set or have length zero.")
  if(length(u1) != length(u2)) stop("Lengths of 'u1' and 'u2' do not match.")
  if(length(u1)<2) stop("Number of observations has to be at least 2.")
  if(!(family %in% c(0,1,2,3,4,5,6,7,9,13,14,16,23,24,26,33,34,36))) stop("Copula family not implemented.")

  if(max.df<=1) stop("The upper bound for the degrees of freedom parameter has to be larger than 1.")
  if(!is.list(max.BB)) stop("'max.BB' has to be a list.")
  if(max.BB$BB1[1] < 0.001) stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
  if(max.BB$BB1[2] < 1.001) stop("The upper bound for the second parameter of the BB1 copula should be greater than 1.001 (lower bound for estimation).")
  if(max.BB$BB7[1] < 1.001) stop("The upper bound for the first parameter of the BB7 copula should be greater than 1.001 (lower bound for estimation).")
  if(max.BB$BB7[2] < 0.001) stop("The upper bound for the second parameter of the BB7 copula should be greater than 0.001 (lower bound for estimation).")


  if(method!="mle" && method!="itau") stop("Estimation method has to be either 'mle' or 'itau'.")

  if(is.logical(se)==FALSE) stop("'se' has to be a logical variable (TRUE or FALSE).")
  
  if(method=="itau" && family %in% c(2,7,8,9))
  {
  message("For two parameter copulas the estimation method 'itau' cannot be used. The method is automatically set to 'mle'.")
  method="mle"
  }

  if(family!=0)
  {
	#tau <- cor(u1,u2,method="kendall")
	tau <- fasttau(u1,u2)
  }
  
    theta=0
    if(family==0) # independent
    {
	theta=0;
    }
    else if(family==1) ## Gaussian 
    {
        theta <- sin(tau*pi/2)
    }
    else if(family==3 || family==13) ## Clayton
    {
	if(tau<=0) stop("Clayton copula cannot be used for negatively dependent data.")
	theta <- max(0,2*tau/(1-tau))
    }
    else if(family==4 || family==14) ## Gumbel
    {
	if(tau<0) stop("Gumbel copula cannot be used for negatively dependent data.")
	theta <- max(1,1/(1-tau))
    }
    else if(family==5) ## Frank
    {
	theta=Frank.itau.JJ(tau)
    }
    else if(family==6 || family==16) ## Joe
    {
	if(tau<=0) stop("Joe copula cannot be used for negatively dependent data.")
  theta=Joe.itau.JJ(tau)	
    }
    else if(family==23 || family==33)
    {
	if(tau>=0) stop("Rotated Clayton copula cannot be used for positively dependent data.")
	theta <- (2*tau/(1+tau))
    }
    else if(family==24 || family==34)
    {
	if(tau>0) stop("Rotated Gumbel copula cannot be used for positively dependent data.")
	theta <- -(1/(1+tau))
    }
    else if(family==26 || family==36)
    {
	if(tau>=0) stop("Rotated Joe copula cannot be used for positively dependent data.")
	theta=-Joe.itau.JJ(-tau)
    }

    se1=0
    if(method=="itau" && se==TRUE)
    {
    p = 2
  	n = length(u1)
  	ec = numeric(n)
  	u = cbind(u1,u2)
  	v = matrix(0,n,p*(p-1)/2)

	if(family == 1) tauder = function(x) 2/(pi * sqrt(1-x^2))
	else if(family %in% c(3,13,23,33)) tauder = function(x) 2*(2+x)^(-2)
	else if(family %in% c(4,14,24,34)) tauder=function(x) x^(-2)
	else if(family == 5){
  	  tauder=function(x){  		
  		 f = function(x) x/(exp(x) - 1)
  		 4/x^2 - 8/x^3 * integrate(f, lower = 0+.Machine$double.eps^0.5, upper = x)$value +4/(x*(exp(x)-1))  		
  	  }
	}  	    
	else if(family %in% c(6,16,26,36)){
  	  tauder=function(x){
            euler=0.5772156649015328606
            -((-2+2*euler+2*log(2)+digamma(1/x)+digamma(1/2*(2+x)/x)+x)/(-2+x)^2)+((-trigamma(1/x)/x^2+trigamma(1/2*(2+x)/x)*(1/(2+x)-(2+x)/(2*x^2))+1)/(-2+x))
        }  
	}
	
	l <- 1
  	for (j in 1:(p-1)) {
  		for (i in (j+1):p) {
  			for (k in 1:n) 
  			ec[k] <- sum(u[,i] <= u[k,i] & u[,j] <= u[k,j])/n 
  			v[,l] <- 2 * ec - u[,i] - u[,j]
  			l <- l + 1
  		}
  	}
  	
	if(family == 0)
 		D = 0
	else if(family %in% c(1,3,13,4,14,5,6,16))
  		D = 1 / tauder(theta)
	else if(family %in% c(23,33,24,34,26,36))
  		D = 1 / tauder(-theta)    
	
    
	se1 = as.numeric(sqrt(16/n * var(v %*% D)))
    }


	if(method=="mle")
	{
	theta1=0
	delta=0

	if(family!=2 && family!=7 && family!=8 && family!=9)
	{
		theta1=theta
	}
	if(family==2) ## t 
	{
		theta1 <- sin(tau*pi/2)
		delta <- (max.df-1)/2       # Nehme die Mitte zwischen 1 und max.df So kann man mit dem Startwert auch nicht außerhalb des vom User gesetzten Bereiches sein.
	}
	else if(family==7)	## BB1
	{
		if(tau<0)
		{
		print("The BB1 copula cannot be used for negatively dependent data.")
		delta=1.001
		theta1=0.001
		}
		else
		{
		#delta <- 1.5
		#theta1=0.5	
    delta=(max.BB$BB1[2]-1.001)/2
    theta1=(max.BB$BB1[1]-0.001)/2
		}
	}
	else if(family==8)	## BB2
	{
		if(tau<0) 
		{
		print("The BB2 copula cannot be used for negatively dependent data.")
		delta=1.001
		theta1=0.001
		}
		else
		{
		#delta <- 1.5
		#theta1=0.5
    delta=(max.BB$BB2[2]-1.001)/2
    theta1=(max.BB$BB2[1]-0.001)/2
		}
	}
	else if(family==9)	## BB7
	{
		if(tau<0) 
		{
		print("The BB7 copula cannot be used for negatively dependent data.")
		delta=0.001
		theta=1.001
		}
		else
		{
		delta <- 0.5
		theta1=1.5
    delta=(max.BB$BB7[2]-0.001)/2
    theta1=(max.BB$BB7[1]-1.001)/2
		}
	}

	out=MLE_intern(cbind(u1,u2),c(theta1, delta),family=family,se,max.df,max.BB)
	theta=out$par
	if(se==TRUE)
		se1=out$se
	}
	

out2=list()
if(length(theta)==2)
  {
    out2$par=theta[1]
    out2$par2=theta[2]
  }
  else
  {
    out2$par=theta
    out2$par2=0
  }
if(se==TRUE)
{
  if(length(se1)==2)
  {
	  out2$se=se1[1]
    out2$se2=se1[2]
  }
  else
  {
    out2$se=se1
    out2$se2=0
  }
}

    return(out2)
}



Frank.itau.JJ<-function(tau){
a<-1
if(tau<0){
a<- -1
tau<- -tau}
f = function(x) {
x/(exp(x) - 1)}
tauF=function(x) 1 - 4/x + 4/x^2 * integrate(f, lower = 0+.Machine$double.eps^0.5, upper = x)$value
v<-uniroot(function(x) tau - tauF(x) ,lower=0,upper=100, tol = .Machine$double.eps^0.5)$root
return(a*v)
}



Joe.itau.JJ<-function(tau){
if(tau<0){
return(1.000001)
}else{
	
	tauF=function(a){
		euler=0.5772156649015328606
		1+((-2+2*euler+2*log(2)+digamma(1/a)+digamma(1/2*(2+a)/a)+a)/(-2+a))
	}
	
	
v<-uniroot(function(x) tau - tauF(x) ,lower=1,upper=100, tol = .Machine$double.eps^0.5)$root
return(v)}
}




#############################################################
# bivariate MLE function  
# 
#------------------------------------------------------------
# INPUT:
#   data	Data for which to estimate parameter
#   start.parm	Start parameter for the MLE
#   Maxiter	max number of iterations
#   se		TRUE or FALSE
# OUTPUT:
#   out     Estimated Parameters and standard error (if se==TRUE)
#--------------------------------------------------------------
# Author: Ulf Schepsmeier
# Date: 2011-02-04
# Version: 1.1
#---------------------------------------------------------------

MLE_intern <-
function(data,start.parm,family,se=FALSE,max.df=30,max.BB=list(BB1=c(5,6),BB7=c(5,6)),weights=NULL)
{

  n = dim(data)[1]

  if(family %in% c(7,8,9)){
			
	t_LL = function(param)
	{
								
  	if(is.null(weights))
  	{
               ll = .C("LL_mod",as.integer(family),as.integer(n),as.double(data[,2]),as.double(data[,1]),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
          }
  	else
  	{
               ll = .C("LL_mod_seperate",as.integer(family),as.integer(n),as.double(data[,2]), as.double(data[,1]), as.double(param[1]),as.double(param[2]), as.double(rep(0,n)),PACKAGE='CDVine')[[7]]%*%weights
          }
  	 
  	if(is.infinite(ll) || is.na(ll)  || ll< -10^300) ll = -10^300
  			
  	return(ll)
  	}
		
		if(family == 7){
      low = c(0.001,1.001)
      up = max.BB$BB1
		}else if(family == 8){
      low = c(0.001,0.001)
      up = c(5,5)
		}else if(family == 9){ 
      low = c(1.001,0.001)
      up = max.BB$BB7
    }
    
    if(se == TRUE){      
      optimout = optim(par=start.parm,fn=t_LL,method="L-BFGS-B",lower=low,control=list(fnscale=-1,maxit = 500),hessian=TRUE)
    }else{
      optimout = optim(par=start.parm,fn=t_LL,method="L-BFGS-B",lower=low,control=list(fnscale=-1,maxit = 500))
    }

  }else if(family == 2){
  
    t_LL = function(param){
								
			if(is.null(weights))
			{
			     ll = .C("LL_mod",as.integer(family),as.integer(n),as.double(data[,2]),as.double(data[,1]),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
			}
			else
			{
			     ll = .C("LL_mod_seperate",as.integer(family),as.integer(n),as.double(data[,2]), as.double(data[,1]), as.double(param[1]),as.double(param[2]), as.double(rep(0,n)),PACKAGE='CDVine')[[7]]%*%weights
			}
			
      if(is.infinite(ll) || is.na(ll)  || ll< -10^300) ll = -10^300
			
      return(ll)
		}
		
    if(se == TRUE){      
      optimout = optim(par=start.parm,fn=t_LL,method="L-BFGS-B",control=list(fnscale=-1,maxit = 500),hessian=TRUE,lower=c(-0.9999,1.0001),upper=c(0.9999,max.df))
    }else{
      optimout = optim(par=start.parm,fn=t_LL,method="L-BFGS-B",control=list(fnscale=-1,maxit = 500),lower=c(-0.9999,1.0001),upper=c(0.9999,max.df))
    }
    
    if(optimout$par[2] >= (max.df-0.0001)) warning(paste("Degrees of freedom of the t-copula estimated to be larger than ",max.df,". Consider using the Gaussian copula instead.",sep=""))
    						
  }else{
		
    t_LL = function(param){
		
		if(is.null(weights))
		{
			ll = .C("LL_mod",as.integer(family),as.integer(n),as.double(data[,2]), as.double(data[,1]), as.double(param),as.double(0), as.double(0),PACKAGE='CDVine')[[7]]
		}
		else
		{
			ll = .C("LL_mod_seperate",as.integer(family),as.integer(n),as.double(data[,2]), as.double(data[,1]), as.double(param[1]),as.double(param[2]), as.double(rep(0,n)),PACKAGE='CDVine')[[7]]%*%weights
		}				
      if(is.infinite(ll) || is.na(ll) || ll< -10^300) ll = -10^300
      
      return(ll)
    }    

    low = -Inf
    up = Inf

		if(family == 1){
      low = -0.9999
      up = 0.9999
    }else if(family %in% c(3,13)){
      low = 0.0001
    }else if(family %in% c(4,14)){
      low = 1.0001
    }else if(family %in% c(6,16)){
      low = 1.0001
    }else if(family %in% c(23,33)){
      up = -0.0001
    }else if(family %in% c(24,34)){
      up = -1.0001
    }else if(family %in% c(26,36)){
      up = -1.0001
    }   
    	    
    if(se == TRUE){      
      optimout = optim(par=start.parm[1],fn=t_LL,method="L-BFGS-B",control=list(fnscale=-1,maxit = 500),lower=low,upper=up,hessian=TRUE)
    }else{
      optimout = optim(par=start.parm[1],fn=t_LL,method="L-BFGS-B",control=list(fnscale=-1,maxit = 500),lower=low,upper=up)
    }

  } 
  
	out=list()
			
	if(se == TRUE){ 

    if(family %in% c(2,7,8,9)){
    	  
      out$par = optimout$par   
  
      if(det(optimout$hessian)==0){
  		  var = diag(1,dim(optimout$hessian)[1])
      }else{
  		  var = (-solve(optimout$hessian))
      } 
      
      out$se = sqrt(diag(var))  
    
      if(family == 2 && out$par[2] >= (max.df-0.0001)) out$se[2] = NA 
    
    }else{   

  	  out$par = optimout$par[1]
  
      if(optimout$hessian==0){
  		  var = 1
      }else{
  		  var = -1/optimout$hessian
      } 
 
      out$se = as.numeric(sqrt(var))
         
    }																							  
																										  		
  }else{
    
    if(family %in% c(2,7,8,9)){
    
      out$par = optimout$par
    
    }else{
    
      out$par[1] = optimout$par[1]
    
    }
  
  }

  return(out)
}


fasttau<- function(x, y)
{
  m=length(x)
  n=length(y)
  if(m == 0 || n == 0) stop("both 'x' and 'y' must be non-empty")
  if(m != n) stop("'x' and 'y' must have the same length")
  out <- .C("ktau",
		x=as.double(x), 
		y=as.double(y),
		N=as.integer(n),
		tau=as.double(0), 
		S=as.double(0),
		D=as.double(0),
		T=as.integer(0),
		U=as.integer(0),
		V=as.integer(0),
		PACKAGE='CDVine')
  ktau=out$tau

return(ktau)
}
