CDVineAIC <-
function(data, family, par, par2=rep(0,dim(data)[2]*(dim(data)[2]-1)/2), type)
{
  if(is.vector(data)){
    data = t(as.matrix(data))
  }else{
    data=as.matrix(data)
  }
  if(type == "CVine") type = 1
  else if(type == "DVine") type = 2
  if(type != 1 & type != 2) stop("Vine model not implemented.")    
    
  d=dim(data)[2]
  T=dim(data)[1]

	# Sicherheitsabfragen
	if(d<2) stop("Dimension has to be at least 2.")

	if(length(family) != d*(d-1)/2) stop("Number of copula families incorrect.")
	if(length(par) != d*(d-1)/2) stop("Number of copula parameters incorrect.")
	if(length(par2) != d*(d-1)/2) stop("Number of second copula parameters incorrect.")

	for(i in 1:(d*(d-1)/2))
	{
	if(!(family[i] %in% c(0,1,2,3,4,5,6,7,9,13,14,16,23,24,26,33,34,36))) stop("Copula family not implemented.")
	# Parameterbereiche abfragen
  	if((family[i]==1 || family[i]==2) && abs(par[i])>=1) stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
	if(family[i]==2 && par2[i]<=1) stop("The degrees of freedom parameter of the t-copula has to be larger than 1.")
	if((family[i]==3 || family[i]==13) && par[i]<=0) stop("The parameter of the Clayton copula has to be positive.")
  	if((family[i]==4 || family[i]==14) && par[i]<1) stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
  	if((family[i]==6 || family[i]==16) && par[i]<=1) stop("The parameter of the Joe copula has to be in the interval (1,oo).")	
  	if(family[i]==5 && par[i]==0) stop("The parameter of the Frank copula has to be unequal to 0.")
  	if(family[i]==7 && par[i]<=0) stop("The first parameter of the BB1 copula has to be positive.")
  	if(family[i]==7 && par2[i]<1) stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
  	if(family[i]==9 && par[i]<1) stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
  	if(family[i]==9 && par2[i]<=0) stop("The second parameter of the BB7 copula has to be positive.")
  	if((family[i]==23 || family[i]==33) && par[i]>=0) stop("The parameter of the rotated Clayton copula has to be negative.")
  	if((family[i]==24 || family[i]==34) && par[i]>-1) stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
  	if((family[i]==26 || family[i]==36) && par[i]>=-1) stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
	}
  
  npar = sum(family >= 1, na.rm=TRUE) + sum(family == 2,na.rm=TRUE) + sum(family == 7,na.rm=TRUE) + sum(family == 8,na.rm=TRUE) + sum(family == 9,na.rm=TRUE)
  npar_pair = (family>=1)+(family%in%c(2,7,8,9))

  like = CDVineLogLik(data, family, par, par2, type)

  AIC = -2*like$loglik + 2*npar
  pair.AIC = -2*like$ll + 2*npar_pair
  
  return(list(AIC=AIC,pair.AIC=pair.AIC))
}

CDVineBIC <-
function(data, family, par, par2=rep(0,dim(data)[2]*(dim(data)[2]-1)/2), type)
{
  if(is.vector(data)){
    data = t(as.matrix(data))
  }else{
    data=as.matrix(data)
  }
  if(type == "CVine") type = 1
  else if(type == "DVine") type = 2
  if(type != 1 & type != 2) stop("Vine model not implemented.")  
    
  d=dim(data)[2]
  T=dim(data)[1]

	# Sicherheitsabfragen
	if(d<2) stop("Dimension has to be at least 2.")

	if(length(family) != d*(d-1)/2) stop("Number of copula families incorrect.")
	if(length(par) != d*(d-1)/2) stop("Number of copula parameters incorrect.")
	if(length(par2) != d*(d-1)/2) stop("Number of second copula parameters incorrect.")

	for(i in 1:(d*(d-1)/2))
	{
	if(!(family[i] %in% c(0,1,2,3,4,5,6,7,9,13,14,16,23,24,26,33,34,36))) stop("Copula family not implemented.")
	# Parameterbereiche abfragen
  	if((family[i]==1 || family[i]==2) && abs(par[i])>=1) stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
	if(family[i]==2 && par2[i]<=1) stop("The degrees of freedom parameter of the t-copula has to be larger than 1.")
	if((family[i]==3 || family[i]==13) && par[i]<=0) stop("The parameter of the Clayton copula has to be positive.")
  	if((family[i]==4 || family[i]==14) && par[i]<1) stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
  	if((family[i]==6 || family[i]==16) && par[i]<=1) stop("The parameter of the Joe copula has to be in the interval (1,oo).")	
  	if(family[i]==5 && par[i]==0) stop("The parameter of the Frank copula has to be unequal to 0.")
  	if(family[i]==7 && par[i]<=0) stop("The first parameter of the BB1 copula has to be positive.")
  	if(family[i]==7 && par2[i]<1) stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
  	if(family[i]==9 && par[i]<1) stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
  	if(family[i]==9 && par2[i]<=0) stop("The second parameter of the BB7 copula has to be positive.")
  	if((family[i]==23 || family[i]==33) && par[i]>=0) stop("The parameter of the rotated Clayton copula has to be negative.")
  	if((family[i]==24 || family[i]==34) && par[i]>-1) stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
  	if((family[i]==26 || family[i]==36) && par[i]>=-1) stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
	}
  
  npar = sum(family >= 1, na.rm=TRUE) + sum(family == 2,na.rm=TRUE) + sum(family == 7,na.rm=TRUE) + sum(family == 8,na.rm=TRUE) + sum(family == 9,na.rm=TRUE)
  npar_pair = (family>=1)+(family%in%c(2,7,8,9))

  like = CDVineLogLik(data, family, par, par2, type)

  BIC = -2*like$loglik + log(T)*npar
  pair.BIC = -2*like$ll + log(T)*npar_pair
  
  return(list(BIC=BIC,pair.BIC=pair.BIC))
}

