BiCopPar2Tau<-function(family,par,par2=0)
{
	if(!(family %in% c(0,1,2,3,4,5,6,7,9,13,14,16,23,24,26,33,34,36))) stop("Copula family not implemented.")
	if((family==1 || family==2) && abs(par[1])>=1) stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
  #no check for degrees of freedom parameter required!
	if((family==3 || family==13) && par<=0) stop("The parameter of the Clayton copula has to be positive.")
	if((family==4 || family==14) && par<1) stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
	if((family==6 || family==16) && par<=1) stop("The parameter of the Joe copula has to be in the interval (1,oo).")	
	if(family==5 && par==0) stop("The parameter of the Frank copula has to be unequal to 0.")
	if(family==7 && par<=0) stop("The first parameter of the BB1 copula has to be positive.")
	if(family==7 && par2<1) stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
	if(family==9 && par<1) stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
	if(family==9 && par2<=0) stop("The second parameter of the BB7 copula has to be positive.")
	if((family==23 || family==33) && par>=0) stop("The parameter of the rotated Clayton copula has to be negative.")
	if((family==24 || family==34) && par>-1) stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
	if((family==26 || family==36) && par>=-1) stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
  
  if(family==0)
	{
		tau=0
	}
	else if(family==1 | family==2)
	{
		tau=2/pi*asin(par)
	}
	else if(family==3 || family==13)
	{
		tau=par/(par+2)
	}
	else if(family==4 || family==14)
	{
		tau=1-1/par
	}
	else if(family==5)
	{
		f=function(x) {x/(exp(x)-1)}
		if(par>0)
		{
		tau=1-4/par+4/par^2*integrate(f,lower=0, upper=par)$value
		}
		else
		{
		tau=1-4/par-4/par^2*integrate(f,lower=par, upper=0)$value
		}
	}
	else if(family==6 || family==16)
	{
		if(par==2)
			tau=0.35
		else
		{
			euler=0.57721
			tau=1+((-2+2*euler+2*log(2)+digamma(1/par)+digamma(1/2*(2+par)/par)+par)/(-2+par))
		}
	}
	else if(family==7)
	{
		theta=par
		delta=par2
		tau=1-2/(delta*(theta+2))
	}
	else if(family==9)
	{
		theta=par
		delta=par2
		tau=1-2/(delta*(2-theta))+4/(theta^2*delta)*gamma(delta+2)*gamma((2-2*theta)/(theta)+1)/gamma(delta+3+(2-2*theta)/(theta))

	}
	else if(family==23 || family==33)
	{
		tau=par/(-par+2)
	}
	else if(family==24 || family==34)
	{
		tau=-1-1/par
	}
	else if(family==26 || family==36)
	{
		if(par==2)
			tau=0.35
		else
		{
			euler=0.57721
			tau=-1-((-2+2*euler+2*log(2)+digamma(1/-par)+digamma(1/2*(2-par)/(-par))-par)/(-2-par))
		}
	}

return(tau)
}