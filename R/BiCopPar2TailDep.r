BiCopPar2TailDep<-function(family,par,par2=0)
{
	if(!(family %in% c(0,1,2,3,4,5,6,7,9,13,14,16,23,24,26,33,34,36))) stop("Copula family not implemented.")
	if((family==1 || family==2) && abs(par[1])>=1) stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
  if(family==2 && par2<=1) stop("The degrees of freedom parameter of the t-copula has to be larger than 1.")	
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
  
  if(family==0 | family==1 | family==5 | family%in%c(23,24,26,33,34,36))
	{
		lower=0
		upper=0
	}
	else if(family==2)
	{
		lower=2*pt((-sqrt(par2+1)*sqrt((1-par)/(1+par))),df=par2+1)
		upper=lower
	}
	else if(family==3)
	{
		lower=2^(-1/par)
		upper=0
	}
	else if(family==4 | family==6)
	{
		lower=0
		upper=2-2^(1/par)
	}
	else if(family==7)
	{
    lower=2^(-1/(par*par2))
    upper=2-2^(1/par2)
	}
	else if(family==9)
	{
    lower=2^(-1/par2)
    upper=2-2^(1/par)
	}
	else if(family==13)
	{
		lower=0
		upper=2^(-1/par)
	}
	else if(family==14 | family==16)
	{
		lower=2-2^(1/par)
		upper=0
	}

  return(list(lower=lower,upper=upper))
}