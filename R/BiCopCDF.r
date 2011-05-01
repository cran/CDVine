BiCopCDF<-function(u1,u2,family,par,par2=0)
{
  if(is.null(u1)==TRUE || is.null(u2)==TRUE) stop("u1 and/or u2 are not set or have length zero.")
	if(length(u1)!=length(u2)) stop("Lengths of 'u1' and 'u2' do not match.")
	if(!(family %in% c(0,1,2,3,4,5,6,7,9,13,14,16,23,24,26,33,34,36))) stop("Copula family not implemented.")
	if(c(2,7,8,9) %in% family && par2==0) stop("For t-, BB1 and BB7 copulas, 'par2' must be set.")
	if(c(1,3,4,5,6,13,14,16,23,24,26,33,34,36) %in% family && length(par)<1) stop("'par' not set.")
	
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

  res = rep(NA, length(u1))
  
  if(family == 0){
    res = u1*u2
  }else if(family == 1){
    cdf = function(u,v) pmvnorm(upper=c(qnorm(u),qnorm(v)), corr=matrix(c(1,par,par,1),2,2))
    res = mapply(cdf, u1, u2, SIMPLIFY=TRUE)
  }else if(family == 2){
    cdf = function(u,v) pmvt(upper=c(qt(u,df=par2),qt(v,df=par2)), corr=matrix(c(1,par,par,1),2,2), df=par2)
    res = mapply(cdf, u1, u2, SIMPLIFY=TRUE)   
  }else if(family %in% c(3,4,5,6,7,8,9)){
    cdf = function(u,v) .C("copCdf",as.double(u),as.double(v),as.integer(1),as.double(c(par,par2)),as.integer(family),as.double(0),PACKAGE='CDRVine')[[6]]
    res = mapply(cdf, u1, u2, SIMPLIFY=TRUE)  
  }else if(family %in% c(13,14,16)){
    cdf = function(u,v) u+v-1+.C("copCdf",as.double(1-u),as.double(1-v),as.integer(1),as.double(c(par,par2)),as.integer(family-10),as.double(0),PACKAGE='CDRVine')[[6]]
    res = mapply(cdf, u1, u2, SIMPLIFY=TRUE)  
  }else if(family %in% c(23,24,26)){
    cdf = function(u,v) v-.C("copCdf",as.double(1-u),as.double(v),as.integer(1),as.double(c(-par,par2)),as.integer(family-20),as.double(0),PACKAGE='CDRVine')[[6]]
    res = mapply(cdf, u1, u2, SIMPLIFY=TRUE) 
  }else if(family %in% c(33,34,36)){
    cdf = function(u,v) u-.C("copCdf",as.double(u),as.double(1-v),as.integer(1),as.double(c(-par,par2)),as.integer(family-30),as.double(0),PACKAGE='CDRVine')[[6]]
    res = mapply(cdf, u1, u2, SIMPLIFY=TRUE)
  }
 
  return(res)
}
