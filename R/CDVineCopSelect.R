CDVineCopSelect <- function(data,familyset=NA,type,selectioncrit="AIC",indeptest=FALSE,level=0.05){

	d = dim(data)[2]
	n = nrow(data)

	if(n<2) stop("Number of observations has to be at least 2.")
	if(d<3) stop("Dimension has to be at least 3.")
  if(!is.na(familyset[1])) for(i in 1:length(familyset)) if(!(familyset[i] %in% c(0,1,2,3,4,5,6,7,9,13,14,16,23,24,26,33,34,36))) stop("Copula family not implemented.")
  if(selectioncrit != "AIC" && selectioncrit != "BIC") stop("Selection criterion not implemented.")
  if(level < 0 & level > 1) stop("Significance level has to be between 0 and 1.")
  
  if(type == "CVine") type = 1
  else if(type == "DVine") type = 2
  if(type != 1 & type != 2) stop("Vine model not implemented.") 

  rhoMat = matrix(0,nrow=d-1,ncol=d-1)
  nuMat = matrix(0,nrow=d-1,ncol=d-1)  
  w = matrix(0,d-1,d-1)

if(type==1)# C-Vine
  {
	v = array(0,c(d-1,d-1,n))

	for(i in 1:(d-1))
	{
      par.out <- BiCopSelect(data[,1],data[,i+1],familyset,selectioncrit,indeptest,level)
			w[1,i] <- par.out$family
      rhoMat[1,i]  <- par.out$par
      nuMat[1,i]  <- par.out$par2
		
		v[1,i,] = .C("Hfunc1",as.integer(w[1,i]),as.integer(n),as.double(data[,i+1]), as.double(data[,1]), as.double(rhoMat[1,i]),as.double(nuMat[1,i]), as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
	}
	for(j in 2:(d-1))
	{
		for(i in 1:(d-j))
		{
        par.out <- BiCopSelect(v[j-1,1,],v[j-1,i+1,],familyset,selectioncrit,indeptest,level)
				w[j,i] <- par.out$family 
        rhoMat[j,i]  <- par.out$par
				nuMat[j,i]  <- par.out$par2

			if(j<(d-1))
			{
			v[j,i,] = .C("Hfunc1",as.integer(w[j,i]),as.integer(n),as.double(v[j-1,i+1,]), as.double(v[j-1,1,]), as.double(rhoMat[j,i]),as.double(nuMat[j,i]), as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
			}
		}
	}
  }
else # D-Vine
{
    v = array(0,c(d,2*d-4,n))

    for(i in 1:(d-1))
		{
      par.out <- BiCopSelect(data[,i],data[,i+1],familyset,selectioncrit,indeptest,level)
			w[1,i] <- par.out$family
			rhoMat[1,i]  <- par.out$par
			nuMat[1,i] <- par.out$par2

    }
    v[1,1,] = .C("Hfunc2",as.integer(w[1,1]),as.integer(n),as.double(data[,1]), as.double(data[,2]), as.double(rhoMat[1,1]),as.double(nuMat[1,1]), as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
    if(d>3)
    {
      for(k in 1:(d-3))
      {
        v[1,2*k,]  =  .C("Hfunc1",as.integer(w[1,k+1]),as.integer(n),as.double(data[,k+2]),as.double(data[,k+1]),as.double(rhoMat[1,k+1]),as.double(nuMat[1,k+1]),as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
        v[1,2*k+1,] = .C("Hfunc2",as.integer(w[1,k+1]),as.integer(n),as.double(data[,k+1]),as.double(data[,k+2]),as.double(rhoMat[1,k+1]),as.double(nuMat[1,k+1]),as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
      }
    }
    v[1,2*d-4,] = .C("Hfunc1",as.integer(w[1,d-1]),as.integer(n),as.double(data[,d]),as.double(data[,d-1]),as.double(rhoMat[1,d-1]),as.double(nuMat[1,d-1]),as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
    for(j in 2:(d-1))
    {
		for(i in 1:(d-j))
		{
        par.out <- BiCopSelect(v[j-1,2*i-1,],v[j-1,2*i,],familyset,selectioncrit,indeptest,level)
				w[j,i] <- par.out$family
				rhoMat[j,i]  <- par.out$par
				nuMat[j,i] <- par.out$par2

		}
      v[j,1,] = .C("Hfunc2",as.integer(w[j,1]),as.integer(n),as.double(v[j-1,1,]),as.double(v[j-1,2,]),as.double(rhoMat[j,1]),as.double(nuMat[j,1]),as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
      if(d > 4 & (d-j-2)>0)
      {
        for(i in 1:(d-j-2))
        {
          v[j,2*i,] = .C("Hfunc1",as.integer(w[j,i+1]),as.integer(n),as.double(v[j-1,2*i+2,]),as.double(v[j-1,2*i+1,]),as.double(rhoMat[j,i+1]),as.double(nuMat[j,i+1]),as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
          v[j,2*i+1,] = .C("Hfunc2",as.integer(w[j,i+1]),as.integer(n),as.double(v[j-1,2*i+1,]),as.double(v[j-1,2*i+2,]),as.double(rhoMat[j,i+1]),as.double(nuMat[j,i+1]),as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
        }
      }
      v[j,2*d-2*j-2,] = .C("Hfunc1",as.integer(w[j,d-j]),as.integer(n),as.double(v[j-1,2*d-2*j,]),as.double(v[j-1,2*d-2*j-1,]),as.double(rhoMat[j,d-j]),as.double(nuMat[j,d-j]),as.double(rep(0,n)),PACKAGE='CDVine')[[7]]
    }
}

theta0 = rep(0,d*(d-1)/2)
fam0 = rep(0,d*(d-1)/2)
k = 1
for(j in 1:(d-1)){
  for(i in 1:(d-j)){
    fam0[k] = w[j,i]
    theta0[k] = rhoMat[j,i]
    k = k+1
  }
}

tt = sum(fam0==2)+sum(fam0==7)+sum(fam0==8)+sum(fam0==9)
nu0=rep(0,tt)
k=1
for(j in 1:(d-1))
    {
      for(i in 1:(d-j))
      {
	if(w[j,i]==2 | w[j,i]==7 | w[j,i]==8 | w[j,i]==9)
	{
		nu0[k] = nuMat[j,i]
		k = k+1
	}
      }
    }

nu1=numeric()
kk=1
dd=d*(d-1)/2
for(k in 1:dd)
{
	if(fam0[k]==2 || fam0[k]==7 || fam0[k]==8 || fam0[k]==9)
	{
		nu1[k]=nu0[kk]
		kk=kk+1
	}
	else
	{
		nu1[k]=0
	}
}

  out=list(family=fam0, par=theta0, par2=nu1)
  return(out)
}

