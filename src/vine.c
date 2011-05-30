/*
** vine.c - C code of the package CDVine  
** 
** with contributions from Carlos Almeida, Aleksey Min, 
** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
** 
** A first version was based on code
** from Daniel Berg <daniel at danielberg.no>
** provided by personal communication. 
**
*/

#include "vine.h"


#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4
                                                                                                
///////////////////////////////////////////////////////////////////////////////
//  Function that allocates space and creates a double matrix.
//  Input: Dimension of the matrix to be created//  Output: Pointer to the created matrix.
///////////////////////////////////////////////////////////////////////////////
double **create_matrix(int rows, int columns)
{
  double **a;
  int i=0;
  a = (double**) Calloc(rows, double*);
  for(i=0;i<rows;i++) a[i] = (double*) Calloc(columns,double);
  return a;
}

///////////////////////////////////////////////////////////////////////////////
//  Function that frees the space that a double matrix has been allocated.
//  Input: Dimension of the matrix and a pointer to the matrix.
//  Output: Void.
///////////////////////////////////////////////////////////////////////////////
void free_matrix(double **a, int rows)
{
  int i=0;
  for(i=0;i<rows;i++) Free(a[i]);
  Free(a);
}

///////////////////////////////////////////////////////////////////////////////
//  Function that allocates space and creates an int matrix.
//  Input: Dimension of the matrix to be created.
//  Output: Pointer to the created matrix.
///////////////////////////////////////////////////////////////////////////////
int **create_intmatrix(int rows, int columns)
{
  int **a;
  int i=0;
  a = (int**) Calloc(rows,int*);
  for(i=0;i<rows;i++) a[i] = (int*) Calloc(columns,int);
  return a;
}

///////////////////////////////////////////////////////////////////////////////
//  Function that frees the space that an int matrix has been allocated.
//  Input: Dimension of the matrix and a pointer to the matrix.
//  Output: Void.
///////////////////////////////////////////////////////////////////////////////
void free_intmatrix(int **a, int rows)
{
  int i=0;
  for(i=0;i<rows;i++) Free(a[i]);
  Free(a);
}

///////////////////////////////////////////////////////////////////////////////
//  Function that allocates space and creates a 3-d double array.
//  Input: Dimensions of the array to be created.
//  Output: Pointer to the created array.
///////////////////////////////////////////////////////////////////////////////
double ***create_3darray(int d1, int d2, int d3)
{
  double ***a;
  int i=0,j=0;
  a = (double ***) Calloc(d1,double*);
  for(i=0;i<d1;i++)
  {  
    a[i] = (double**) Calloc(d2, double*);
    for(j=0;j<d2;j++)
    {
      a[i][j] = (double*) Calloc(d3,double);
    }
  }
  return a;
}



///////////////////////////////////////////////////////////////////////////////
//  Function that frees the space that a 3-d double array has been allocated.
//  Input: Dimensions of the array and a pointer to the array.
//  Output: Void.
///////////////////////////////////////////////////////////////////////////////
void free_3darray(double ***a, int d1, int d2)
{
  int i=0,j=0;
  for(i=0;i<d1;i++)
  {  
    for(j=0;j<d2;j++)
    {
      Free(a[i][j]);
    }
    Free(a[i]);
  }
  Free(a);
}






//////////////////////////////////////////////////////////////////////////////
// Print error text and return to R
//////////////////////////////////////////////////////////////////////////////
void printError(char *text, char filename[200])
{
  Rprintf(text);
  Rprintf(": %s ", filename);
  Rprintf(" !!!\n");
  exit(0);
}

//////////////////////////////////////////////////////////////////////////////
// Compute gamma division in a more stable way than gamma(x1)/gamma(x2)
// Input:
// x1 - Divisor
// x2 - Denominator
//-------------------
// a1 - modulus of a
// a2 - integer part of a
// b1 - modulus of b
// b2 - integer part of b
//-------------------
// We are after gamma(x1)/gamma(x2) and this computation will hopefully make it more numerically 
// stable. gamma(x1)/gamma(x2) will sometimes be INF/INF.
//////////////////////////////////////////////////////////////////////////////
double StableGammaDivision(double x1, double x2)
{
  int i;
  double a1, a2, b1, b2, sum=1.0;
  a1 = fmod(MAX(x1,x2),1.0);
  a2 = MAX(x1,x2)-a1;
  b1 = fmod(MIN(x1,x2),1.0);
  b2 = MIN(x1,x2)-b1;
  if(a1==0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=b2 ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
  }
  else if(a1>0.0 && b1==0.0)
  {
    for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=(int)b2 ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gammafn(a1);
  }
  else if(a1==0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum /= gammafn(b1);
  }
  else if(a1>0.0 && b1>0.0)
  {
    for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
    for(i=((int)b2+1) ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    sum *= gammafn(a1)/gammafn(b1);
  }
  if(x2 > x1) sum = 1.0/sum;
  return sum;
}


//////////////////////////////////////////////////////////////
// Generatorfunction of BB1, BB2 and BB7
// Input:
// u		variable
// n		number of iterations
// param	vector of parameter (theta, delta)
// copula	copula family (7=BB1, 8=BB2, 9=BB7)
// out		outout
//////////////////////////////////////////////////////////////
void gen(double* u, int* n, double* param, int* copula, double* out)
{
	int j;
	double *h;
	h = Calloc(*n,double);

	for(j=0;j<*n;j++)
	{
		if(u[j]==0) h[j] = 0;
		else if (u[j]==1) h[j] = u[j];
		else
		{
			if(*copula==3)	//Clayton
			{
				h[j] = 1/param[0]*(pow(u[j],(-param[0]))-1);
			}
			if(*copula==4)	//Gumbel
			{
				h[j] = pow((-log(u[j])),param[0]);
			}
			if(*copula==5)	//Frank
			{
				h[j] = -log((exp(-param[0]*u[j])-1)/(exp(-param[0])-1));
			}
			if(*copula==6)	//Joe
			{
				h[j] = -log(1-pow((1-u[j]),param[0]));
			}			
			if(*copula==7)	//BB1
			{
				h[j] = pow((pow(u[j],(-param[0]))-1),param[1]);
			}
			else if(*copula==8) //BB2
			{
				h[j] = exp(param[1]*(pow(u[j],-param[0])-1))-1;
			}
			else if(*copula==9)	//BB7
			{
				h[j] = pow(1-pow(1-u[j],param[0]),-param[1])-1;
			}
		}
		out[j]=h[j];
	}
	Free(h);
}


//////////////////////////////////////////////////////////////
// Inverse generator of BB1, BB2 and BB7
// Input:
// u		variable
// n		number of iterations
// param	vector of parameter (theta, delta)
// copula	copula family (7=BB1, 8=BB2, 9=BB7)
// out		outout
//////////////////////////////////////////////////////////////
void genInv(double* u, int* n, double* param, int* copula, double* out)
{
	int j;
	double *h;
	h = Calloc(*n,double);

	for(j=0;j<*n;j++)
	{
		if(u[j]==0) h[j] = 0;
		else if (u[j]==1) h[j] = u[j];
		else
		{
			if(*copula==3)	//Clayton
			{
				h[j] = pow((1+param[0]*u[j]),(-1/param[0]));
			}
			if(*copula==4)	//Gumbel
			{
				h[j] = exp(-pow(u[j],1/param[0]));
			}
			if(*copula==5)	//Frank
			{
				h[j] = -1/param[0]*log(1-exp(-u[j])*(1-exp(-param[0])));
			}
			if(*copula==6)	//Joe
			{
				h[j] = 1-pow((1-exp(-u[j])),1/param[0]);
			}
			if(*copula==7)	//BB1
			{
				h[j] = pow(1+pow(u[j],1/param[1]),(-1/param[0]));
			}
			else if(*copula==8) //BB2
			{
				h[j] = pow(1+1/param[1]*log(1+u[j]),-1/param[0]);
			}
			else if(*copula==9)	//BB7
			{
				h[j] = 1-pow(1-pow(1+u[j],-1/param[1]),(1/param[0]));
			}
		}
		out[j]=h[j];
	}
	Free(h);
}


//////////////////////////////////////////////////////////////
// First derivative of the generator of BB1, BB2 and BB7
// Input:
// u		variable
// n		number of iterations
// param	vector of parameter (theta, delta)
// copula	copula family (7=BB1, 8=BB2, 9=BB7)
// out		outout
//////////////////////////////////////////////////////////////
void genDrv(double* u, int* n, double* param, int* copula, double* out)
{
	int j;
	double *h;
	h = Calloc(*n,double);
	
	for(j=0;j<*n;j++)
	{
		if(u[j]==0) h[j] = 0;
		else if (u[j]==1) h[j] = u[j];
		else
		{
			if(*copula==7)	//BB1
			{
				h[j] = -(param[0]*param[1])*pow(pow(u[j],-param[0])-1,param[1]-1)*pow(u[j],-1-param[0]);
			}
			else if(*copula==8) //BB2
			{
				h[j] = -param[0]*param[1]*pow(u[j],-param[0]-1)*exp(param[1]*(pow(u[j],-param[0])-1));
			}
			else if(*copula==9)	//BB7
			{
				h[j] = -(param[0]*param[1])*pow(1-u[j],param[0]-1)*pow(1-pow(1-u[j],param[0]),-1-param[1]);
			}
		}
		out[j]=h[j];
	}
	Free(h);
}


//////////////////////////////////////////////////////////////
// Second derivative of the generator of BB1, BB2 and BB7
// Input:
// u		variable
// n		number of iterations
// param	vector of parameter (theta, delta)
// copula	copula family (7=BB1, 8=BB2, 9=BB7)
// out		outout
//////////////////////////////////////////////////////////////
void genDrv2(double* u, int* n, double* param, int* copula, double* out)
{
	int j;
	double *h;
	h = Calloc(*n,double);
	
	for(j=0;j<*n;j++)
	{
		if(u[j]==0) h[j] = 0;
		else if (u[j]==1) h[j] = u[j];
		else
		{
			if(*copula==7)	//BB1
			{
				h[j] = param[0]*param[1]*pow(u[j],-2-param[0])*pow(pow(u[j],-param[0])-1,param[1]-2)*((1+param[0]*param[1])*pow(u[j],-param[0])-param[0]-1);
			}
			else if(*copula==8) //BB2
			{
				h[j] = param[0]*param[1]*pow(u[j],-param[0]-2)*exp(param[1]*(pow(u[j],-param[0])-1))*(param[0]+1+param[1]*param[0]*pow(u[j],-param[0]));
			}
			else if(*copula==9)	//BB7
			{
				h[j] = param[0]*param[1]*pow(1-u[j],param[0]-2)*pow(1-pow(1-u[j],param[0]),-2-param[1])*((1+param[0]*param[1])*pow(1-u[j],param[0])+param[0]-1);
			}
		}
		out[j]=h[j];
	}
	Free(h);
}


//////////////////////////////////////////////////////////////
// Copula of BB1, BB2 and BB7
// Input:
// u		variable 1
// v		variable 2
// n		number of iterations
// param	vector of parameter (theta, delta)
// copula	copula family (7=BB1, 8=BB2, 9=BB7)
// out		outout
//////////////////////////////////////////////////////////////
void copCdf(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double *out1;
	double *out2;
	double *out3;
	out1 = Calloc(*n,double);
	out2 = Calloc(*n,double);
	out3 = Calloc(*n,double);
	gen(u, n, param, copula, out1);
	gen(v, n, param, copula, out2);
	for(j=0;j<*n;j++)
	{
	out3[j]=out1[j]+out2[j];
	}
	genInv(out3 , n, param, copula, out);
	Free(out1);
	Free(out2);
	Free(out3);
}


//////////////////////////////////////////////////////////////
// Copula density of BB1, BB2 and BB7
// Input:
// u		variable 1
// v		variable 2
// n		number of iterations
// param	vector of parameter (theta, delta)
// copula	copula family (7=BB1, 8=BB2, 9=BB7)
// out		outout
//////////////////////////////////////////////////////////////
void copPdf(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double *out1, *out2, *out3, *out4, *out5;
	out1 = Calloc(*n,double);
	out2 = Calloc(*n,double);
	out3 = Calloc(*n,double);
	out4 = Calloc(*n,double);
	out5 = Calloc(*n,double);
	copCdf(u,v,n,param,copula,out1);
	genDrv2(out1,n,param,copula,out2);
	genDrv(u,n,param,copula,out3);
	genDrv(v,n,param,copula,out4);
	genDrv(out1,n,param,copula,out5);
	for(j=0;j<*n;j++)
	{
	out[j]=-(out2[j]*out3[j]*out4[j])/pow(out5[j],3);
	}
	Free(out1);
	Free(out2);
	Free(out3);
	Free(out4);
	Free(out5);
}


//////////////////////////////////////////////////////////////
// h-function of BB1, BB2 and BB7
// Input:
// u		variable 1
// v		variable 2
// n		number of iterations
// param	vector of parameter (theta, delta)
// copula	copula family (7=BB1, 8=BB2, 9=BB7)
// out		output
//////////////////////////////////////////////////////////////
void hfuncbb(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double *out1, *out2, *out3;
	out1 = Calloc(*n,double);
	out2 = Calloc(*n,double);
	out3 = Calloc(*n,double);
	copCdf(u,v,n,param,copula,out1);
	genDrv(v,n,param,copula,out2);
	genDrv(out1,n,param,copula,out3);
	for(j=0;j<*n;j++)
	{
	out[j]=out2[j]/out3[j];
	}
	Free(out1);
	Free(out2);
	Free(out3);
}


/////////////////////////////////////////////////////////////////
// derivative of the h-function of BB2
// Input:
// u		variable 1
// v		variable 2
// n		number of iterations
// param	vector of parameter (theta, delta)
// copula	copula family (7=BB1, 8=BB2, 9=BB7)
// out		output
//////////////////////////////////////////////////////////////
void hdbb(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double *out1, *out2, *out3, *out4, *out5, *out6;
	out1 = Calloc(*n,double);
	out2 = Calloc(*n,double);
	out3 = Calloc(*n,double);
	out4 = Calloc(*n,double);
	out5 = Calloc(*n,double);
	out6 = Calloc(*n,double);
	copCdf(u,v,n,param,copula,out1);
	genDrv2(v,n,param,copula,out2);
	genDrv2(out1,n,param,copula,out3);
	genDrv(out1,n,param,copula,out4);
	genDrv(v,n,param,copula,out5);
	hfuncbb(u, v, n, param, copula, out6);
	for(j=0;j<*n;j++)
	{
	out[j]=(out2[j]*out4[j]-out5[j]*out3[j]*out6[j])/pow(out4[j],2);
	}
	Free(out1);
	Free(out2);
	Free(out3);
	Free(out4);
	Free(out5);
	Free(out6);
}


// Since the h function is not symmetric in case of double Gumbel and double Clayton we have two implement both separately,
// i.e. Hfunc1 and Hfunc2
void  Hfunc1(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out)
{
  double *negv, *negu;
  negv = (double *) malloc(*n* sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta;
  int nfamily;
  ntheta = -*theta;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  if(((*family==23) | (*family==24) | (*family==26)))
    {
	  nfamily=(*family)-20;
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc (&nfamily, n, u, negv, &ntheta, nu, out);
    }
  else if(((*family==33) | (*family==34) | (*family==36)))
	{
	  nfamily=(*family)-30;
      for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
      Hfunc(&nfamily, n, negu, v, &ntheta, nu, out);
		for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
	}
  else {
    Hfunc (family, n, u, v, theta, nu, out);
  }
  free(negv);
  free(negu);
}

void  Hfunc2(int* family,int* n,double* v,double* u,double* theta,double* nu,double* out)
{
  double *negv, *negu;
  negv = (double *) malloc(*n * sizeof(double));
  negu = (double *) malloc(*n * sizeof(double));
  double ntheta;
  int nfamily;
  ntheta = -*theta;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  if(((*family==23) | (*family==24) | (*family==26)))
    {
	  nfamily=(*family)-20;
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc(&nfamily, n, negv, u, &ntheta, nu, out);
		for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
  else if(((*family==33) | (*family==34) | (*family==36)))
	{
	  nfamily=(*family)-30;
      for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
      Hfunc(&nfamily, n, v, negu, &ntheta, nu, out);
	  //*out = 1 - *out;
	}
  else
    { 
      Hfunc(family, n, v, u, theta, nu, out);
    }
  free(negv);
  free(negu);
}



//////////////////////////////////////////////////////////////
// Function to compute h-function for vine simulation and estimation
// Input:
// family   copula family (0=independent,  1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1, 8=BB7)
// n        number of iterations
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// nu       degrees-of-freedom for the students copula
// out      output
//////////////////////////////////////////////////////////////
void Hfunc(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
  int j;
  double *h;
  h = Calloc(*n,double);
  double x;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  for(j=0;j<*n;j++)
  {
    if((v[j]==0) | ( u[j]==0)) h[j] = 0;
    else if (v[j]==1) h[j] = u[j];
    else
      {
    	if(*family==0) //independent
    	{
    	  h[j] = u[j]; 
    	}
    	else if(*family==1) //gaussian
    	{
    	  x = (qnorm(u[j],0.0,1.0,1,0) - *theta*qnorm(v[j],0.0,1.0,1,0))/sqrt(1.0-pow(*theta,2.0));
    	  if (isfinite(x))
		h[j] = pnorm(x,0.0,1.0,1,0);
    	  else if ((qnorm(u[j],0.0,1.0,1,0) - *theta*qnorm(v[j],0.0,1.0,1,0)) < 0)
		h[j] = 0;
    	  else 
		h[j] = 1;
    	}
    	else if(*family==2) //student
    	{
    	  double t1, t2, mu, sigma2;
    	  t1 = qt(u[j],*nu,1,0); t2 = qt(v[j],*nu,1,0); mu = *theta*t2; sigma2 = ((*nu+t2*t2)*(1.0-*theta*(*theta)))/(*nu+1.0);
    	  h[j] = pt((t1-mu)/sqrt(sigma2),*nu+1.0,1,0);
    	}
    	else if(*family==3) //clayton
    	  { 
			if(*theta == 0) h[j] = u[j] ;
			if(*theta < XEPS) h[j] = u[j] ;
			else
			{ 
		    x = pow(u[j],-*theta)+pow(v[j],-*theta)-1.0 ;
		    h[j] =   pow(v[j],-*theta-1.0)*pow(x,-1.0-1.0/(*theta)); // pow(v[j],-*theta-1.0)*pow(pow(u[j],-*theta)+pow(v[j],-*theta)-1.0,-1.0-1.0/(*theta));
		    if(*theta < 0)
		      {
				if(x < 0) h[j] = 0; 
		      }
			}
    	  }
    	else if(*family==4) //gumbel
    	{
    	  	if(*theta == 1) h[j] = u[j] ; 
			else
			{
		    h[j] = -(exp(-pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)))*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)-1.0)*pow(-log(v[j]),*theta))/(v[j]*log(v[j]));
			}
    	}
    	else if(*family==5) //frank
    	{
			if(*theta==0) h[j]=u[j];
			else
			{
    		h[j] = -(exp(*theta)*(exp(*theta*u[j])-1.0))/(exp(*theta*v[j]+*theta*u[j])-exp(*theta*v[j]+*theta)-exp(*theta*u[j]+*theta)+exp(*theta));
    		}
		}
    	else if(*family==6) //joe
    	{
			if(*theta==1) h[j]=u[j];
			else
			{
			h[j] = pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
    	  //h[j] = 1.0-pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta));
			}
    	}
		else if(*family==7)	//BB1
		{
			double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
			if(*nu==1) 
				{
				if(*theta==0) h[j]=u[j];
				else h[j]=pow(pow(u[j],-*theta)+pow(v[j],-*theta)-1,-1/(*theta)-1)*pow(v[j],-*theta-1);
				}
			else
			{
			hfuncbb(&u[j],&v[j],&T,param,family,&h[j]);
			}
			Free(param);
		}
		else if(*family==8) //BB2
		{
			double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
			
			hfuncbb(&u[j],&v[j],&T,param,family,&h[j]);

			Free(param);
		}
		else if(*family==9)	//BB7
		{
			double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
			if(*theta==1)
			{
				if(*nu==0) h[j]=u[j];
				else h[j]=pow(pow(u[j],-*nu)+pow(v[j],-*nu)-1,-1/(*nu)-1)*pow(v[j],-*nu-1);
			}
			else
			{
			hfuncbb(&u[j],&v[j],&T,param,family,&h[j]);
			}
			Free(param);
		}
		else if(*family==13) //rotated clayton (180°)
		  {
				u[j]=1-u[j]; 
				v[j]=1-v[j];
		    x = pow(u[j],-*theta)+pow(v[j],-*theta)-1.0 ;
		    h[j] =   pow(v[j],-*theta-1.0)*pow(x,-1.0-1.0/(*theta)); // pow(v[j],-*theta-1.0)*pow(pow(u[j],-*theta)+pow(v[j],-*theta)-1.0,-1.0-1.0/(*theta));
				h[j]= 1-h[j];
				u[j]=1-u[j];
				v[j]=1-v[j];
		  }
		else if(*family==14) //rotated gumbel (180°)
			{
					v[j]= 1-v[j];
					u[j]= 1-u[j];
					h[j]= -(exp(-pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)))*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)-1.0)*pow(-log(v[j]),*theta))/(v[j]*	log(v[j]));
					h[j]= 1-h[j];
					u[j]=1-u[j];
					v[j]=1-v[j];
			}
		else if(*family==16)
		  {
				v[j]= 1-v[j];
				u[j]= 1-u[j];
				h[j] = pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
    			h[j]= 1-h[j];
				u[j]=1-u[j];
				v[j]=1-v[j];
		  }
		
    }	
    out[j] = MAX(MIN(h[j],UMAX),UMIN);
  }	
  Free(h);
}



//////////////////////////////////////////////////////////////
// New function to compute log-likelihood for bivariate copula (for the rotated copulas
// Input:
// family    copula family (0=independent, 1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
// n         sample size
// u         first variable of data set
// v         second variable of data set
// theta     dependency parameter
// nu        degrees-of-freedom for students copula
// loglik    output
//////////////////////////////////////////////////////////////

void LL_mod(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
  double* negv;
  double* negu;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta;
  int nfamily;
  ntheta = -*theta;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  if(((*family==23) | (*family==24) | (*family==26)) )	// 90° rotated copulas
    {
	  nfamily = (*family)-20;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      LL(&nfamily, n, u,  negv, &ntheta, nu, loglik);
    }
  else if(((*family==33) | (*family==34) | (*family==36)) )	// 270° rotated copulas
    {
	  nfamily = (*family)-30;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      LL(&nfamily, n, negu,  v, &ntheta, nu, loglik);
    }
  else {
    LL(family, n, u,  v, theta, nu, loglik);
  }
  free(negv);
  free(negu);
}

void LL_mod2(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
  double* negv;
  double* negu;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta;
  int nfamily;
  ntheta = -*theta;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  if(((*family==23) | (*family==24) | (*family==26)))	// 90° rotated copulas
    {
	  nfamily = (*family)-20;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      LL(&nfamily, n, negu,  v, &ntheta, nu, loglik);
    }
  else if(((*family==33) | (*family==34) | (*family==36)))	// 270° rotated copulas
    {
	  nfamily = (*family)-30;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      LL(&nfamily, n, u,  negv, &ntheta, nu, loglik);
    }
  else {
    LL(family, n, u,  v, theta, nu, loglik);
  }
  free(negv);
  free(negu);
}


//////////////////////////////////////////////////////////////
// Function to compute log-likelihood for bivariate copula
// Input:
// family    copula family (0=independent, 1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
// n         sample size
// u         first variable of data set
// v         second variable of data set
// theta     dependency parameter
// nu        degrees-of-freedom for students copula
// loglik    output
//////////////////////////////////////////////////////////////
void LL(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
  int j;
  double *dat, rho, ll=0.0, t1=0.0, t2=0.0, f;
  //Allocate memory:
  dat = Calloc(2,double);

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  //Compute log-likelihood:
  if(*family==0) //independent
    ll = 0;
  else if(*family==1) //Gaussian
  {
    rho=*theta;
    for(j=0;j<*n;j++)
    {
      dat[0]=u[j]; dat[1]=v[j];
      t1 = qnorm(dat[0],0.0,1.0,1,0); t2 = qnorm(dat[1],0.0,1.0,1,0);
      f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
      if(log(f)>XINFMAX) ll += log(XINFMAX);
//       else if(f==0.0) ll += -100.0;
      else ll += log(f);
    }
  }
  else if(*family==2) //Student
  {
    rho=*theta;
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      t1 = qt(dat[0],*nu,1,0); t2 = qt(dat[1],*nu,1,0);
      f = StableGammaDivision((*nu+2.0)/2.0,*nu/2.0)/(*nu*pi*sqrt(1.0-pow(rho,2.0))*dt(t1,*nu,0)*dt(t2,*nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(*nu*(1.0-pow(rho,2.0))),-(*nu+2.0)/2.0);
      if(log(f)>XINFMAX) ll += log(XINFMAX);
//       else if(f==0.0) ll += -100.0;
      else ll += log(f);
    }
  }
  else if(*family==3) //Clayton
  {
    if(*theta == 0) ll = 0; // ca : 0 if the parameter defines independent
	  else if(*theta < XEPS) ll = 0;
				// copula !!!! 
    else
      {
	  for(j=0;j<*n;j++)
	  {
	    dat[0] = u[j]; dat[1] = v[j];
	    //f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
		 f=log(1+*theta)-(1+*theta)*log(dat[0]*dat[1])-(2+1/(*theta))*log(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0);
	    //f = MAX(f,0);
	    //if(log(f)>XINFMAX) ll += log(XINFMAX);
	    //       else if(f==0.0) ll += -100.0;
	    //else ll += log(f);
		ll +=f;
	  }
      }
  }
  else if(*family==4) //Gumbel
  {
    for(j=0;j<*n;j++)
    {
//       f=-(pow(-log(u[j]),*theta)*exp(-1.0/pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)))*(*theta*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta))+pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)-1.0))*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),-2.0/(*theta)-2.0)*pow(-log(v[j]),*theta))/(u[j]*log(u[j])*v[j]*log(v[j]));
      dat[0] = u[j]; dat[1] = v[j];
      t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
      //t2 = exp(-pow(t1,1.0/(*theta)));
      //f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
	  f= -pow(t1,1/(*theta))+(2/(*theta)-2)*log(t1)+(*theta-1)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log(1+(*theta-1)*pow(t1,-1.0/(*theta)));
      //if(log(f)>XINFMAX) ll += log(XINFMAX);
//       else if(f==0.0) ll += -100.0;
      //else ll += log(f);
	  ll += f;
    }
  }
  else if(*family==5) // Frank
  {
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      f = (*theta*(exp(*theta)-1.0)*exp(*theta*dat[1]+*theta*dat[0]+*theta))/pow(exp(*theta*dat[1]+*theta*dat[0])-exp(*theta*dat[1]+*theta)-exp(*theta*dat[0]+*theta)+exp(*theta),2.0);
/*      t1 = dat[1]*exp(-*theta*dat[1])*(exp(-*theta*dat[0])-1.0);
      t2 = -*theta*dat[0]*dat[1]*exp(-*theta*dat[0])*exp(-*theta*dat[1]);
      n1 = (exp(-*theta)-1.0)+(exp(-*theta*dat[0])-1.0)*(exp(-*theta*dat[1])-1.0);
      n2 = -*theta*dat[0]*exp(-*theta*dat[0])*(exp(-*theta*dat[1])-1.0);
      f = (t2*n1 - n2*t1)/pow(n1,2.0);*/
      if(log(f)>XINFMAX) ll += log(XINFMAX);
//       else if(f==0.0) ll += -100.0;
      else ll += log(f);
    }
  }
  else if(*family==6)	//Joe
  {
	  for(j=0;j<*n;j++)
		{
		//dat[0]=u[j]; dat[1]=v[j];
		f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
		if(log(f)>XINFMAX) ll += log(XINFMAX);
//      else if(f==0.0) ll += -100.0;
		else ll += log(f);
		}
  }
  else if(*family==7)	//BB1
  {
		//double *param, *fuc;
		double part1, part2, part3, part4;
		//param=Calloc(2,double);
		//param[0]=*theta;
		//param[1]=*nu;
		//fuc = Calloc(*n,double);
		//copPdf(u, v, n, param, family, fuc);
		for(j=0;j<*n;j++)
		{
			part1=pow(1+pow(pow(pow(u[j],-*theta)-1,*nu)+pow(pow(v[j],-*theta)-1,*nu),1/(*nu)),-1/(*theta)-2);
			part2=pow(pow(pow(u[j],-*theta)-1,*nu)+pow(pow(v[j],-*theta)-1,*nu),2/(*nu)-2);
			part3=(*theta)*(*nu)+1+(*theta)*(*nu-1)*pow(pow(pow(u[j],-*theta)-1,*nu)+pow(pow(v[j],-*theta)-1,*nu),-1/(*nu));
			part4=pow(pow(u[j],-*theta)-1,*nu-1)*pow(u[j],-*theta-1)*pow(pow(v[j],-*theta)-1,*nu-1)*pow(v[j],-*theta-1);
			if(!isfinite(part1) || isnan(part1))
			{
				part1=1;
			}
			if(!isfinite(part2) || isnan(part2))
			{
				part2=1;
			}
			if(!isfinite(part3) || isnan(part3))
			{
				part3=1;
			}
			if(!isfinite(part4) || isnan(part4))
			{
				part4=1;
			}
			ll+=log(part1)+log(part2)+log(part3)+log(part4);
			//if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			//else if(fuc[j]==0.0 || fuc[j]<0.0) ll += -100.0;
			//else ll += log(fuc[j]);
		}
		//Free(fuc); 
		//Free(param);
	}
	else if(*family==8)	//BB2
	{
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		copPdf(u, v, n, param, family, fuc);
		for(j=0;j<*n;j++)
		{
			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			//else if(fuc[j]==0.0 || fuc[j]<0.0) ll += -100.0;
			else ll += log(fuc[j]);
		}
		Free(fuc); Free(param);
	}
	else if(*family==9)	//BB7
	{
		//double h, dvh, duh, duvh, duS, dvS;
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		copPdf(u, v, n, param, family, fuc);
		for(j=0;j<*n;j++)
		{
			/*
			h=1-pow(pow(1-pow(1-u[j],*theta),-*nu)-pow(1-pow(1-v[j],*theta),-*nu)-1,1/(*nu));
			duh=-*theta*pow(pow(1-pow(1-u[j],*theta),-*nu)-pow(1-pow(1-v[j],*theta),-*nu)-1,1/(*nu-1))*pow(1-pow(1-u[j],*theta),-*nu-1)*pow(1-u[j],*theta-1);
			dvh=-*theta*pow(pow(1-pow(1-u[j],*theta),-*nu)-pow(1-pow(1-v[j],*theta),-*nu)-1,1/(*nu-1))*pow(1-pow(1-v[j],*theta),-*nu-1)*pow(1-v[j],*theta-1);
			duS=-*theta*(*nu)*pow(1-pow(1-u[j],*theta),-*nu-1)*pow(1-u[j],*theta-1);
			dvS=-*theta*(*nu)*pow(1-pow(1-v[j],*theta),-*nu-1)*pow(1-v[j],*theta-1);
			duvh=1/(*nu)*(-1/(*nu)-1)*pow(pow(1-pow(1-u[j],*theta),-*nu)-pow(1-pow(1-v[j],*theta),-*nu)-1,1/(*nu)-2)*duS*dvS;

			f=(-1/(*theta))*(1/(*nu)-1)*pow(h,1/(*theta)-2)*dvh*duh-1/(*theta)*pow(h,1/(*theta)-1)*duvh;
			ll+=log(f);
			*/
			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			//else if(fuc[j]==0.0 || fuc[j]<0.0) ll += -100.0;
			else ll += log(fuc[j]);
		}
		Free(fuc); Free(param);
  }
  else if(*family==13) //rotated Clayton (180°)
  {
		if(*theta == 0) ll = 0; // ca : 0 if the parameter defines independent copula !!!! 
	  else if(*theta < XEPS) ll = 0;
		else
		{
		for(j=0;j<*n;j++)
		  {
			dat[0] = 1-u[j]; dat[1] = 1-v[j];
			f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
			f = MAX(f,0);
			if(log(f)>XINFMAX) ll += log(XINFMAX);
			else ll += log(f);
		  }
      }
  }
  else if(*family==14) //rotated Gumbel (180°)
  {
			for(j=0;j<*n;j++)
			{
				dat[0] = 1-u[j]; dat[1] = 1-v[j];
				t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
				t2 = exp(-pow(t1,1.0/(*theta)));
				f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
				if(log(f)>XINFMAX) ll += log(XINFMAX);
				else ll += log(f);
			}
  }
  else if(*family==16) //rotated Joe (180°)
  {
	for(j=0;j<*n;j++)
		{
		u[j]=1-u[j]; v[j]=1-v[j];
		f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
		if(log(f)>XINFMAX) ll += log(XINFMAX);
		else ll += log(f);
		u[j]=1-u[j]; v[j]=1-v[j];
		}
  }
  else 
  {
	  Rprintf("%d\n",*family);
	  printError("Error in LL\t:","Unknown copula family"); 
  }
  //Free memory:
  Free(dat);
  //Write to output vector:
  *loglik = ll;
}



//////////////////////////////////////////////////////////////
// Function to compute likelihood for bivariate copula
// Input:
// family    copula family (0=independent, 1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
// n         sample size
// u         first variable of data set
// v         second variable of data set
// theta     dependency parameter
// nu        degrees-of-freedom for students copula
// coplik    output
//////////////////////////////////////////////////////////////
void copLik_mod(int* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik)
{
	double* negv;
  double* negu;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta;
  int nfamily;
  ntheta = -*theta;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  if(((*family==23) | (*family==24) | (*family==26)) )	// 90° rotated copulas
    {
	  nfamily = (*family)-20;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      copLik(&nfamily, n, negu,  v, &ntheta, nu, coplik);
    }
  else if(((*family==33) | (*family==34) | (*family==36)) )	// 270° rotated copulas
    {
	  nfamily = (*family)-30;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      copLik(&nfamily, n, u,  negv, &ntheta, nu, coplik);
    }
  else {
    copLik(family, n, u,  v, theta, nu, coplik);
  }
  free(negv);
  free(negu);
}


void copLik(int* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik)
{
  int j;
  double *dat, rho, lik=1.0, t1=0.0, t2=0.0, f;
  //Allocate memory:
  dat = Calloc(2,double);

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  //Compute likelihood:
  if(*family==0) //independent
    lik = 1.0;
  else if(*family==1) //Gaussian
  {
    rho=*theta;
    for(j=0;j<*n;j++)
    {
      dat[0]=u[j]; dat[1]=v[j];
      t1 = qnorm(dat[0],0.0,1.0,1,0); t2 = qnorm(dat[1],0.0,1.0,1,0);
      f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
      lik *= f;
    }
  }
  else if(*family==2) //Student
  {
    rho=*theta;
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      t1 = qt(dat[0],*nu,1,0); t2 = qt(dat[1],*nu,1,0);
      f = StableGammaDivision((*nu+2.0)/2.0,*nu/2.0)/(*nu*pi*sqrt(1.0-pow(rho,2.0))*dt(t1,*nu,0)*dt(t2,*nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(*nu*(1.0-pow(rho,2.0))),-(*nu+2.0)/2.0);
      lik *= f;
    }
  }
  else if(*family==3) //Clayton
  {
    if(*theta == 0) lik = 1.0; // ca : 0 if the parameter defines independent
	if(*theta < XEPS) lik = 1.0;
				// copula !!!!
    else
      {
	for(j=0;j<*n;j++)
	  {
	    dat[0] = u[j]; dat[1] = v[j];
	    f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
	    f = MAX(f,0);
	    lik *= f;
	  }
      }
  }
  else if(*family==4) //Gumbel
  {
    for(j=0;j<*n;j++)
    {
//       f=-(pow(-log(u[j]),*theta)*exp(-1.0/pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)))*(*theta*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta))+pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)-1.0))*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),-2.0/(*theta)-2.0)*pow(-log(v[j]),*theta))/(u[j]*log(u[j])*v[j]*log(v[j]));
      dat[0] = u[j]; dat[1] = v[j];
      t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
      t2 = exp(-pow(t1,1.0/(*theta)));
      f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
      lik *= f;
    }
  }
  else if(*family==5) // Frank
  {
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      f = (*theta*(exp(*theta)-1.0)*exp(*theta*dat[1]+*theta*dat[0]+*theta))/pow(exp(*theta*dat[1]+*theta*dat[0])-exp(*theta*dat[1]+*theta)-exp(*theta*dat[0]+*theta)+exp(*theta),2.0);
/*      t1 = dat[1]*exp(-*theta*dat[1])*(exp(-*theta*dat[0])-1.0);
      t2 = -*theta*dat[0]*dat[1]*exp(-*theta*dat[0])*exp(-*theta*dat[1]);
      n1 = (exp(-*theta)-1.0)+(exp(-*theta*dat[0])-1.0)*(exp(-*theta*dat[1])-1.0);
      n2 = -*theta*dat[0]*exp(-*theta*dat[0])*(exp(-*theta*dat[1])-1.0);
      f = (t2*n1 - n2*t1)/pow(n1,2.0);*/
      lik *= f;
    }
  }
  else if(*family==6)	//Joe
  {
	  for(j=0;j<*n;j++)
		{
		//dat[0]=u[j]; dat[1]=v[j];
		f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
		lik *= f;
		}
  }
  else if(*family==7)	//BB1
  {
		//double *param, *fuc;
		double part1, part2, part3, part4;
		//param=Calloc(2,double);
		//param[0]=*theta;
		//param[1]=*nu;
		//fuc = Calloc(*n,double);
		//copPdf(u, v, n, param, family, fuc);
		for(j=0;j<*n;j++)
		{
			part1=pow(1+pow(pow(pow(u[j],-*theta)-1,*nu)+pow(pow(v[j],-*theta)-1,*nu),1/(*nu)),-1/(*theta)-2);
			part2=pow(pow(pow(u[j],-*theta)-1,*nu)+pow(pow(v[j],-*theta)-1,*nu),2/(*nu)-2);
			part3=(*theta)*(*nu)+1+(*theta)*(*nu-1)*pow(pow(pow(u[j],-*theta)-1,*nu)+pow(pow(v[j],-*theta)-1,*nu),-1/(*nu));
			part4=pow(pow(u[j],-*theta)-1,*nu-1)*pow(u[j],-*theta-1)*pow(pow(v[j],-*theta)-1,*nu-1)*pow(v[j],-*theta-1);
			if(!isfinite(part1) || isnan(part1))
			{
				part1=1;
			}
			if(!isfinite(part2) || isnan(part2))
			{
				part2=1;
			}
			if(!isfinite(part3) || isnan(part3))
			{
				part3=1;
			}
			if(!isfinite(part4) || isnan(part4))
			{
				part4=1;
			}
			lik *= part1*part2*part3*part4;
			//if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			//else if(fuc[j]==0.0 || fuc[j]<0.0) ll += -100.0;
			//else ll += log(fuc[j]);
		}
		//Free(fuc);
		//Free(param);
  }
  else if(*family==9)	//BB7
  {
		//double h, dvh, duh, duvh, duS, dvS;
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		copPdf(u, v, n, param, family, fuc);
		for(j=0;j<*n;j++)
		{
			/*
			h=1-pow(pow(1-pow(1-u[j],*theta),-*nu)-pow(1-pow(1-v[j],*theta),-*nu)-1,1/(*nu));
			duh=-*theta*pow(pow(1-pow(1-u[j],*theta),-*nu)-pow(1-pow(1-v[j],*theta),-*nu)-1,1/(*nu-1))*pow(1-pow(1-u[j],*theta),-*nu-1)*pow(1-u[j],*theta-1);
			dvh=-*theta*pow(pow(1-pow(1-u[j],*theta),-*nu)-pow(1-pow(1-v[j],*theta),-*nu)-1,1/(*nu-1))*pow(1-pow(1-v[j],*theta),-*nu-1)*pow(1-v[j],*theta-1);
			duS=-*theta*(*nu)*pow(1-pow(1-u[j],*theta),-*nu-1)*pow(1-u[j],*theta-1);
			dvS=-*theta*(*nu)*pow(1-pow(1-v[j],*theta),-*nu-1)*pow(1-v[j],*theta-1);
			duvh=1/(*nu)*(-1/(*nu)-1)*pow(pow(1-pow(1-u[j],*theta),-*nu)-pow(1-pow(1-v[j],*theta),-*nu)-1,1/(*nu)-2)*duS*dvS;

			f=(-1/(*theta))*(1/(*nu)-1)*pow(h,1/(*theta)-2)*dvh*duh-1/(*theta)*pow(h,1/(*theta)-1)*duvh;
			lik *= f;
			*/
			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			lik *= fuc[j];
		}
		Free(fuc); Free(param);
  }
  else if(*family==13) //rotated Clayton (180°)
  {
		if(*theta == 0) lik = 1.0; // ca : 0 if the parameter defines independent copula !!!! 
		else
		{
		for(j=0;j<*n;j++)
		  {
			dat[0] = 1-u[j]; dat[1] = 1-v[j];
			f = (1.0+*theta)*pow(dat[0]*dat[1],-1.0-*theta)*pow(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0,-2.0-1.0/(*theta));
			lik *= f;
		  }
      }
  }
  else if(*family==14) //rotated Gumbel (180°)
  {
			for(j=0;j<*n;j++)
			{
				//       f=-(pow(-log(u[j]),*theta)*exp(-1.0/pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)))*(*theta*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta))+pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),1.0/(*theta)-1.0))*pow(pow(-log(v[j]),*theta)+pow(-log(u[j]),*theta),-2.0/(*theta)-2.0)*pow(-log(v[j]),*theta))/(u[j]*log(u[j])*v[j]*log(v[j]));
				dat[0] = 1-u[j]; dat[1] = 1-v[j];
				t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
				t2 = exp(-pow(t1,1.0/(*theta)));
				f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
				//if(log(f)>XINFMAX) ll += log(XINFMAX);
				//       else if(f==0.0) ll += -100.0;
				//else
				lik *= f;
			}
  }
  else if(*family==16) //rotated Joe (180°)
  {
	for(j=0;j<*n;j++)
		{
		dat[0]=1-u[j]; dat[1]=1-v[j];
		f = pow(pow(1-dat[0],*theta)+pow(1-dat[1],*theta)-pow(1-dat[0],*theta)*pow(1-dat[1],*theta),1/(*theta)-2)*pow(1-dat[0],*theta-1)*pow(1-dat[1],*theta-1)*(*theta-1+pow(1-dat[0],*theta)+pow(1-dat[1],*theta)-pow(1-dat[0],*theta)*pow(1-dat[1],*theta));
		lik *= f;
		}
  }
  else printError("Error in copLik\t:","Unknown copula family");
  //Free memory:
  Free(dat);
  //Write to output vector:
  *coplik = lik;
}



///////////////////////////////////////////////////////////////
// Function that returns the derivative of H
// Input:
// family   copula family (4=gumbel, 5=frank, 6=joe, 7=BB1, 8=BB7)
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// nu       degrees-of-freedom for the students copula
// out      output
//////////////////////////////////////////////////////////////
void Hderiv(int* family, double* u, double* v, double* theta, double* nu, double* out)
{
  double hd=0.0;
  if(*family==4) //gumbel
  {
    hd = (pow(-log(*u),*theta)*exp(-pow(pow(-log(*v),*theta)+pow(-log(*u),*theta),1.0/(*theta)))*(pow(pow(-log(*v),*theta)+pow(-log(*u),*theta),1.0/(*theta))+*theta-1.0)*pow(pow(-log(*v),*theta)+pow(-log(*u),*theta),1.0/(*theta)-2.0)*pow(-log(*v),*theta))/(*u*log(*u)*(*v)*log(*v));
//     hd = exp(-pow(pow(-log(*u),*theta)+pow(-log(*v),*theta),1.0/(*theta)))/(*u)/(*v)*pow(log(*u)*log(*v),*theta-1.0)*pow(pow(-log(*u),*theta)+pow(-log(*v),*theta),-2.0+2.0/(*theta))*(1.0+(*theta-1.0)*pow(pow(-log(*u),*theta)+pow(-log(*v),*theta),-1.0/(*theta)));
  }
  else if(*family==5) //frank
  {
    hd = (*theta*(exp(*theta)-1.0)*exp(*theta*(*u+*v+1.0)))/pow(exp(*theta*(*u+*v))-exp(*theta*(*v+1.0))-exp(*theta*(*u+1.0))+exp(*theta),2.0);
//     hd = -( 1.0/(*theta) * ( (exp(-*theta * (*u)) - 1.0) * (exp(-*theta * (*v)) * (*theta) * (*theta)) / (exp(-*theta) - 1.0) / (1.0 + (exp(-*theta * (*u)) - 1.0) * (exp(-*theta * (*v)) - 1.0) / (exp(-*theta) - 1.0)) - (exp(-*theta * (*u)) - 1.0) * (exp(-*theta * (*v)) * (*theta))/(exp(-*theta) - 1.0) * ((exp(-*theta * (*u)) - 1.0) * (exp(-*theta * (*v)) * (*theta))/(exp(-*theta) - 1.0)) / pow(1.0 + (exp(-*theta * (*u)) - 1.0) * (exp(-*theta * (*v)) - 1.0) / (exp(-*theta) - 1.0),2.0) ) );
  }
  else if(*family==6) //joe
  {
	  hd = (1-*theta)*pow(pow(1-*u,*theta)+pow(1-*v,*theta)-pow(1-*u,*theta)*pow(1-*v,*theta),1/(*theta)-2)*pow(1-*v,*theta-1)*(1-pow(1-*u,*theta)) - (*theta-1)*pow(pow(1-*u,*theta)+pow(1-*v,*theta)-pow(1-*u,*theta)*pow(1-*v,*theta),1/(*theta)-1)*pow(1-*v,*theta-2)*(1-pow(1-*u,*theta));
    //hd = -( pow(1.0 - (1.0 - pow(1.0 - *u,*theta)) * (1.0 - pow(1.0 - *v,*theta)),1.0/(*theta)-1.0) * ((1.0 - pow(1.0 - *u,*theta)) * pow(1.0 - *v,*theta-2.0) * (*theta - 1.0)) + pow(1.0 - (1.0 - pow(1.0 - *u,*theta)) * (1.0 - pow(1.0 - *v,*theta)),1.0/(*theta)-2.0) * ((1.0 - *theta) * (1.0 - pow(1.0 - *u,*theta)) * pow(1.0 - *v,*theta-1.0)) * ((1.0 - pow(1.0 - *u,*theta)) * pow(1.0 - *v,*theta-1.0)) );  
  }
  else if(*family==7) //BB1
	{
		double h;
		double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
			if(*nu==1) 
				{
				if(*theta==0) h=*u;
				else h=pow(pow(*u,-*theta)+pow(*v,-*theta)-1,-1/(*theta)-1)*pow(*v,-*theta-1);
				}
			else
			{
			hfuncbb(u,v,&T,param,family,&h);
			}
			Free(param);

			h = MAX(MIN(h,UMAX),UMIN);
		hd = h*( (*theta-1)*pow(1+pow(pow(pow(*u,-*theta)-1,*nu)+pow(pow(*v,-*theta)-1,*nu),1/(*nu)),-1) * pow(pow(pow(*u,-*theta)-1,*nu)+pow(pow(*v,-*theta)-1,*nu),1/(*nu)-1)*pow(pow(*v,-*theta)-1,*nu-1)*pow(*v,-*theta-1) + (*nu-1)*(*theta)*pow(pow(pow(*u,-*theta)-1,*nu)+pow(pow(*v,-*theta)-1,*nu),-1)*pow(pow(*v,-*theta)-1,*nu-1)*pow(*v,-*theta-1) + (1-*nu)*(*theta)*pow(pow(*v,-*theta)-1,-1)*pow(*v,-*theta-1) + (-*theta-1)*pow(*v,-1));
	}
	else if(*family==8) //BB2
	{
		double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
			hdbb(u,v,&T,param,family,&hd);
			Free(param);
	}
	else if(*family==9) //BB7
	{
		double h;
		double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
			if(*theta==1)
			{
				if(*nu==0) h=*u;
				else h=pow(pow(*u,-*nu)+pow(*v,-*nu)-1,-1/(*nu)-1)*pow(*v,-*nu-1);
			}
			else
			{
			hfuncbb(u,v,&T,param,family,&h);
			}
			Free(param);
			h = MAX(MIN(h,UMAX),UMIN);
			hd = h*( (1-*theta)*pow(1-pow(pow(1-pow(1-*u,*theta),-*nu)+pow(1-pow(1-*v,*theta),-*nu)-1,-1/(*nu)),-1) * pow(pow(1-pow(1-*u,*theta),-*nu)+pow(1-pow(1-*v,*theta),-*nu)-1,-1/(*nu)-1)*pow(1-pow(1-*v,*theta),-*nu-1)*pow(1-*v,*theta-1) + (*nu+1)*(*theta)*pow(pow(1-pow(1-*u,*theta),-*nu)+pow(1-pow(1-*v,*theta),-*nu)-1,-1)*pow(1-pow(1-*v,*theta),-*nu-1)*pow(1-*v,*theta-1) + (-*nu-1)*(*theta)*pow(1-pow(1-*v,*theta),-1)*pow(1-*v,*theta-1) + (1-*theta)*pow(1-*v,-1) );	
	}
	else if(*family==14) //rotated gumbel (180°)
	{
		*u = 1-*u;
		*v = 1-*v ;
		hd = (pow(-log(*u),*theta)*exp(-pow(pow(-log(*v),*theta)+pow(-log(*u),*theta),1.0/(*theta)))*(pow(pow(-log(*v),*theta)+pow(-log(*u),*theta),1.0/(*theta))+*theta-1.0)*pow(pow(-log(*v),*theta)+pow(-log(*u),*theta),1.0/(*theta)-2.0)*pow(-log(*v),*theta))/(*u*log(*u)*(*v)*log(*v));
		//     hd = exp(-pow(pow(-log(*u),*theta)+pow(-log(*v),*theta),1.0/(*theta)))/(*u)/(*v)*pow(log(*u)*log(*v),*theta-1.0)*pow(pow(-log(*u),*theta)+pow(-log(*v),*theta),-2.0+2.0/(*theta))*(1.0+(*theta-1.0)*pow(pow(-log(*u),*theta)+pow(-log(*v),*theta),-1.0/(*theta)));
		*u = 1-*u;
		*v = 1-*v;
	}
	else if(*family==16)
	{
		*u = 1-*u;
		*v = 1-*v;
		hd = (1-*theta)*pow(pow(1-*u,*theta)+pow(1-*v,*theta)-pow(1-*u,*theta)*pow(1-*v,*theta),1/(*theta)-2)*pow(1-*v,*theta-1)*(1-pow(1-*u,*theta)) - (*theta-1)*pow(pow(1-*u,*theta)+pow(1-*v,*theta)-pow(1-*u,*theta)*pow(1-*v,*theta),1/(*theta)-1)*pow(1-*v,*theta-2)*(1-pow(1-*u,*theta));
		*u = 1-*u;
		*v = 1-*v;
	}
	
  *out = hd;
}


///////////////////////////////////////////////////////////////
// Function to compute inversion of H numerically through Newton-Raphson
// Input:
// family   copula family (4=gumbel, 5=frank, 6=joe)
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// out      output
//////////////////////////////////////////////////////////////
void HNumInv(int* family, double* u, double* v, double* theta, double* nu, double* out)
{
  int i, it=10000, br=0, in=1;
  double ans=0.0, tol=0.00001, x0=0.000001, x1=0.999999, xl, xh, rts, dx, dxold, fl, fh, f, df, temp;
  
  Hfunc(family,&in,&x0,v,theta,nu,&fl); fl -= *u; 
  Hfunc(family,&in,&x1,v,theta,nu,&fh); fh -= *u;
  if(fabs(fl)<=tol) { ans=x0; br=1; }
  if(fabs(fh)<=tol) { ans=x1; br=1; }
  
  if(!br)
  {
    //if((fl>0.0 && fh>0.0) || (fl<0.0 && fh<0.0)) printError("Error in HNumInv\t:","Root must be bracketed in rtsafe");
	if(fl>0.0 && fh>0.0){ans=x0; br=1;}
    if(fl<0.0 && fh<0.0){ans=x1; br=1;}


    if(fl<0.0)
    {
      xl=x0;
      xh=x1;
    }
    else
    {
      xh=x0;
      xl=x0;
    }
    rts=0.5*(x0+x1);
    dxold=fabs(x1-x0);
    dx=dxold;
    Hfunc(family,&in,&rts,v,theta,nu,&f); f -= *u; 
    Hderiv(family,&rts,v,theta,nu,&df);
    for(i=0;i<it;i++)
    {
      if((((rts-xh)*df-f)*((rts-xl)*df-f)>0.0) || (fabs(2.0*f)>fabs(dxold*df)))
      {
        dxold=dx;
        dx=0.5*(xh-xl);
        rts=xl+dx;
        if(xl==rts) { ans=rts; break; }
      }
      else
      {
        dxold=dx;
        dx=f/df;
        temp=rts;
        rts-=dx;
        if(temp==rts) { ans=rts; break; }
      }
      if(fabs(dx)<tol) { ans=rts; break; }
      Hfunc(family,&in,&rts,v,theta,nu,&f); f -= *u; 
      Hderiv(family,&rts,v,theta,nu,&df);
      if(f<0.0) xl=rts;
      else xh=rts;
    }
  }
  *out = ans;
}

///////////////////////////////////////////////////////////////
// Function to compute inversion of H numerically through Bisection
// Input:
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// out      output
//////////////////////////////////////////////////////////////
void HNumInvBisect(int* family, double* u, double* v, double* theta, double* nu, double* out)
{

  int br=0, in=1;
  double ans=0.0, tol=0.00001, x0=0.000001, x1=0.999999, fl, fh, val;
  
  Hfunc(family,&in,&x0,v,theta,nu,&fl); fl -= *u; 
  Hfunc(family,&in,&x1,v,theta,nu,&fh); fh -= *u;
  if(fabs(fl)<=tol) { ans=x0; br=1; }
  if(fabs(fh)<=tol) { ans=x1; br=1; }

  while(!br){
  
  		ans = (x0+x1)/2.0;
  		Hfunc(family,&in,&ans,v,theta,nu,&val);
  		val -= *u;
  		if(fabs(val)<=tol) br=1;
  		if(fabs(x0-x1)<=1e-10) br=1; //stop if values become too close (avoid infinite loop)
  		/*if(fl > 0.0)   // Hfunc is mon. increasing
  		{
  			if(val < 0.0) {x1 = ans; fh = val;}
  			else {x0 = ans; fl = val;}
 			}
  		if(fl <= 0.0)
  		{*/
  			if(val > 0.0) {x1 = ans; fh = val;}
  			else {x0 = ans; fl = val;}
 			//}
  
  }
	*out = ans;
}  

/////////////////////////////////////////////
// Function to invert h-function for vine simulation and estimation
/////////////////////////////////////////////
void Hinv1(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
  double *negv, *negu;
  negv = (double*) Calloc(*n,double);
  negu = (double*) Calloc(*n,double);
  double ntheta;
  int nfamily;
  ntheta = -*theta;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

   if(((*family ==23) | (*family ==24) | (*family==26)))
      {
	   nfamily=(*family)-20;
		for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
		Hinv(&nfamily,  n,  u,  negv,  &ntheta,  nu,  out);
      }
	else if(((*family==33) | (*family==34) | (*family==36)))
      {
	   nfamily=(*family)-30;
		for (int i = 0; i < *n; i++) {negu[i]=1 - u[i];};
		Hinv(&nfamily,  n,  negu,  v,  &ntheta,  nu,  out);
		//Hinv(&nfamily,  n,  v,  negu,  &ntheta,  nu,  out);
		for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
      }
   else {
     Hinv( family,  n,  u,  v,  theta,  nu,  out);
   }
   Free(negv);
   Free(negu);
}

void Hinv2(int* family, int* n, double* v, double* u, double* theta, double* nu, double* out)
{
  double *negv, *negu;
  negv = (double *) malloc(*n*sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta;
  int nfamily;
  ntheta = -*theta;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

   if(((*family ==23) | (*family ==24) | (*family==26)))
      {
	   nfamily = (*family)-20;
		for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
		Hinv(&nfamily,  n,  negv, u,  &ntheta,  nu,  out);
		for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
      }
	else if(((*family==33) | (*family==34) | (*family==36)))
      {
	   nfamily=(*family)-30;
		for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
		Hinv(&nfamily,  n,  v,  negu,  &ntheta,  nu,  out);
		//*out = 1-*out;
      }
   else {
     Hinv( family,  n,  v,  u,  theta,  nu,  out);
   }
   free(negv);
   free(negu);
}



//////////////////////////////////////////////////////////////
// Function to invert h-function for vine simulation and estimation
// Input:
// family   copula family (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// n        number of iterations
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// nu       degrees-of-freedom for the students copula
//////////////////////////////////////////////////////////////
void Hinv(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out)
{
  int j,in=1;
  double hmax, umax = 1-1e-6;
  double hmin, umin = 1e-6;
  double *hinv;
  hinv = Calloc(*n,double);

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  for(j=0;j<*n;j++)
  {
    Hfunc(family,&in,&umax,v,theta,nu,&hmax);
    Hfunc(family,&in,&umin,v,theta,nu,&hmin);
    if(*u >= hmax) *out = umax;
    else if(*u <= hmin) *out = umin;
	else if(*family==0)
	  {
		hinv[j]=u[j];
	  }
    else if(*family==1) //gaussian
      {
      hinv[j] = pnorm(qnorm(u[j],0.0,1.0,1,0)*sqrt(1.0-pow(*theta,2.0))+*theta*qnorm(v[j],0.0,1.0,1,0),0.0,1.0,1,0);
    }
    else if(*family==2) //student
    {
      double temp1, temp2, mu, var;
      temp1 = qt(u[j],*nu+1.0,1,0); temp2 = qt(v[j],*nu,1,0); mu = *theta*temp2; var=((*nu+(temp2*temp2))*(1.0-(*theta*(*theta))))/(*nu+1.0);
      hinv[j] = pt((sqrt(var)*temp1)+mu,*nu,1,0);
    }
    else if(*family==3) //clayton
    {
	  if(*theta < XEPS) hinv[j]=u[j];
	  else
      hinv[j] = pow(pow(u[j]*pow(v[j],*theta+1.0),-*theta/(*theta+1.0))+1.0-pow(v[j],-*theta),-1.0/(*theta));
    }
    else if(*family==4) //gumbel - must turn to numerical inversion
    {
	  double nu=0.0;
      HNumInv(family,&u[j],&v[j],theta,&nu,&hinv[j]);
    }
    else if(*family==5) //frank - numerical inversion
    {
	  double nu=0.0;
      HNumInv(family,&u[j],&v[j],theta,&nu,&hinv[j]);
    }
    else if(*family==6) //joe - numerical inversion
    {
	  double nu=0.0;
	  if(*theta < 5) HNumInv(family,&u[j],&v[j],theta,&nu,&hinv[j]);
	  else  
      HNumInvBisect(family,&u[j],&v[j],theta,&nu,&hinv[j]);
    }
	else if(*family==7) //BB1
	{
	  HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
	}
	else if(*family==8) //BB2
	{
	  HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
	}
	else if(*family==9) //BB7
	{
	  HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
	}
	else if(*family==13)
	  {
			u[j]=1-u[j];
			v[j]=1-v[j];
			hinv[j] = pow(pow(u[j]*pow(v[j],*theta+1.0),-*theta/(*theta+1.0))+1.0-pow(v[j],-*theta),-1.0/(*theta));
			hinv[j]=1-hinv[j];
			u[j]=1-u[j];
			v[j]=1-v[j];
	  }
	else if(*family==14) //rotated gumbel (180°) - must turn to numerical inversion
		{
			int jj=4;
			u[j]=1-u[j];
			v[j]=1-v[j];
			HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
			hinv[j]=1-hinv[j];
			u[j]=1-u[j];
			v[j]=1-v[j];
		}
	else if(*family==16) //rotated joe (180°) - must turn to numerical inversion
		{
			int jj=6;
			u[j]=1-u[j];
			v[j]=1-v[j];			
  	  if(*theta < 5) HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
  	  else  
        HNumInvBisect(&jj,&u[j],&v[j],theta,nu,&hinv[j]);			
			hinv[j]=1-hinv[j];
			u[j]=1-u[j];
			v[j]=1-v[j];
		}
	
    out[j] = MAX(MIN(hinv[j],UMAX),UMIN); //hinv[j];
  }
  Free(hinv);
}


//////////////////////////////////////////////////////////////
// Function to simulate from a pair-copula construction (vine)
// Input:
// n         sample size
// d         dimension (>= 2)
// type      vine type (1=Canonical vine, 2=D-vine)
// family    copula family (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe, 7=BB1)
// par       parameter values (at least d*(d-1)/2 parameters)
////////////////////////////////////////////////////////////////

void pcc(int* n, int* d, int* family, int* type, double* par, double* nu, double* out)
{
  int i, j, in=1, k, **fam;
  double *w, **v, **theta, **x, **ny;
  /* //Initialize random number generator: */
  GetRNGstate();
  //Allocate memory:
  w = Calloc((*d+1),double);

  v = create_matrix(*d+1,2*(*d)-1);
  theta = create_matrix(*d,*d);
  x = create_matrix(*n+1,*d+1);
  ny = create_matrix(*d,*d);
  fam = create_intmatrix(*d,*d);
  //Initialize dependency parameters
  k = 0;
  for(i=1;i<=*d-1;i++)
  {
    for(j=1;j<=*d-i;j++)
    {
      fam[i][j] = family[k];
      ny[i][j] = nu[k];
      theta[i][j] = par[k];
      k ++;
    }
  }
  //Simulate:
  if(*type==1) //Canonical vine
  {
    for(j=1;j<=*n;j++)
    {
      for(i=1;i<=*d;i++) w[i] = runif(0,1);
      v[1][1] = w[1];
      x[j][1] = v[1][1];
      for(i=2;i<=*d;i++)
      {
        v[i][1] = w[i];
        for(k=i-1;k>=1;k--)
        {
          Hinv1(&fam[k][i-k],&in, &v[i][1],&v[k][k],&theta[k][i-k],&ny[k][i-k],&v[i][1]);
        }
        x[j][i] = v[i][1];
        if(i<*d)
        {
          for(k=1;k<i;k++) Hfunc1(&fam[k][i-k],&in, &v[i][k],&v[k][k],&theta[k][i-k],&ny[k][i-k],&v[i][k+1]);
        }
      }
    }
  }
  else if(*type==2) //D-vine
  {
    for(j=1;j<=*n;j++)
    {
      for(i=1;i<=*d;i++) { w[i] = runif(0,1);} //  printf("%10.8f  \t",w[i]);} ; printf("\n");
      v[1][1] = w[1];
      v[2][1] = w[2];
      //printf("inv: %d,%d :  %5.2f : %10.8f \t",1,0, theta[1][1], v[1][1]);
      Hinv1(&fam[1][1],&in,&w[2],&v[1][1],&theta[1][1],&ny[1][1],&v[2][1]);
      //printf("%10.8f  \n",v[2][1]);
      Hfunc2(&fam[1][1],&in, &v[1][1],&v[2][1],&theta[1][1],&ny[1][1],&v[2][2]);
      for(i=3;i<=*d;i++)
      {
        v[i][1] = w[i];
	
        for(k=i-1;k>=2;k--) { 
	  //printf("inv: %d,%d :  %5.2f : %10.8f \t",k,i-k, theta[k][i-k], v[i-1][2*k-2]);
	  Hinv1(&fam[k][i-k],&in, &v[i][1],&v[i-1][2*k-2],&theta[k][i-k],&ny[k][i-k],&v[i][1]);
	  //printf("%10.8f  \n",v[i][1]);
	}
        Hinv1(&fam[1][i-1],&in, &v[i][1],&v[i-1][1],&theta[1][i-1],&ny[1][i-1],&v[i][1]);
	//printf("inv: %d,%d :  %5.2f : %10.8f \t %10.8f \n",1,i-1, theta[1][i-1], v[i-1][1],v[i][1]);
        // Compute conditional cdf's needed in next step:
        if(i<*d)
        {
          Hfunc2(&fam[1][i-1],&in, &v[i-1][1],&v[i][1],&theta[1][i-1],&ny[1][i-1],&v[i][2]);
          Hfunc1(&fam[1][i-1],&in, &v[i][1],&v[i-1][1],&theta[1][i-1],&ny[1][i-1],&v[i][3]);
          if(i>3)
          {
            for(k=2;k<=(i-2);k++)
            {
              Hfunc2(&fam[k][i-k],&in, &v[i-1][2*k-2],&v[i][2*k-1],&theta[k][i-k],&ny[k][i-k],&v[i][2*k]);
              Hfunc1(&fam[k][i-k],&in, &v[i][2*k-1],&v[i-1][2*k-2],&theta[k][i-k],&ny[k][i-k],&v[i][2*k+1]);
            }
          }
          Hfunc2(&fam[i-1][1],&in, &v[i-1][2*i-4],&v[i][2*i-3],&theta[i-1][1],&ny[i-1][1],&v[i][2*i-2]);
        }
      }
      for(i=1;i<=*d;i++) x[j][i] = v[i][1];
    }
  }
  //Write to output vector:
  k = 0;
  for(i=1;i<=*d;i++)
  {
    for(j=1;j<=*n;j++)
    {
      out[k] = x[j][i];
      k ++;
    }
  }
  PutRNGstate();
  //Free memory:
  Free(w); free_matrix(v,*d+1); free_matrix(theta,*d); free_matrix(ny,*d); free_intmatrix(fam,*d); free_matrix(x,*n+1);
}




//////////////////////////////////////////////////////////////
// Function to compute -log-likelihood for the pair-copula construction (vine)
// Input:
// n        sample size
// d        dimension (>=2)
// type     vine type (1=canonical vine, 2=d-vine)
// family   copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 7=BB1, 8=BB7)
// par      parameter values (at least d*(d-1)/2 parameters
// data     data set for which to compute log-likelihood
// Output:
// out      Loglikelihood
// ll       array with the contribution to LL (for each copula)
// vv       array for the transformation operated (Hfunc)  
/////////////////////////////////////////////////////////////
void VineLogLikm(int* T, int* d, int* type, int* family, double* par, double* data, 
		 double* out, double* ll, double* vv)
{
  int i, j, k, t, kk, **fam;
  double loglik=0.0, sumloglik=0.0, **x, **theta, **nu, ***v;
  
  //Allocate memory:
  x = create_matrix(*d+1,*T);
  // By Ulf Schepsmeier
  if(*type==1) //C-vine
	{
		v = create_3darray(*d-1,*d,*T);
	}
  else //D-vine
	{
		v = create_3darray(*d,2*(*d)-3,*T);
	}
  
  theta = create_matrix(*d,*d);
  nu = create_matrix(*d,*d);
  fam = create_intmatrix(*d+1,*d+1);
  
  //Initialize:
  k = 0;
  for(i=1;i<=*d;i++)
    { 
      for (t=0;t<=*T-1;t++ ) 
	{
	  x[i][t] = data[k];
	  k++;
	}
    }
  k = 0;
  for(i=1;i<=(*d-1);i++)
    {
      for(j=1;j<=(*d-i);j++)
		{
			theta[i][j] = par[k];
			fam[i][j] = family[k];
			nu[i][j] = par[*d*(*d-1)/2+k];
			k++;
		}
    }

  if(*type==1) //C-vine
    {
      // By Ulf Schepsmeier
	  kk=0;
	  //Compute likelihood at level 1:
		for(i=1;i<*d;i++)
		{
			LL_mod2(&fam[1][i],T,x[1],x[i+1],&theta[1][i],&nu[1][i],&loglik);
			sumloglik += loglik;
			ll[kk] = loglik; 
			++kk;
			//Compute variables for next level:
			Hfunc1(&fam[1][i],T,x[i+1],x[1],&theta[1][i],&nu[1][i],v[1][i]);
		}
		//Compute likelihood at next levels:
		for(k=2;k<=(*d-1);k++)
		{
			for(i=1;i<=(*d-k);i++)
			{
				LL_mod2(&fam[k][i],T,v[k-1][1],v[k-1][i+1],&theta[k][i],&nu[k][i],&loglik);
				sumloglik += loglik;
				ll[kk] = loglik; 
				++kk;
			}
			if(k<(*d-1))
			{
				for (i=1;i<=(*d-k);i++)
				{
					Hfunc1(&fam[k][i],T,v[k-1][i+1],v[k-1][1],&theta[k][i],&nu[k][i],v[k][i]);
				}
			}
		}
    }
  else if(*type==2) //D-vine
  {
    kk=0;
    //Compute the likelihood at level 1:
    for(i=1;i<*d;i++)
      {
        LL_mod2(&fam[1][i],T,x[i],x[i+1],&theta[1][i],&nu[1][i],&loglik);
        sumloglik += loglik;
		ll[kk] = loglik; 
		++kk;
      }
    //Compute variables for next level:
    Hfunc2(&fam[1][1],T,x[1],x[2],&theta[1][1],&nu[1][1],v[1][1]);
    for(k=1;k<=(*d-3);k++)
      {
        Hfunc1(&fam[1][k+1],T,x[k+2],x[k+1],&theta[1][k+1],&nu[1][k+1],v[1][2*k]);
        Hfunc2(&fam[1][k+1],T,x[k+1],x[k+2],&theta[1][k+1],&nu[1][k+1],v[1][2*k+1]);
      }
    Hfunc1(&fam[1][*d-1],T,x[*d],x[*d-1],&theta[1][*d-1],&nu[1][*d-1],v[1][2*(*d)-4]);
    //Compute likelihood at next levels:
    for(k=2;k<=(*d-1);k++)
      {
        for(i=1;i<=(*d-k);i++)
    	  {
  			LL_mod2(&fam[k][i],T,v[k-1][2*i-1],v[k-1][2*i],&theta[k][i],&nu[k][i],&loglik);
    	    sumloglik += loglik;
  			ll[kk] = loglik; ++kk;
    	  }
        if(k<(*d-1))
    	  {
    	    Hfunc2(&fam[k][1],T,v[k-1][1],v[k-1][2],&theta[k][1],&nu[k][1],v[k][1]);
    	    if((*d)>4)
    	      {
    			for(i=1;i<=(*d-k-2);i++)
    			{
  					Hfunc1(&fam[k][i+1],T,v[k-1][2*i+2],v[k-1][2*i+1],&theta[k][i+1],&nu[k][i+1],v[k][2*i]);
    				Hfunc2(&fam[k][i+1],T,v[k-1][2*i+1],v[k-1][2*i+2],&theta[k][i+1],&nu[k][i+1],v[k][2*i+1]);
    			}
    	      }
    	    Hfunc1(&fam[k][*d-k],T,v[k-1][2*(*d)-2*k],v[k-1][2*(*d)-2*k-1],&theta[k][*d-k],&nu[k][*d-k],v[k][2*(*d)-2*k-2]);
    	  }
      }
  }
  //Write to output:
  *out = -sumloglik;
  kk=00;
  // By Ulf Schepsmeier
	if(*type==1) //C-Vine
	{
		for(k=1;k<(*d-1);k++)
			for(i=1;i<=(*d-k);i++)
				for(t=0;t<*T;t++)
				{
					vv[kk] = v[k][i][t];
  					++kk;
  				}
		//Free memory:
		free_3darray(v,*d-1,*d);
	}
	else //D-Vine
	{
		for(k=1;k<*d;k++)
			for(i=1;i<=2*(*d-k-1);i++)
				for(t=0;t<*T;t++)
  				{
  					vv[kk] = v[k][i][t];
  					++kk;
  				}
		//Free memory:
		free_3darray(v,*d,2*(*d)-3);
	}
  //Free memory:
  free_matrix(x,*d+1); free_matrix(theta,*d); free_matrix(nu,*d); free_intmatrix(fam,*d);
}



void SimulateVine(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* out)
{
	int i, j, k, m, **fam, **cindirect, **mat, **mmat, **fam2, **cindirect2, **mat2, **mmat2;
	double **theta, **nu, **theta2, **nu2, ***vdirect, ***vindirect;
	
	//Allocate memory
	theta=create_matrix(*d,*d);
	nu=create_matrix(*d,*d);
	fam=create_intmatrix(*d,*d);
	mmat=create_intmatrix(*d,*d);
	cindirect=create_intmatrix(*d,*d);
	mat=create_intmatrix(*d,*d);
	theta2=create_matrix(*d,*d);
	nu2=create_matrix(*d,*d);
	fam2=create_intmatrix(*d,*d);
	mmat2=create_intmatrix(*d,*d);
	cindirect2=create_intmatrix(*d,*d);
	mat2=create_intmatrix(*d,*d);
	vdirect = create_3darray(*d,*d,*T);
	vindirect = create_3darray(*d,*d,*T);

	//Initialize random number generator:
	GetRNGstate();

	//Initialize
	k=0;
 	for(i=0;i<(*d);i++)
    {
        for(j=0;j<(*d);j++)
        {
            theta2[i][j]=par[(i+1)+(*d)*j-1] ;
            nu2[i][j]=par2[(i+1)+(*d)*j-1]    ;
            mmat2[i][j]=maxmat[(i+1)+(*d)*j-1] ;
            mat2[i][j]=matrix[(i+1)+(*d)*j-1] ;
            cindirect2[i][j]=conindirect[(i+1)+(*d)*j-1] ;
            fam2[i][j]=family[(i+1)+(*d)*j-1] ;
		}
    }
	
	// Matrizen rotieren für den Algo
	for(i=0;i<(*d);i++)
	{
		for(j=0;j<(*d);j++)
		{
			theta[(*d-i-1)][(*d-j-1)]=theta2[i][j];
			nu[(*d-i-1)][(*d-j-1)]=nu2[i][j];
			mmat[(*d-i-1)][(*d-j-1)]=mmat2[i][j];
			mat[(*d-i-1)][(*d-j-1)]=mat2[i][j];
			cindirect[(*d-i-1)][(*d-j-1)]=cindirect2[i][j];
			fam[(*d-i-1)][(*d-j-1)]=fam2[i][j];
		}
	}

	free_matrix(theta2,*d);
	free_matrix(nu2,*d); 
	free_intmatrix(fam2,*d);
	free_intmatrix(mmat2,*d); 
	free_intmatrix(cindirect2,*d); 
	free_intmatrix(mat2, *d);
	
	/*
	Declare variable to hold seconds on clock.
*/
//time_t seconds;
/*
Get value from system clock and
place in seconds variable.
*/
//time(&seconds);
/*
Convert seconds to a unsigned
integer.
*/
//srand((unsigned int) seconds);


	// Der eigentliche Algo
	for(j=0;j<*T;j++)
    {
		for(i=0;i<*d;i++) 	vdirect[i][i][j] = runif(0,1);
		vindirect[0][0][j] = vdirect[0][0][j];
	}

		for(i=1;i<*d;i++)
		{
			for(k=(i-1);k>(-1);k--)
			{
				m = mmat[k][i];
				if(mat[k][i]==m)
				{
					Hinv1(&fam[k][i],T,vdirect[k+1][i],vdirect[k][m-1],&theta[k][i],&nu[k][i],vdirect[k][i]);
				}
				else
				{
					Hinv1(&fam[k][i],T,vdirect[k+1][i],vindirect[k][m-1],&theta[k][i],&nu[k][i],vdirect[k][i]);
				}

				if(i+1<(*d))
				{
					if(cindirect[k+1][i]==1)
					{
						if(mat[k][i]==m)
						{
							Hfunc2(&fam[k][i],T,vdirect[k][m-1],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k+1][i]);
						}
						else
						{
							Hfunc2(&fam[k][i],T,vindirect[k][m-1],vdirect[k][i],&theta[k][i],&nu[k][i],vindirect[k+1][i]);
						}
					}
				}
			}
		}
	

	k=0;
	for(i=0;i<(*d);i++)
	{
		for(j=0;j<(*T);j++)
		{
			out[k]=vdirect[0][i][j];
			k++;
		}
	}

	//Free memory:
	free_matrix(theta,*d);
	free_matrix(nu,*d); 
	free_intmatrix(fam,*d);
	free_intmatrix(mmat,*d); 
	free_intmatrix(cindirect,*d); 
	free_intmatrix(mat, *d);
	free_3darray(vdirect,*d,*d);
	free_3darray(vindirect,*d,*d);
	PutRNGstate();
}



void ktau(double *X, double *Y, int *N, double *tau, double *S, double *D, int *T, int *U, int *V)
{
  // Defining variables
  int K, L, I, J, Iend, Jend;
  int i, j, m;
  //double *Y2=(double *)malloc((*N)*sizeof(double));
  double *Y2 = (double*) Calloc(*N, double);
  //double *X2=(double *)malloc((*N)*sizeof(double));
  double *X2 = (double*) Calloc(*N, double);
  double *xptr,*yptr; // HJ addition for swapping
  boolean Iflag, Jflag, Xflag;
  *S = 0.; *D = 0.; *T = 0; *U = 0; *V = 0;

  /* 1.1 Sort X and Y in X order */
  /* Break ties in X according to Y */
  K=1;
  do
  {
    L=0;
    do
    {
      I = L;
      J = (I+K)<(*N)?(I+K):(*N);
      Iend = J;
      Jend = (J+K)<(*N)?(J+K):(*N);
      do
      {
	Iflag = (I < Iend);
	Jflag = (J < Jend);
	Xflag = ((X[I] > X[J]) | ((X[I] == X[J]) & (Y[I] > Y[J])));
	if((Iflag & !Jflag) | (Iflag & Jflag & !Xflag))
	{
	  X2[L] = X[I];
	  Y2[L] = Y[I];
	  I++;
	  L++;
	};
	if((!Iflag & Jflag) | (Iflag & Jflag & Xflag))
	{
	  X2[L] = X[J];
	  Y2[L] = Y[J];
	  J++;
	  L++;
	};
      } while(Iflag | Jflag);
    } while(L < *N);
    // Swap lists
    xptr=X; X=X2; X2=xptr;
    yptr=Y; Y=Y2; Y2=yptr;
#ifdef OLD
    for(i = 0; i < *N; i++)
    { Xtem = X[i]; Ytem = Y[i];
      X[i] = X2[i]; Y[i] = Y2[i];
      X2[i] = Xtem; Y2[i] = Ytem;
    };
#endif
    K *= 2;
  } while (K < *N);

  /* 1.2 Count pairs of tied X, T */
  j = 1;
  m = 1;
  for(i = 1; i < *N; i++)
    if(X[i] == X[i-1])
    {
      j++;
      if(Y[i] == Y[i-1])
	m++;
    }
    else if(j > 1)
    {
      *T += j * (j - 1) / 2;
      if(m > 1)
	*V += m * (m - 1) / 2;
      j = 1;
      m = 1;
    };
  *T += j * (j - 1) / 2;
  *V += m * (m - 1) / 2;

  /* 2.1 Sort Y again and count exchanges, S */
  /* Keep original relative order if tied */

  K=1;
  do
  {
    L=0;
    do
    {
      I = L;
      J = (I+K)<(*N)?(I+K):(*N);
      Iend = J;
      Jend = (J+K)<(*N)?(J+K):(*N);
      do
      {
	Iflag = (I < Iend);
	Jflag = (J < Jend);
	Xflag = (Y[I] > Y[J]);
	if((Iflag & !Jflag) | (Iflag & Jflag & !Xflag))
	{
	  X2[L] = X[I];
	  Y2[L] = Y[I];
	  I++;
	  L++;
	};
	if((!Iflag & Jflag) | (Iflag & Jflag & Xflag))
	{
	  X2[L] = X[J];
	  Y2[L] = Y[J];
	  *S += Iend - I;
	  J++;
	  L++;
	};
      } while(Iflag | Jflag);
    } while(L < *N);
    
    // Swap lists
    xptr=X; X=X2; X2=xptr;
    yptr=Y; Y=Y2; Y2=yptr;
#ifdef OLD
    for(i = 0; i < *N; i++)
    { Xtem = X[i]; Ytem = Y[i];
      X[i] = X2[i]; Y[i] = Y2[i];
      X2[i] = Xtem; Y2[i] = Ytem;
    };
#endif
    K *= 2;
  } while (K < *N);

  /* 2.2 Count pairs of tied Y, U */
  j=1;
  for(i = 1; i < *N; i++)
    if(Y[i] == Y[i-1])
      j++;
    else if(j > 1)
    {
      *U += j * (j - 1) / 2;
      j = 1;
    };
  *U += j * (j - 1) / 2;


  /* 3. Calc. Kendall's Score and Denominator */
  *D = 0.5 * *N * (*N - 1);
  *S = *D - (2. * *S + *T + *U - *V);
  //if(*T > 0 | *U > 0) // adjust for ties
    *D = sqrt((*D - *T) * (*D - *U));
  *tau = *S / *D;


  Free(Y2);
  Free(X2);
}
