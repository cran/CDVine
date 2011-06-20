/*
** hfunc.c - C code of the package CDRVine  
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


// h-func for BB1

void pcondbb1(double* u, double* v, int* n, double* param, double* out)
{
	int i;
	double th, de;
	double *t1, *t2, *t3, *t16, *t17, *t4, *t5, *t6, *t7, *t9, *t10, *t12, *t13, *t20;
	t1 = Calloc(*n,double);
	t2 = Calloc(*n,double);
	t3 = Calloc(*n,double);
	t16 = Calloc(*n,double);
	t17 = Calloc(*n,double);
	t4 = Calloc(*n,double);
	t5 = Calloc(*n,double);
	t6 = Calloc(*n,double);
	t7 = Calloc(*n,double);
	t9 = Calloc(*n,double);
	t10 = Calloc(*n,double);
	t12 = Calloc(*n,double);
	t13 = Calloc(*n,double);
	t20 = Calloc(*n,double);

	th = param[0];
	de = param[1];
	for(i=0;i<*n;i++)
		{
		  t1[i] = pow(u[i],-th);
		  t2[i] = t1[i]-1.;
		  t3[i] = pow(t2[i],de);
		  t16[i] = 1./u[i];
		  t17[i] = 1./t2[i];
		  t4[i] = pow(v[i],-th);
		  t5[i] = t4[i]-1.;
		  t6[i] = pow(t5[i],de);
		  t7[i] = t3[i]+t6[i];
		  t9[i] = pow(t7[i],1/de);
		  t10[i] = 1.0+t9[i];
		  t12[i] = pow(t10[i],-1/th);
		  t13[i] = t12[i]*t9[i];
		  t20[i] = 1./t10[i];
		  out[i] = t13[i]*t3[i]*t1[i]*t16[i]*t17[i]/t7[i]*t20[i];
		}

  Free(t1);Free(t2);Free(t3);Free(t16);Free(t17);Free(t4);Free(t5);Free(t6);Free(t7);
  Free(t9);Free(t10);Free(t13);Free(t20);
}



void pcondbb6(double* u, double* v, int* n, double* param, double* out)
{
	int i;
	double th, de;
	double *t1, *t2, *t3, *t4, *t5, *t12, *t16, *t6, *t7, *t8,*t9, *t10, *t11, *t13, *t14, *t15, *t17;
	t1 = Calloc(*n,double);
	t2 = Calloc(*n,double);
	t3 = Calloc(*n,double);
	t4 = Calloc(*n,double);
	t5 = Calloc(*n,double);
	t12 = Calloc(*n,double);
	t16 = Calloc(*n,double);
	t6 = Calloc(*n,double);
	t7 = Calloc(*n,double);
	t8 = Calloc(*n,double);
	t9 = Calloc(*n,double);
	t10 = Calloc(*n,double);
	t11 = Calloc(*n,double);
	t13 = Calloc(*n,double);
	t14 = Calloc(*n,double);
	t15 = Calloc(*n,double);
	t17 = Calloc(*n,double);

	th = param[0];
	de = param[1];

	for(i=0;i<*n;i++)
	{
		  t1[i] = 1.0-u[i];
		  t2[i] = pow(t1[i],th);
		  t3[i] = 1.0-t2[i];
		  t4[i] = log(t3[i]);
		  t5[i] = pow(-t4[i],de);
		  t12[i] = 1/de;
		  t16[i] = 1/th;
		  t6[i] = 1.0-v[i];
		  t7[i] = pow(t6[i],th);
		  t8[i] = 1.0-t7[i];
		  t9[i] = log(t8[i]);
		  t10[i] = pow(-t9[i],de);
		  t11[i] = t5[i]+t10[i];
		  t13[i] = pow(t11[i],t12[i]);
		  t14[i] = exp(-t13[i]);
		  t15[i] = 1.0-t14[i];
		  t17[i] = pow(t15[i],t16[i]);

		  out[i] = -t17[i]*t13[i]*t5[i]*t2[i]/t1[i]/t3[i]/t4[i]/t11[i]*t14[i]/t15[i];
	}

	Free(t1);Free(t2);Free(t3);Free(t4);Free(t5);Free(t12);Free(t16);
	Free(t6);Free(t7);Free(t8);Free(t9);Free(t10);
	Free(t11);Free(t13);Free(t14);Free(t15);Free(t17);
}


void pcondbb7(double* u, double* v, int* n, double* param, double* out)
{
	int i;
	double th, de;
	double *t1, *t2, *t3, *t4, *t6, *t8, *t9, *t11, *t12, *t13, *t14;
	t1 = Calloc(*n,double);
	t2 = Calloc(*n,double);
	t3 = Calloc(*n,double);
	t4 = Calloc(*n,double);
	t6 = Calloc(*n,double);
	t8 = Calloc(*n,double);
	t9 = Calloc(*n,double);
	t11 = Calloc(*n,double);
	t12 = Calloc(*n,double);
	t13 = Calloc(*n,double);
	t14 = Calloc(*n,double);

	th = param[0];
	de = param[1];

	for(i=0;i<*n;i++)
	{
		  t1[i] = 1.0-u[i];
		  t2[i] = pow(t1[i],1.0*th);
		  t3[i] = 1.0-t2[i];
		  t4[i] = pow(t3[i],-1.0*de);
		  t6[i] = pow(1.0-v[i],1.0*th);
		  t8[i] = pow(1.0-t6[i],-1.0*de);
		  t9[i] = t4[i]+t8[i]-1.0;
		  t11[i] = pow(t9[i],-1.0/de);
		  t12[i] = 1.0-t11[i];
		  t14[i] = pow(t12[i],1.0/th);

		  out[i] = t14[i]*t11[i]*t4[i]*t2[i]/t1[i]/t3[i]/t9[i]/t12[i];
	}

	Free(t1);Free(t2);Free(t3);Free(t6);
	Free(t8);Free(t9);Free(t11);Free(t12);
	Free(t13);Free(t14);
}


void pcondbb8(double* u, double* v, int* n, double* param, double* out)
{
	int i;
	double th, de;
	double *t2, *t3, *t12, *t16, *t6, *t7, *t8, *t10, *t11, *t13, *t15, *t17;
	t2 = Calloc(*n,double);
	t3 = Calloc(*n,double);
	t12 = Calloc(*n,double);
	t16 = Calloc(*n,double);
	t6 = Calloc(*n,double);
	t7 = Calloc(*n,double);
	t8 = Calloc(*n,double);
	t10 = Calloc(*n,double);
	t11 = Calloc(*n,double);
	t13 = Calloc(*n,double);
	t15 = Calloc(*n,double);
	t17 = Calloc(*n,double);

	th = param[0];
	de = param[1];

	for(i=0;i<*n;i++)
	{
		  t2[i] = 1.0-de*u[i];
		  t3[i] = pow(t2[i],th);
		  t10[i] = 1.0-de;
		  t11[i] = pow(t10[i],th);
		  t12[i] = 1.0-t11[i];
		  t13[i] = 1/t12[i];
		  t16[i] = 1/th;
		  t6[i] = 1.0-de*v[i];
		  t7[i] = pow(t6[i],th);
		  t8[i] = 1.0-t7[i];
		  t15[i] = 1.0-(1.0-t3[i])*t8[i]*t13[i];
		  t17[i] = pow(t15[i],t16[i]);

		  out[i] = t17[i]*t3[i]/t2[i]*t8[i]*t13[i]/t15[i];
	}

	Free(t2);Free(t3);Free(t12);Free(t16);
	Free(t6);Free(t7);Free(t8);Free(t10);
	Free(t11);Free(t13);Free(t15);Free(t17);
}



// Since the h function is not symmetric in case of double Gumbel and double Clayton we have two implement both separately,
// i.e. Hfunc1 and Hfunc2
void  Hfunc1(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out)
{
  double *negv, *negu;
  negv = (double *) malloc(*n* sizeof(double));
  negu = (double *) malloc(*n*sizeof(double));
  double ntheta, nnu;
  int nfamily;
  ntheta = -*theta;
  nnu = -*nu;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30)))
    {
	  nfamily=(*family)-20;
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc (&nfamily, n, u, negv, &ntheta, &nnu, out);
    }
  else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40)))
	{
	  nfamily=(*family)-30;
      for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
      Hfunc(&nfamily, n, negu, v, &ntheta, &nnu, out);
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
  double ntheta, nnu;
  int nfamily;
  ntheta = -*theta;
  nnu = -*nu;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

  if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30)))
    {
	  nfamily=(*family)-20;
      for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
      Hfunc(&nfamily, n, negv, u, &ntheta, &nnu, out);
		for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
    }
  else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40)))
	{
	  nfamily=(*family)-30;
      for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
      Hfunc(&nfamily, n, v, negu, &ntheta, &nnu, out);
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
		    h[j] =   pow(v[j],-*theta-1.0)*pow(x,-1.0-1.0/(*theta));
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
			else if(*theta==0) 
				{
				h[j]=-(exp(-pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)))*pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(v[j]),*nu))/(v[j]*log(v[j]));
				}
			else
			{
			pcondbb1(&v[j],&u[j],&T,param,&h[j]);
			}
			Free(param);
		}
		else if(*family==8) //BB6
		{
			double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
      if(*theta==1) 
				{
				if(*nu==1) h[j]=u[j];
				else h[j]=-(exp(-pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)))*pow(pow(-log(v[j]),*nu)+pow(-log(u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(v[j]),*nu))/(v[j]*log(v[j]));
				}     
      else if(*nu==1) 
				{
				h[j]=pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
    	}
			else
			{
			pcondbb6(&v[j],&u[j],&T,param,&h[j]);
			}
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
			else if(*nu==0)
			{
			  h[j] = pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
			}
			else
			{
			pcondbb7(&v[j],&u[j],&T,param,&h[j]);
			}
			Free(param);
		}
		else if(*family==10) //BB8
		  {
			double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
      if(*nu==0)
			{
				h[j]=u[j];
			}			
      else if(*nu==1)
			{
				if(*theta==1) h[j]=u[j];
				else h[j]=pow(pow(1.0-u[j],*theta) + pow(1.0-v[j],*theta) - pow(1.0-u[j],*theta)*pow(1.0-v[j],*theta),1.0/(*theta)-1) * pow(1.0-v[j],*theta-1.0)*(1-pow(1-u[j],*theta));
    	}
			else
			{
			pcondbb8(&v[j],&u[j],&T,param,&h[j]);
			}
			Free(param);		  
		  }
		else if(*family==13) //rotated clayton (180°)
		  {
			if(*theta == 0) h[j] = u[j] ;
			if(*theta < XEPS) h[j] = u[j] ;
			else
			{
				u[j]=1-u[j]; 
				v[j]=1-v[j];
		    x = pow(u[j],-*theta)+pow(v[j],-*theta)-1.0 ;
		    h[j] =   pow(v[j],-*theta-1.0)*pow(x,-1.0-1.0/(*theta)); // pow(v[j],-*theta-1.0)*pow(pow(u[j],-*theta)+pow(v[j],-*theta)-1.0,-1.0-1.0/(*theta));
				h[j]= 1-h[j];
				u[j]=1-u[j];
				v[j]=1-v[j];
		  }
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
		else if(*family==17) //rotated BB1
		  {
			  double* param;
				param = Calloc(2,double);
				param[0]=*theta;
				param[1]=*nu;
				int T=1;
				if(*nu==1) 
				{
				if(*theta==0) h[j]=u[j];
				else
          { 
            h[j]=pow(pow(1-u[j],-*theta)+pow(1-v[j],-*theta)-1,-1/(*theta)-1)*pow(1-v[j],-*theta-1);
				    h[j]= 1-h[j];
          }
        }
	   else if(*theta==0) 
				{
				h[j]=-(exp(-pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)))*pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(1-v[j]),*nu))/((1-v[j])*log(1-v[j]));
				h[j]= 1-h[j];
        }
				else
				{
					v[j]= 1-v[j];
					u[j]= 1-u[j];
					pcondbb1(&v[j],&u[j],&T,param,&h[j]);
					u[j]=1-u[j];
					v[j]=1-v[j];
					h[j]= 1-h[j];
				}
				Free(param);
		  }
		else if(*family==18) //rotated BB6
		  {
			double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
      if(*theta==1) 
				{
				if(*nu==1) h[j]=u[j];
				else
          {
            h[j]=-(exp(-pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)))*pow(pow(-log(1-v[j]),*nu)+pow(-log(1-u[j]),*nu),1.0/(*nu)-1.0)*pow(-log(1-v[j]),*nu))/((1-v[j])*log(1-v[j]));
				    h[j]= 1-h[j];
          }
        }     
      else if(*nu==1) 
			{
				h[j]=pow(pow(u[j],*theta) + pow(v[j],*theta) - pow(u[j],*theta)*pow(v[j],*theta),1.0/(*theta)-1) * pow(v[j],*theta-1.0)*(1-pow(u[j],*theta));
    	  h[j]= 1-h[j];
      }
			else
			{
          v[j]= 1-v[j];
					u[j]= 1-u[j];
					pcondbb6(&v[j],&u[j],&T,param,&h[j]);
					u[j]=1-u[j];
					v[j]=1-v[j];
					h[j]= 1-h[j];	
			}
			Free(param);	  		  
		  }
		else if(*family==19) //rotated BB7
		  {
				double* param;
				param = Calloc(2,double);
				param[0]=*theta;
				param[1]=*nu;
				int T=1;
				if(*theta==1)
				{
					if(*nu==0) h[j]=u[j];
					else{
            h[j]=pow(pow(1-u[j],-*nu)+pow(1-v[j],-*nu)-1,-1/(*nu)-1)*pow(1-v[j],-*nu-1);
            h[j]= 1-h[j];
          }
				}
			else if(*nu==0)
			{
			  h[j] = pow(pow(u[j],*theta) + pow(v[j],*theta) - pow(u[j],*theta)*pow(v[j],*theta),1.0/(*theta)-1) * pow(v[j],*theta-1.0)*(1-pow(u[j],*theta));
			  h[j]= 1-h[j];
      }				
				else
				{
					v[j]= 1-v[j];
					u[j]= 1-u[j];
					pcondbb7(&v[j],&u[j],&T,param,&h[j]);
					u[j]=1-u[j];
					v[j]=1-v[j];
					h[j]= 1-h[j];
				}
				Free(param);
		  }
		else if(*family==20) //rotated BB8
		  {
			double* param;
			param = Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			int T=1;
      if(*nu==0)
			{
				h[j]=u[j];
			}			
      else if(*nu==1)
			{
				if(*theta==1) h[j]=u[j];
				else{
          h[j]=pow(pow(u[j],*theta) + pow(v[j],*theta) - pow(u[j],*theta)*pow(v[j],*theta),1.0/(*theta)-1) * pow(v[j],*theta-1.0)*(1-pow(u[j],*theta));
    	    h[j]= 1-h[j];
        }
      }
			else
			{
					v[j]= 1-v[j];
					u[j]= 1-u[j];
          pcondbb8(&v[j],&u[j],&T,param,&h[j]);
					u[j]=1-u[j];
					v[j]=1-v[j];
					h[j]= 1-h[j];
			}
			Free(param);		  
		  }	
    }	
    out[j] = MAX(MIN(h[j],UMAX),UMIN);
  }	
  Free(h);
}




///////////////////////////////////////////////////////////////
// Function to compute inversion of H numerically through Bisection
// Input:
// u        variable for which h-function computes conditional distribution function
// v        variable on which h-function conditions
// theta    parameter for the copula family
// out      output
//////////////////////////////////////////////////////////////
void HNumInv(int* family, double* u, double* v, double* theta, double* nu, double* out)
{

  int br=0, in=1;
  double ans=0.0, tol=0.000001, x0=UMIN, x1=UMAX, fl, fh, val;
  
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

  			if(val > 0.0) {x1 = ans; fh = val;}
  			else {x0 = ans; fl = val;}
  
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
  double ntheta, nnu;
  int nfamily;
  ntheta = -*theta;
  nnu = -*nu;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

   if(((*family ==23) | (*family ==24) | (*family==26) | (*family ==27) | (*family ==28) | (*family==29) | (*family==30)))
      {
	   nfamily=(*family)-20;
		for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
		Hinv(&nfamily,  n,  u,  negv,  &ntheta,  &nnu,  out);
      }
	else if(((*family==33) | (*family==34) | (*family==36) | (*family ==37) | (*family ==38) | (*family==39) | (*family==40)))
      {
	   nfamily=(*family)-30;
		for (int i = 0; i < *n; i++) {negu[i]=1 - u[i];};
		Hinv(&nfamily,  n,  negu,  v,  &ntheta,  &nnu,  out);
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
  double ntheta, nnu;
  int nfamily;
  ntheta = -*theta;
  nnu = -*nu;

	for(int i=0;i<*n;i++)
	{
		if(u[i]<UMIN) u[i]=UMIN;
		else if(u[i]>UMAX) u[i]=UMAX;
		if(v[i]<UMIN) v[i]=UMIN;
		else if(v[i]>UMAX) v[i]=UMAX;
	}

   if(((*family ==23) | (*family ==24) | (*family==26) | (*family ==27) | (*family ==28) | (*family==29) | (*family==30)))
      {
	   nfamily = (*family)-20;
		for (int i = 0; i < *n; ++i) {negv[i]=1 - v[i];}
		Hinv(&nfamily,  n,  negv, u,  &ntheta,  &nnu,  out);
		for (int i = 0; i < *n; i++) {out[i]=1-out[i];};
      }
	else if(((*family==33) | (*family==34) | (*family==36) | (*family ==37) | (*family ==38) | (*family==39) | (*family==40)))
      {
	   nfamily=(*family)-30;
		for (int i = 0; i < *n; ++i) {negu[i]=1 - u[i];}
		Hinv(&nfamily,  n,  v,  negu,  &ntheta,  &nnu,  out);
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
  int j;
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
	if(*family==0)
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
	  HNumInv(family,&u[j],&v[j],theta,&nu,&hinv[j]);
	  }
	else if(*family==7) //BB1
	{
	  HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
	}
	else if(*family==8) //BB6
	{
	  HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
	}
	else if(*family==9) //BB7
	{
	  HNumInv(family,&u[j],&v[j],theta,nu,&hinv[j]);
	}
	else if(*family==10) //BB8
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
  	  HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);			
			hinv[j]=1-hinv[j];
			u[j]=1-u[j];
			v[j]=1-v[j];
		}
	else if(*family==17) //rotated BB1 (180°) - must turn to numerical inversion
	{
			int jj=7;
			u[j]=1-u[j];
			v[j]=1-v[j];
			HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
			hinv[j]=1-hinv[j];
			u[j]=1-u[j];
			v[j]=1-v[j];
	}
	else if(*family==18) //rotated BB6 (180°) - must turn to numerical inversion
	{
			int jj=8;
			u[j]=1-u[j];
			v[j]=1-v[j];
			HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
			hinv[j]=1-hinv[j];
			u[j]=1-u[j];
			v[j]=1-v[j];
	}
	else if(*family==19) //rotated BB7 (180°) - must turn to numerical inversion
	{
			int jj=9;
			u[j]=1-u[j];
			v[j]=1-v[j];
			HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
			hinv[j]=1-hinv[j];
			u[j]=1-u[j];
			v[j]=1-v[j];
	}
	else if(*family==20) //rotated BB8 (180°) - must turn to numerical inversion
	{
			int jj=10;
			u[j]=1-u[j];
			v[j]=1-v[j];
			HNumInv(&jj,&u[j],&v[j],theta,nu,&hinv[j]);
			hinv[j]=1-hinv[j];
			u[j]=1-u[j];
			v[j]=1-v[j];
	}
	
    out[j] = MAX(MIN(hinv[j],UMAX),UMIN); 
  }
  Free(hinv);
}

