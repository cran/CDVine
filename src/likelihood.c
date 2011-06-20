/*
** likelihood.c - C code of the package CDRVine  
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


void archCDF(double* u, double* v, int* n, double* param, int* copula, double* out)
{
	int j;
	double *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *t9, *t10, *t11, *t12, *t13, *t14;
	t1 = Calloc(*n,double);
	t2 = Calloc(*n,double);
	t3 = Calloc(*n,double);
	t4 = Calloc(*n,double);
	t5 = Calloc(*n,double);
	t6 = Calloc(*n,double);
	t7 = Calloc(*n,double);
	t8 = Calloc(*n,double);
	t9 = Calloc(*n,double);
	t10 = Calloc(*n,double);
	t11 = Calloc(*n,double);
	t12 = Calloc(*n,double);
	t13 = Calloc(*n,double);
	t14 = Calloc(*n,double);

	for(j=0;j<*n;j++)
	{
		if(u[j]>UMAX && v[j]>UMAX){ out[j]=1;}
		else if(u[j]>UMAX){ out[j]=v[j];}
		else if(v[j]>UMAX){ out[j]=u[j];}
		else if(u[j]<UMIN || v[j]<UMIN){ out[j]=0;}
		else
		{
		if(*copula==3)	//Clayton
		{
			t1[j] = pow(u[j],-param[0]);
			t2[j] = pow(v[j],-param[0]);
			t3[j] = t1[j]+t2[j]-1;
			out[j] = pow(t3[j],-1/param[0]);
		}
		else if(*copula==4)	//Gumbel
		{
			t1[j] = -log(u[j]);
			t2[j] = -log(v[j]);
			t3[j] = pow(t1[j],param[0]);
			t4[j] = pow(t2[j],param[0]);
			t5[j] = t3[j]+t4[j];
			t6[j] = -pow(t5[j],1/param[0]);
			out[j] = exp(t6[j]);
		}
		else if(*copula==5)	//Frank
		{
			double nu;
			t1[j] = -param[0]*u[j];
			t2[j] = -param[0]*v[j];
			t3[j] = exp(t1[j]);
			t4[j] = exp(t2[j]);
			t5[j] = 1-t3[j];
			t6[j] = 1-t4[j];
			nu = 1-exp(-param[0]);
			t7[j] = t5[j]*t6[j];
			t8[j] = nu-t7[j];
			out[j] = -1/param[0]*log(t8[j]/nu);
		}
		else if(*copula==6)	//Joe
		{
			t1[j] = 1-u[j];
			t2[j] = 1-v[j];
			t3[j] = pow(t1[j],param[0]);
			t4[j] = pow(t2[j],param[0]);
			t5[j] = t3[j]*t4[j];
			out[j] = 1-pow(t3[j]+t4[j]-t5[j],1/param[0]);
		}
		else if(*copula==7)	//BB1
		{
			t1[j] = pow(u[j],-param[0]);
			t2[j] = pow(v[j],-param[0]);
			t3[j] = t1[j]-1;
			t4[j] = t2[j]-1;
			t5[j] = pow(t3[j],param[1]);
			t6[j] = pow(t4[j],param[1]);
			t7[j] = t5[j]+t6[j];
			t8[j] = pow(t7[j],1/param[1]);
			out[j] = pow(1+t8[j],-1/param[0]);
		}
		else if(*copula==8)	//BB6
		{
			t1[j] = 1-u[j];
			t2[j] = 1-v[j];
			t3[j] = pow(t1[j],param[0]);
			t4[j] = pow(t2[j],param[0]);
			t5[j] = 1-t3[j];
			t6[j] = 1-t4[j];
			t7[j] = -log(t5[j]);
			t8[j] = -log(t6[j]);
			t9[j] = pow(t7[j],param[1]);
			t10[j] = pow(t8[j],param[1]);
			t11[j] = t9[j]+t10[j];
			t12[j] = pow(t11[j],1/param[1]);
			t13[j] = exp(-t12[j]);
			t14[j] = 1-t13[j];
			out[j] = 1-pow(t14[j],1/param[0]);
		}
		else if(*copula==9)	//BB7
		{
			t1[j] = 1-u[j];
			t2[j] = 1-v[j];
			t3[j] = pow(t1[j],param[0]);
			t4[j] = pow(t2[j],param[0]);
			t5[j] = 1-t3[j];
			t6[j] = 1-t4[j];
			t7[j] = pow(t5[j],-param[1]);
			t8[j] = pow(t6[j],-param[1]);
			t9[j] = t7[j]+t8[j]-1;
			t10[j] = pow(t9[j],-1/param[1]);
			t11[j] = 1-t10[j];
			t12[j] = pow(t11[j],1/param[0]);
			out[j] = 1-t12[j];
		}
		else if(*copula==10)    //BB8
        {
            double nu;
            t1[j] = param[1]*u[j];
            t2[j] = param[1]*v[j];
            t3[j] = 1-t1[j];
            t4[j] = 1-t2[j];
            t5[j] = pow(t3[j],param[0]);
            t6[j] = pow(t4[j],param[0]);
            t7[j] = 1-t5[j];
            t8[j] = 1-t6[j];
            nu = 1-param[1];
            nu = pow(nu,param[0]);
            nu = 1-nu;
			nu = 1/nu;
            t9[j] = 1-nu*t7[j]*t8[j];
            t10[j] = pow(t9[j],1/param[0]);
            out[j] = 1/param[1]*(1-t10[j]);
        }

		}
	}
	Free(t1); Free(t2); Free(t3); Free(t4); Free(t5); Free(t6); Free(t7); Free(t8);
	Free(t9); Free(t10); Free(t11); Free(t12); Free(t13); Free(t14);
}



void dbb1(double* u, double* v, int* n, double* param, double* out)
{
	int i;
	double th, de;
	double *t1, *t2, *t3, *t16, *t17, *t38, *t39, *t4, *t5, *t6, *t7, *t8,*t9, *t10, *t12, *t13, *t20, *t24, *t25, *t27, *t29, *t32, *t33, *t34, *t36, *t43, *t59;
	t1 = Calloc(*n,double);
	t2 = Calloc(*n,double);
	t3 = Calloc(*n,double);
	t16 = Calloc(*n,double);
	t17 = Calloc(*n,double);
	t38 = Calloc(*n,double);
	t39 = Calloc(*n,double);
	t4 = Calloc(*n,double);
	t5 = Calloc(*n,double);
	t6 = Calloc(*n,double);
	t7 = Calloc(*n,double);
	t8 = Calloc(*n,double);
	t9 = Calloc(*n,double);
	t10 = Calloc(*n,double);
	t12 = Calloc(*n,double);
	t13 = Calloc(*n,double);
	t20 = Calloc(*n,double);
	t24 = Calloc(*n,double);
	t25 = Calloc(*n,double);
	t27 = Calloc(*n,double);
	t29 = Calloc(*n,double);
	t32 = Calloc(*n,double);
	t33 = Calloc(*n,double);
	t34 = Calloc(*n,double);
	t36 = Calloc(*n,double);
	t43 = Calloc(*n,double);
	t59 = Calloc(*n,double);

	th = param[0];
	de = param[1];

	for(i=0;i<*n;i++)
	{
		t1[i] = pow(u[i],(-th));
		t2[i] = t1[i]-1.0;
		t3[i] = pow(t2[i],de);
		t16[i] = 1./u[i];
		t17[i] = 1./t2[i];
		t38[i] = t1[i]*t16[i];
		t39[i] = t38[i]*t17[i];
		t4[i] = pow(v[i],(-th));
		t5[i] = t4[i]-1.0;
		t6[i] = pow(t5[i],de);
		t7[i] = t3[i]+t6[i];
		t9[i] = pow(t7[i],(1./de));
		t10[i] = 1.0+t9[i];
		t12[i] = pow(t10[i],(-1./th));
		t13[i] = t12[i]*t9[i];
		t20[i] = 1./t10[i];
		t24[i] = t9[i]*t9[i];
		t25[i] = t12[i]*t24[i];
		t27[i] = 1./v[i];
		t29[i] = 1./t5[i];
		t32[i] = t7[i]*t7[i];
		t33[i] = 1./t32[i];
		t34[i] = t10[i]*t10[i];
		t36[i] = t33[i]/t34[i];
		t43[i] = t4[i]*th;
		t59[i] = t43[i]*t27[i]*t29[i];

		out[i] = t25[i]*t6[i]*t27[i]*t4[i]*t29[i]*t36[i]*t3[i]*t39[i]-t13[i]*t6[i]*t43[i]*t27[i]*t29[i]*t33[i]*t3[i]*t38[i]*t17[i]*t20[i]+
        t13[i]*t3[i]*t38[i]*t17[i]*t33[i]*t20[i]*t6[i]*de*t59[i]+t25[i]*t3[i]*t39[i]*t36[i]*t6[i]*t59[i];
	}

	Free(t1);Free(t2);Free(t3);Free(t16);Free(t17);Free(t38);Free(t39);Free(t4);Free(t5);Free(t6);
	Free(t7);Free(t8);Free(t9);Free(t10);Free(t12);Free(t13);Free(t20);Free(t24);Free(t25);Free(t27);
	Free(t29);Free(t32);Free(t33);Free(t34);Free(t36);Free(t43);Free(t59);
}


void dbb6(double* u, double* v, int* n, double* param, double* out)
{
	int i;
	double th, de;
	double *t1, *t2, *t3, *t4, *t5, *t12, *t16, *t32, *t38, *t39, *t40, *t47, *t50, *t61, *t90, *t6, *t7, *t8,*t9, *t10, *t11, *t13, *t14, *t15, *t17, *t35, *t36, *t37, *t42, *t48, *t53, *t56, *t57, *t59, *t78, *t80, *t87, *t93;
	t1 = Calloc(*n,double);
	t2 = Calloc(*n,double);
	t3 = Calloc(*n,double);
	t4 = Calloc(*n,double);
	t5 = Calloc(*n,double);
	t12 = Calloc(*n,double);
	t16 = Calloc(*n,double);
	t32 = Calloc(*n,double);
	t38 = Calloc(*n,double);
	t39 = Calloc(*n,double);
	t40 = Calloc(*n,double);
	t47 = Calloc(*n,double);
	t50 = Calloc(*n,double);
	t61 = Calloc(*n,double);
	t90 = Calloc(*n,double);
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
	t35 = Calloc(*n,double);
	t36 = Calloc(*n,double);
	t37 = Calloc(*n,double);
	t42 = Calloc(*n,double);
	t48 = Calloc(*n,double);
	t53 = Calloc(*n,double);
	t56 = Calloc(*n,double);
	t57 = Calloc(*n,double);
	t59 = Calloc(*n,double);
	t78 = Calloc(*n,double);
	t80 = Calloc(*n,double);
	t87 = Calloc(*n,double);
	t93 = Calloc(*n,double);

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
		  t32[i] = de-1.0;
		  t38[i] = 2.0*de;
		  t39[i] = -1.0+t38[i];
		  t40[i] = pow(-t4[i],t39[i]);
		  t47[i] = 3.0*de-1.0;
		  t50[i] = pow(-t4[i],t32[i]);
		  t61[i] = pow(-t4[i],t47[i]);
		  t90[i] = pow(-t4[i],t38[i]);
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
		  t35[i] = pow(t11[i],-2.0*t32[i]*t12[i]);
		  t36[i] = t35[i]*th;
		  t37[i] = exp(t13[i]);
		  t42[i] = pow(-t9[i],t39[i]);
		  t48[i] = pow(-t9[i],t47[i]);
		  t53[i] = t13[i]*de;
		  t56[i] = pow(-t9[i],t32[i]);
		  t57[i] = t37[i]*t50[i]*t56[i];
		  t59[i] = t13[i]*th;
		  t78[i] = t37[i]-1.0;
		  t80[i] = pow(t78[i]*t14[i],t16[i]);
		  t87[i] = t78[i]*t78[i];
		  t93[i] = pow(-t9[i],t38[i]);

		  out[i] = (2.0*t36[i]*t37[i]*t40[i]*t42[i]+t36[i]*t37[i]*t48[i]*t50[i]+t53[i]*th*t57[i]-t59[i]*t57[i]+
					t36[i]* t37[i]*t61[i]*t56[i]-2.0*t35[i]*t40[i]*t42[i]-t35[i]*t61[i]*t56[i]-t53[i]*th*t50[i]*t56[i]+t59[i]*t50[i]*t56[i]-
					t35[i]*t48[i]*t50[i]) *t80[i]*t7[i]*t2[i]/t3[i]/t8[i]/t87[i]/(t90[i]+2.0*t5[i]*t10[i]+t93[i])/t1[i]/t6[i];
	}

	Free(t1);Free(t2);Free(t3);Free(t4);Free(t5);Free(t12);Free(t16);Free(t32);Free(t38);Free(t39);
	Free(t40);Free(t47);Free(t50);Free(t61);Free(t90);Free(t6);Free(t7);Free(t8);Free(t9);Free(t10);
	Free(t11);Free(t13);Free(t14);Free(t15);Free(t17);Free(t35);Free(t36);Free(t37);Free(t42);Free(t48);
	Free(t53);Free(t56);Free(t57);Free(t59);Free(t78);Free(t80);Free(t87);Free(t93);
}


void dbb7(double* u, double* v, int* n, double* param, double* out)
{
	int i;
	double th, de;
	double *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *t9, *t11, *t12, *t14, *t15, *t16, *t18, *t20, *t23, *t24, *t25, *t27, *t30, *t31, *t32, *t35, *t37, *t42, *t54;
	t1 = Calloc(*n,double);
	t2 = Calloc(*n,double);
	t3 = Calloc(*n,double);
	t4 = Calloc(*n,double);
	t5 = Calloc(*n,double);
	t6 = Calloc(*n,double);
	t7 = Calloc(*n,double);
	t8 = Calloc(*n,double);
	t9 = Calloc(*n,double);
	t11 = Calloc(*n,double);
	t12 = Calloc(*n,double);
	t14 = Calloc(*n,double);
	t15 = Calloc(*n,double);
	t16 = Calloc(*n,double);
	t18 = Calloc(*n,double);
	t20 = Calloc(*n,double);
	t23 = Calloc(*n,double);
	t24 = Calloc(*n,double);
	t25 = Calloc(*n,double);
	t27 = Calloc(*n,double);
	t30 = Calloc(*n,double);
	t31 = Calloc(*n,double);
	t32 = Calloc(*n,double);
	t35 = Calloc(*n,double);
	t37 = Calloc(*n,double);
	t42 = Calloc(*n,double);
	t54 = Calloc(*n,double);

	th = param[0];
	de = param[1];

	for(i=0;i<*n;i++)
	{
		t1[i] = 1.0-u[i];
		t2[i] = pow(t1[i],th);
		t3[i] = 1.0-t2[i];
		t4[i] = pow(t3[i],-de);
		t5[i] = 1.0-v[i];
		t6[i] = pow(t5[i],th);
		t7[i] = 1.0-t6[i];
		t8[i] = pow(t7[i],-de);
		t9[i] = t4[i]+t8[i]-1.0;
		t11[i] = pow(t9[i],-1.0/de);
		t12[i] = 1.0-t11[i];
		t14[i] = pow(t12[i],1.0/th);
		t15[i] = t11[i]*t11[i];
		t16[i] = t14[i]*t15[i];
		t18[i] = 1./t5[i];
		t20[i] = 1./t7[i];
		t23[i] = t9[i]*t9[i];
		t24[i] = 1./t23[i];
		t25[i] = t12[i]*t12[i];
		t27[i] = t24[i]/t25[i];
		t30[i] = t2[i]/t1[i];
		t31[i] = 1./t3[i];
		t32[i] = t30[i]*t31[i];
		t35[i] = t14[i]*t11[i];
		t37[i] = t6[i]*th;
		t42[i] = 1./t12[i];
		t54[i] = t37[i]*t18[i]*t20[i];

		out[i] = -t16[i]*t8[i]*t6[i]*t18[i]*t20[i]*t27[i]*t4[i]*t32[i] + t35[i]*t8[i]*t37[i]*t18[i]*t20[i]*t24[i]*t4[i]*t30[i]*t31[i]*t42[i]+
   t35[i]*t4[i]*t30[i]*t31[i]*t24[i]*t42[i]*t8[i]*de*t54[i]+t16[i]*t4[i]*t32[i]*t27[i]*t8[i]*t54[i];
	}

	Free(t1);Free(t2);Free(t3);Free(t4);Free(t5);Free(t6);Free(t7);Free(t8);Free(t9);
	Free(t11);Free(t12);Free(t14);Free(t15);Free(t16);Free(t18);Free(t20);Free(t23);Free(t24);Free(t25);
	Free(t27);Free(t30);Free(t31);Free(t32);Free(t35);Free(t37);Free(t42);Free(t54);
}


void dbb8(double* u, double* v, int* n, double* param, double* out)
{
	int i;
	double th, de;
	double *t2, *t3, *t4, *t12, *t16, *t6, *t7, *t8, *t10, *t11, *t13, *t15, *t17, *t33, *t38, *t39, *t49, *t59, *t69, *t25, *t26, *t29, *t44, *t45, *t50, *t54, *t62, *t67;
	t2 = Calloc(*n,double);
	t3 = Calloc(*n,double);
	t4 = Calloc(*n,double);
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
	t33 = Calloc(*n,double);
	t38 = Calloc(*n,double);
	t39 = Calloc(*n,double);
	t49 = Calloc(*n,double);
	t59 = Calloc(*n,double);
	t69 = Calloc(*n,double);
	t25 = Calloc(*n,double);
	t26 = Calloc(*n,double);
	t29 = Calloc(*n,double);
	t44 = Calloc(*n,double);
	t45 = Calloc(*n,double);
	t50 = Calloc(*n,double);
	t54 = Calloc(*n,double);
	t62 = Calloc(*n,double);
	t67 = Calloc(*n,double);

	th = param[0];
	de = param[1];

	for(i=0;i<*n;i++)
	{
		  t2[i] = 1.0-de*u[i];
		  t3[i] = pow(t2[i],th);
		  t4[i] = 1.0-t3[i];
		  t10[i] = 1.0-de;
		  t11[i] = pow(t10[i],th);
		  t12[i] = 1.0-t11[i];
		  t13[i] = 1/t12[i];
		  t16[i] = 1/th;
		  t33[i] = th*t3[i];
		  t38[i] = 2.0*th;
		  t39[i] = pow(t10[i],t38[i]);
		  t49[i] = pow(t2[i],t38[i]);
		  t59[i] = pow(t10[i],3.0*th);
		  t69[i] = t12[i]*t12[i];
		  t6[i] = 1.0-de*v[i];
		  t7[i] = pow(t6[i],th);
		  t8[i] = 1.0-t7[i];
		  t15[i] = 1.0-(1.0-t3[i])*t8[i]*t13[i];
		  t17[i] = pow(t15[i],t16[i]);
		  t25[i] = t3[i]*t7[i];
		  t26[i] = t11[i]-t7[i]-t3[i]+t25[i];
		  t29[i] = pow(-t26[i]/t12[i],t16[i]);
		  t44[i] = pow(t6[i],t38[i]);
		  t45[i] = t3[i]*t44[i];
		  t50[i] = t49[i]*t7[i];
		  t54[i] = t49[i]*t44[i];
		  t62[i] = -2.0*t25[i]*t11[i]+t25[i]-t33[i]*t7[i]+3.0*t33[i]*t7[i]*t11[i]-3.0*t33[i]*t7[i]*t39[i]+t25[i]*t39[i]+
			  2.0* t45[i]*t11[i]-t45[i]*t39[i]+2.0*t50[i]*t11[i]-t50[i]*t39[i]-2.0*t54[i]*t11[i]+t54[i]*t39[i]+t54[i]-
			  t50[i]-t45[i]+t33[i]*t7[i]*t59[i];
		  t67[i] = t26[i]*t26[i];
		  out[i] = -de*t29[i]*t62[i]/t6[i]/t2[i]/t67[i]/t69[i];
	}

	Free(t2);Free(t3);Free(t12);Free(t16);
	Free(t6);Free(t7);Free(t8);Free(t10);
	Free(t11);Free(t13);Free(t15);Free(t17);
	Free(t33);Free(t38);Free(t39);Free(t49);
	Free(t59);Free(t69);Free(t25);Free(t26);
	Free(t29);Free(t44);Free(t45);Free(t50);
	Free(t54);Free(t62);Free(t67);Free(t4);
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

  if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30)) )	// 90° rotated copulas
    {
	  nfamily = (*family)-20;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      LL(&nfamily, n, u,  negv, &ntheta, &nnu, loglik);
    }
  else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40)) )	// 270° rotated copulas
    {
	  nfamily = (*family)-30;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      LL(&nfamily, n, negu,  v, &ntheta, &nnu, loglik);
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

  if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30)))	// 90° rotated copulas
    {
	  nfamily = (*family)-20;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      LL(&nfamily, n, negu,  v, &ntheta, &nnu, loglik);
    }
  else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40)))	// 270° rotated copulas
    {
	  nfamily = (*family)-30;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      LL(&nfamily, n, u,  negv, &ntheta, &nnu, loglik);
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
      else ll += log(f);
    }
  }
  else if(*family==3) //Clayton
  {
    if(*theta == 0) ll = 0; 
	  else if(*theta < XEPS) ll = 0;
    else
      {
	  for(j=0;j<*n;j++)
	  {
	    dat[0] = u[j]; dat[1] = v[j];
		f=log(1+*theta)-(1+*theta)*log(dat[0]*dat[1])-(2+1/(*theta))*log(pow(dat[0],-*theta)+pow(dat[1],-*theta)-1.0);
		ll +=f;
	  }
      }
  }
  else if(*family==4) //Gumbel
  {
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
	  f= -pow(t1,1/(*theta))+(2/(*theta)-2)*log(t1)+(*theta-1)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log(1+(*theta-1)*pow(t1,-1.0/(*theta)));

	  ll += f;
    }
  }
  else if(*family==5) // Frank
  {
    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      f = (*theta*(exp(*theta)-1.0)*exp(*theta*dat[1]+*theta*dat[0]+*theta))/pow(exp(*theta*dat[1]+*theta*dat[0])-exp(*theta*dat[1]+*theta)-exp(*theta*dat[0]+*theta)+exp(*theta),2.0);
      if(log(f)>XINFMAX) ll += log(XINFMAX);
      else ll += log(f);
    }
  }
  else if(*family==6)	//Joe
  {
	  for(j=0;j<*n;j++)
		{
		f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
		if(log(f)>XINFMAX) ll += log(XINFMAX);
		else ll += log(f);
		}
  }
  else if(*family==7)	//BB1
  {
		if(*theta == 0){
      for(j=0;j<*n;j++)
      {
         dat[0] = u[j]; dat[1] = v[j];
         t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
  	     f= -pow(t1,1/(*nu))+(2/(*nu)-2)*log(t1)+(*nu-1)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log(1+(*nu-1)*pow(t1,-1.0/(*nu)));
  	     ll += f;
      }		
		}else{
		
			double *param, *fuc;
			param=Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			fuc = Calloc(*n,double);
			dbb1(u, v, n, param, fuc);
			for(j=0;j<*n;j++)
			{
				if(!isfinite(fuc[j]) || isnan(fuc[j]))
				{
					fuc[j]=1;
				}

				if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
				else ll += log(fuc[j]);
			}
			Free(fuc); Free(param);
		}
	}
	else if(*family==8)	//BB6
	{
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		dbb6(u, v, n, param, fuc);
		for(j=0;j<*n;j++)
		{
			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			else ll += log(fuc[j]);
		}
		Free(fuc); Free(param);
	}
	else if(*family==9)	//BB7
	{
		if(*nu==0)
		{
			  for(j=0;j<*n;j++)
				{
				f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
				if(log(f)>XINFMAX) ll += log(XINFMAX);
				else ll += log(f);
				}		
		}
		else
		{
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		dbb7(u, v, n, param, fuc);
		for(j=0;j<*n;j++)
		{	
			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			else ll += log(fuc[j]);
		}
		Free(fuc); Free(param);
		}
  }
  else if(*family==10) //BB8
  {
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		dbb8(u, v, n, param, fuc);
		for(j=0;j<*n;j++)
		{
			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			else ll += log(fuc[j]);
		}
		Free(fuc); Free(param);

  }
  else if(*family==13) //rotated Clayton (180°)
  {
		if(*theta == 0) ll = 0; 
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
  else if(*family==17) //rotated BB1
  {
		if(*theta == 0){
      for(j=0;j<*n;j++)
      {
         dat[0] = 1-u[j]; dat[1] = 1-v[j];
         t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
  	     f= -pow(t1,1/(*nu))+(2/(*nu)-2)*log(t1)+(*nu-1)*log(log(dat[0])*log(dat[1]))-log(dat[0]*dat[1])+log(1+(*nu-1)*pow(t1,-1.0/(*nu)));
  	     ll += f;
      }		
		}else{	  
    
		double *param, *fuc;
			param=Calloc(2,double);
			param[0]=*theta;
			param[1]=*nu;
			fuc = Calloc(*n,double);
			int k=1;
			for(j=0;j<*n;j++)
			{
				dat[0] = 1-u[j]; dat[1] = 1-v[j];

				dbb1(&dat[0], &dat[1], &k, param, &fuc[j]);

				if(!isfinite(fuc[j]) || isnan(fuc[j]))
				{
					fuc[j]=1;
				}

				if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
				else ll += log(fuc[j]);
			}
			Free(fuc); Free(param);
		}
  }
  else if(*family==18) //rotated BB6
  {
	  double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		int k=1;
		for(j=0;j<*n;j++)
			{
				dat[0] = 1-u[j]; dat[1] = 1-v[j];

				dbb6(&dat[0], &dat[1], &k, param, &fuc[j]);

			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			else ll += log(fuc[j]);
		}
		Free(fuc); Free(param);
  }
  else if(*family==19) //rotated BB7
  {
		if(*nu==0){
  	  for(j=0;j<*n;j++)
  		{
  		f = pow(pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta),1/(*theta)-2)*pow(u[j],*theta-1)*pow(v[j],*theta-1)*(*theta-1+pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta));
  		if(log(f)>XINFMAX) ll += log(XINFMAX);
  		else ll += log(f);
  		}		
		}else{
	  double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		int k=1;
		
		for(j=0;j<*n;j++)
		{
			dat[0] = 1-u[j]; dat[1] = 1-v[j];
			dbb7(&dat[0], &dat[1], &k, param, &fuc[j]);

			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			else ll += log(fuc[j]);
		}
		Free(fuc); Free(param);
		}
  }
  else if(*family==20) //rotated BB8
  {
	  double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		int k=1;
		
		for(j=0;j<*n;j++)
		{
			dat[0] = 1-u[j]; dat[1] = 1-v[j];
			dbb7(&dat[0], &dat[1], &k, param, &fuc[j]);

			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			if(log(fuc[j])>XINFMAX) ll += log(XINFMAX);
			else ll += log(fuc[j]);
		}
		Free(fuc); Free(param);
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

  if(((*family==23) | (*family==24) | (*family==26) | (*family==27) | (*family==28) | (*family==29) | (*family==30)) )	// 90° rotated copulas
    {
	  nfamily = (*family)-20;
      for (int i = 0; i < *n; ++i) {negu[i] = 1 - u[i];}
      copLik(&nfamily, n, negu,  v, &ntheta, &nnu, coplik);
    }
  else if(((*family==33) | (*family==34) | (*family==36) | (*family==37) | (*family==38) | (*family==39) | (*family==40)) )	// 270° rotated copulas
    {
	  nfamily = (*family)-30;
      for (int i = 0; i < *n; ++i) {negv[i] = 1 - v[i];}
      copLik(&nfamily, n, u,  negv, &ntheta, &nnu, coplik);
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
    if(*theta == 0) lik = 1.0; 
	if(*theta < XEPS) lik = 1.0;
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
      lik *= f;
    }
  }
  else if(*family==6)	//Joe
  {
	  for(j=0;j<*n;j++)
		{
		f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
		lik *= f;
		}
  }
  else if(*family==7)	//BB1
  {
		if(*theta == 0){

    for(j=0;j<*n;j++)
    {
      dat[0] = u[j]; dat[1] = v[j];
      t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
      t2 = exp(-pow(t1,1.0/(*nu)));
      f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*nu))*pow(log(dat[0])*log(dat[1]),*nu-1.0)*(1.0+(*nu-1.0)*pow(t1,-1.0/(*nu)));
      lik *= f;	
    }
			
		}else{
		
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		dbb1(u, v, n, param, fuc);
		for(j=0;j<*n;j++)
		{
			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			lik *= fuc[j];
		}
		Free(fuc); Free(param);
		}
  }
  else if(*family==8)	//BB6
	{
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		dbb6(u, v, n, param, fuc);
		for(j=0;j<*n;j++)
		{
			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			lik *= fuc[j];
		}
		Free(fuc); Free(param);
	}
  else if(*family==9)	//BB7
  {
		if(*nu == 0){
  	  for(j=0;j<*n;j++)
  		{
  		f = pow(pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta),1/(*theta)-2)*pow(1-u[j],*theta-1)*pow(1-v[j],*theta-1)*(*theta-1+pow(1-u[j],*theta)+pow(1-v[j],*theta)-pow(1-u[j],*theta)*pow(1-v[j],*theta));
  		lik *= f;
  		}		
		}else{
    
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		dbb7(u, v, n, param, fuc);
		for(j=0;j<*n;j++)
		{
			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			lik *= fuc[j];
		}
		Free(fuc); Free(param);
		}
  }
  else if(*family==10)	//BB8
	{
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		dbb8(u, v, n, param, fuc);
		for(j=0;j<*n;j++)
		{
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
		if(*theta == 0) lik = 1.0; 
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
				dat[0] = 1-u[j]; dat[1] = 1-v[j];
				t1 = pow(-log(dat[0]),*theta)+pow(-log(dat[1]),*theta);
				t2 = exp(-pow(t1,1.0/(*theta)));
				f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*theta))*pow(log(dat[0])*log(dat[1]),*theta-1.0)*(1.0+(*theta-1.0)*pow(t1,-1.0/(*theta)));
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
  else if(*family==17)	//rotated BB1
	{
		if(*theta == 0){

    for(j=0;j<*n;j++)
    {		
      dat[0] = 1-u[j]; dat[1] = 1-v[j];
      t1 = pow(-log(dat[0]),*nu)+pow(-log(dat[1]),*nu);
      t2 = exp(-pow(t1,1.0/(*nu)));
      f = t2/(dat[0]*dat[1])*pow(t1,-2.0+2.0/(*nu))*pow(log(dat[0])*log(dat[1]),*nu-1.0)*(1.0+(*nu-1.0)*pow(t1,-1.0/(*nu)));
      lik *= f;		
		}
		
		}else{
		
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		int k=1;
		for(j=0;j<*n;j++)
			{
				dat[0] = 1-u[j]; dat[1] = 1-v[j];

				dbb1(&dat[0], &dat[1], &k, param, &fuc[j]);

			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			lik *= fuc[j];
		}
		Free(fuc); Free(param);
		}
	}
  else if(*family==18)	//rotated BB6
	{
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		int k=1;
		for(j=0;j<*n;j++)
			{
				dat[0] = 1-u[j]; dat[1] = 1-v[j];

				dbb6(&dat[0], &dat[1], &k, param, &fuc[j]);

			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			lik *= fuc[j];
		}
		Free(fuc); Free(param);
	}
  else if(*family==19)	//rotated BB7
	{
		if(*nu == 0){
  	  for(j=0;j<*n;j++)
  		{
  		f = pow(pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta),1/(*theta)-2)*pow(u[j],*theta-1)*pow(v[j],*theta-1)*(*theta-1+pow(u[j],*theta)+pow(v[j],*theta)-pow(u[j],*theta)*pow(v[j],*theta));
  		lik *= f;
  		}		
		}else{		
    
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		int k=1;
		
		for(j=0;j<*n;j++)
		{
			dat[0] = 1-u[j]; dat[1] = 1-v[j];

			dbb7(&dat[0], &dat[1], &k, param, &fuc[j]);

			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			lik *= fuc[j];
		}
		Free(fuc); Free(param);
		}
	}
  else if(*family==20)	//rotated BB8
	{
		double *param, *fuc;
		param=Calloc(2,double);
		param[0]=*theta;
		param[1]=*nu;
		fuc = Calloc(*n,double);
		int k=1;
		
		for(j=0;j<*n;j++)
		{
			dat[0] = 1-u[j]; dat[1] = 1-v[j];

			dbb8(&dat[0], &dat[1], &k, param, &fuc[j]);

			if(!isfinite(fuc[j]) || isnan(fuc[j]))
			{
				fuc[j]=1;
			}

			lik *= fuc[j];
		}
		Free(fuc); Free(param);
	}

  else printError("Error in copLik\t:","Unknown copula family");
  //Free memory:
  Free(dat);
  //Write to output vector:
  *coplik = lik;
}



/////////////////////////////////////////////////////////////
// log likelihood einzeln 
//
/////////////////////////////////////////////////////////////

void LL_mod_seperate(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik)
{
    int kk=1;
    for(int i=0; i<(*n); i++){
        LL_mod(family,&kk,&u[i],&v[i],theta,nu,&loglik[i]);
    };
}
