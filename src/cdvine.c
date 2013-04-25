/*
** cdvine.c - C code of the package CDRVine  
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
  int i=0, j=0, in=1, k=0, **fam;
  double *w, **v, t=0.0, **theta, **x, **ny;

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
      x[j][1] = w[1];
      for(i=2;i<=*d;i++)
      {
        t = w[i];
        for(k=i-1;k>=1;k--)
        {
		  Hinv1(&fam[k][i-k],&in, &t,&w[k],&theta[k][i-k],&ny[k][i-k],&t);
        }
        x[j][i] = t;
      }
    }
  }
  else if(*type==2) //D-vine
  {
    for(j=1;j<=*n;j++)
    {
      for(i=1;i<=*d;i++) { w[i] = runif(0,1);}
      v[1][1] = w[1];
      v[2][1] = w[2];
      Hinv1(&fam[1][1],&in,&w[2],&v[1][1],&theta[1][1],&ny[1][1],&v[2][1]);
      Hfunc2(&fam[1][1],&in, &v[1][1],&v[2][1],&theta[1][1],&ny[1][1],&v[2][2]);
      for(i=3;i<=*d;i++)
      {
        v[i][1] = w[i];
	
        for(k=i-1;k>=2;k--) { 
	  Hinv1(&fam[k][i-k],&in, &v[i][1],&v[i-1][2*k-2],&theta[k][i-k],&ny[k][i-k],&v[i][1]);
	}
        Hinv1(&fam[1][i-1],&in, &v[i][1],&v[i-1][1],&theta[1][i-1],&ny[1][i-1],&v[i][1]);
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
  int i=0, j=0, k=0, t=0, kk=0, **fam;
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

//////////////////////////////////////////////////////////////
// Function to compute -log-likelihood for the pair-copula construction (vine) 
// Input:
// n        sample size
// d        dimension (>=2)
// type     vine type (1=canonical vine, 2=d-vine)
// family   copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank)
// par      parameter values (at least d*(d-1)/2 parameters
// mpar     index of modified parameter (related to previous computation)
// data     data set for which to compute log-likelihood
// ll       array with the stored contribution of the likelihood in a previous computation
// vv       3d array  array with the stored transformations in a previous computation
// Output:
// ll       array with the contribution to LL (for each copula)
// vv       array for the transformation operated (Hfunc)  
/////////////////////////////////////////////////////////////
void VineLogLikmP(int* T, int* d, int* type, int* family, double* par, int* mpar, double* data, 
		  double* out, double* ll, double* vv)
{
  int i=0, j=0, ii=0, jj=0, k=0, t=0, kk=0 ,**fam, **ind;
  double sumloglik=0.0,loglik=0.0,  **x, **theta, **nu, ***v;
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
  fam = create_intmatrix(*d,*d);
  ind = create_intmatrix(*d,*d);
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
  jj = *d;
  ii = *d;
  kk=00;

  // By Ulf Schepsmeier
	if(*type==1) //C-Vine
	{
		for(i=1;i<=(*d-1);i++)
		{
		  for(j=1;j<=(*d-1);j++)
			{
				ind[i][j] = 0;
				if(j <= *d-i)
				{
					++k;
					if (k == *mpar)
					{
						ii=i;
						jj=j;
					}
					if(ii+jj-i>0)
					{
						if(i>=ii && j==ii+jj-i)
						{
							ind[i][j]=1;
						}
					}
					else if(ii+jj-i<=0)
					{
						ind[i][j]=1;
					}
				}
			}
		}

		for(k=1;k<*d-1;k++)
			for(i=1;i<=*d-k;i++)
				for(t=0;t<*T;t++)
				{
					v[k][i][t] = vv[kk];
  					++kk;
  				}
	}
	else //D-Vine
	{
	  for(i=1;i<=(*d-1);i++)
		{
		  for(j=1;j<=(*d-1);j++)
			{
			ind[i][j] = 0;
			if(j <= *d-i)
			{
				++k;
				if (k == *mpar)
				{
					ii=i;
					jj=j;
				}            
				if(i >= ii &&  j >= ii+jj-i &&  j <= *d-i && j <= jj)
				{
					ind[i][j]=1;
				}
			//	    //printf("%d ", ind[i][j]);
			}
			}
		  //      //printf("\n");
		}

		for(k=1;k<*d;k++)
			for(i=1;i<=2*(*d-k-1);i++)
				for(t=0;t<*T;t++)
				{
					v[k][i][t] = vv[kk];
					++kk;
				}
	}
  k=0;
  for(i=1;i<=(*d-1);i++)
    {
      for(j=1;j<=(*d-i);j++)
	{
	  theta[i][j] = par[k];
	  fam[i][j] = family[k];
	  nu[i][j] = par[*d*(*d-1)/2+k];
	  k ++;
      }
    }
 
  

  if(*type==1) //C-vine
    {
      // By Ulf Schepsmeier
	  kk=0;
	  //Compute likelihood at level 1:
		for(i=1;i<*d;i++)
		{
			if(ind[1][i]==1) 
			{
				LL_mod2(&fam[1][i],T,x[1],x[i+1],&theta[1][i],&nu[1][i],&loglik);
				ll[kk] = loglik;
				//Compute variables for next level:
				Hfunc1(&fam[1][i],T,x[i+1],x[1],&theta[1][i],&nu[1][i],v[1][i]);
			}
			sumloglik += ll[kk];
			++kk;
		}
		//Compute likelihood at next levels:
		for(k=2;k<=(*d-1);k++)
		{
			for(i=1;i<=(*d-k);i++)
			{
				if(ind[k][i]==1)
				{
					LL_mod2(&fam[k][i],T,v[k-1][1],v[k-1][i+1],&theta[k][i],&nu[k][i],&loglik);
					ll[kk] = loglik; 
					if(k<(*d-1))
					{
						Hfunc1(&fam[k][i],T,v[k-1][i+1],v[k-1][1],&theta[k][i],&nu[k][i],v[k][i]);
					}
				}
				sumloglik += ll[kk];
				++kk;
			}
		}
    }
  else if(*type==2) //D-vine
  {
    kk=0;
    //Compute the likelihood at level 1:
    for(i=1;i<*d;i++)
      {
        if(ind[1][i]==1) 
		{
			LL_mod2(&fam[1][i],T,x[i],x[i+1],&theta[1][i],&nu[1][i],&loglik);
			ll[kk] = loglik;
		}
		sumloglik += ll[kk];
		++kk;
      }
    //Compute variables for next level:
    if(ind[1][1]==1) Hfunc2(&fam[1][1],T,x[1],x[2],&theta[1][1],&nu[1][1],v[1][1]);
    for(k=1;k<=(*d-3);k++)
      if(ind[1][k+1]==1)
		{
			Hfunc1(&fam[1][k+1],T,x[k+2],x[k+1],&theta[1][k+1],&nu[1][k+1],v[1][2*k]);
			Hfunc2(&fam[1][k+1],T,x[k+1],x[k+2],&theta[1][k+1],&nu[1][k+1],v[1][2*k+1]);
		}
    if(ind[1][*d-1]) Hfunc1(&fam[1][*d-1],T,x[*d],x[*d-1],&theta[1][*d-1],&nu[1][*d-1],v[1][2*(*d)-4]);
    //Compute likelihood at next levels:
    for(k=2;k<=(*d-1);k++)
      {
        for(i=1;i<=(*d-k);i++)
    	  {
	    if(ind[k][i] == 1)
	      {
		LL_mod2(&fam[k][i],T,v[k-1][2*i-1],v[k-1][2*i],&theta[k][i],&nu[k][i],&loglik);
		ll[kk] = loglik;
	      }
	    sumloglik += ll[kk];
	    ++kk;
    	  }
        if(k<(*d-1))
    	  {
    	    if(ind[k][1] == 1) Hfunc2(&fam[k][1],T,v[k-1][1],v[k-1][2],&theta[k][1],&nu[k][1],v[k][1]);
    	    if((*d)>4)
    	      {
    		for(i=1;i<=(*d-k-2);i++)
    		  {
		    if(ind[k][i+1]==1)
		      {
			Hfunc1(&fam[k][i+1],T,v[k-1][2*i+2],v[k-1][2*i+1],&theta[k][i+1],&nu[k][i+1],v[k][2*i]);
			Hfunc2(&fam[k][i+1],T,v[k-1][2*i+1],v[k-1][2*i+2],&theta[k][i+1],&nu[k][i+1],v[k][2*i+1]);
		      }
		  }
    	      }
    	    if(ind[k][*d-k]==1) Hfunc1(&fam[k][*d-k],T,v[k-1][2*(*d)-2*k],v[k-1][2*(*d)-2*k-1],&theta[k][*d-k],&nu[k][*d-k],v[k][2*(*d)-2*k-2]);
    	  }
      }
  }
  //Write to output:
  *out = -sumloglik;
  kk=00;
  // By Ulf Schepsmeier
	if(*type==1) //C-Vine
	{
		for(k=1;k<*d-1;k++)
			for(i=1;i<=*d-k;i++)
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
  free_matrix(x,*d+1); free_matrix(theta,*d); free_matrix(nu,*d); free_intmatrix(fam,*d);free_intmatrix(ind,*d);
}


void SimulateVine(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* out)
{
	int i=0, j=0, k=0, m=0, **fam, **cindirect, **mat, **mmat, **fam2, **cindirect2, **mat2, **mmat2;
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
	
	// Matrizen rotieren f�r den Algo
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
