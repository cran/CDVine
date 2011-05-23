#ifndef VINE_H
#define VINE_H

/*
** vine.c - C code of the package CDRVine    
** 
** with contributions from Carlos Almeida, Aleksey Min, 
** Ulf Schepsmeier, Jakob Stoeber and Eike Brechmann
** 
** A first version was based on code
** from Daniel Berg <daniel at danielberg.no>
** provided by personal communication. 
**
*/

#ifndef   	VINE_H_
# define   	VINE_H_



#endif 	    /* !VINE_H_ */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <setjmp.h>
#include <R.h>
#include <Rdefines.h> 
#include <Rmath.h>
#include <Rinternals.h>
#include "vine.h"


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define pi 3.14159265358979
//#define  XINFMAX 1e308
# define XINFMAX DBL_MAX

/* define boolean type for C */
typedef unsigned int boolean;
#define false 0
#define true (!false)

double **create_matrix(int rows, int columns);
void free_matrix(double **a, int rows);
int **create_intmatrix(int rows, int columns);
void free_intmatrix(int **a, int rows);
double ***create_3darray(int d1, int d2, int d3);
void free_3darray(double ***a, int d1, int d2);
double StableGammaDivision(double x1, double x2);


void Hfunc(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void LL(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik);
void Hderiv(int* family, double* u, double* v, double* theta, double* nu, double* out);
void HNumInv(int* family, double* u, double* v, double* theta, double* nu, double* out);
void HNumInvBisect(int* family, double* u, double* v, double* theta, double* nu, double* out);
void Hinv(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void pcc(int* n, int* d, int* family, int* type, double* par, double* nu, double* out);
void VineLogLikm(int* T, int* d, int* type, int* family, double* par, double* data, double* out, double* ll, double* vv);

void LL_mod(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik);
void LL_mod2(int* family, int* n, double* u, double* v, double* theta, double* nu, double* loglik);
void Hfunc1(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out);
void Hfunc2(int* family,int* n,double* v,double* u,double* theta,double* nu,double* out);
void Hinv1(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void Hinv2(int* family, int* n, double* v, double* u, double* theta, double* nu, double* out);

void copLik(int* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik);
void copLik_mod(int* family, int* n, double* u, double* v, double* theta, double* nu, double* coplik);
void SimulateRVine(int* T, int* d, int* family, int* maxmat, int* matrix, int* conindirect, double* par, double* par2, double* out);

void ktau(double *X, double *Y, int *N, double *tau, double *S, double *D, int *T, int *U, int *V);

#endif
