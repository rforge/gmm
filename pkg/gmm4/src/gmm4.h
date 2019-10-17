#ifndef R_GMM4_H
#define R_GMM4_H

#include <R_ext/RS.h>


void F77_SUB(wu)(double *gt, double *tol, int *maxit,
		  int *n, int *q, int *k, int *conv,
		  double *obj, double *lam);
		 
void F77_SUB(lamcuep)(double *gt, int *n, int *q, int *k,
		       int *maxit, int *conv, double *lam,
		       double *pt, double *obj);

#endif		     