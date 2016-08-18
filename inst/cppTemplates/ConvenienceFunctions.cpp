//#include "ConvenienceFunctions.h"
//#include "UserDefinitionSpace.h"
#include <Rcpp.h>
// check whether each element of an array is between given bounds
int inRange(double x[], double lBound[], double uBound[], int n) {
  int ii = 0;
  while(x[ii] >= lBound[ii] && x[ii] <= uBound[ii] && ii < n) ii++;
  return(ii == n);
}

int isNaN(double x[], int n) {
  int ii = 0;
  //for(ii = 0; ii < n; ii++) Rprintf("x[%i] = %f\n", ii, x[ii]);
  //ii = 0;
  while(x[ii] == x[ii] && ((2.0 * x[ii] != x[ii] && -2.0 * x[ii] != x[ii]) || x[ii] == 0.0) && ii < n) ii++;
  return(ii < n);
}

// cholesky decomposition of a symmetric, positive definite matrix.
// returns a vector of the *Upper Triangular* cholesy factor, leaving the other elements of the array unchanged.
// in other words, U' %*% U \neq A, but lowTri(U') %*% upTri(U) = A.
// both U and A are stacked by COLUMN
void chol_decomp(double U[], double A[], int n) {
  int ii, jj, kk, colI, colJ;
  double tmpSum, tmpInv;
  for(ii = 0; ii < n; ii++) {
    colI = ii*n;
    tmpSum = 0.0;
    //Rprintf("ii = %i\n", ii);
    for(kk = 0; kk < ii; kk++) {
      //Rprintf("ii = %i, kk = %i\n", ii, kk);
      tmpSum += U[colI + kk] * U[colI + kk];
    }
    tmpInv = sqrt(A[colI + ii] - tmpSum);
    U[colI + ii] = tmpInv;
    tmpInv = 1.0/tmpInv;
    for(jj = ii+1; jj < n; jj++) {
      colJ = jj*n;
      tmpSum = 0.0;
      for(kk = 0; kk < ii; kk++) tmpSum += U[colJ + kk] * U[colI + kk];
      U[colJ + ii] = tmpInv * (A[colJ + ii] - tmpSum);
    }
  }
  return;
}

// x = sd * z + mean.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
void xmvn(double x[], double z[], double mean[], double cholSd[], int n) {
  int ii, jj, colI;
  for(ii = 0; ii < n; ii++) {
    colI = n*ii;
    x[ii] = 0;
    for(jj = 0; jj <= ii; jj++) x[ii] += cholSd[colI + jj] * z[jj];
    x[ii] += mean[ii];
  }
  return;
}

// z = sd^{-1} * (x - mean).  only calculates first nMax values of z.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
void zmvn(double z[], double x[], double mean[], double cholSd[], int n, int nMax) {
  int ii, jj, colI;
  double tmpSum;
  for(ii = 0; ii < nMax; ii++) z[ii] = x[ii] - mean[ii];
  // forward substitution
  for(ii = 0; ii < nMax; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) tmpSum += cholSd[colI + jj] * z[jj];
    z[ii] = (z[ii] - tmpSum)/cholSd[colI + ii];
  }
  return;
}

// log-normal density evaluation.  z[] is required as temporary storage of residuals.
// i.e., z = sd^{-1} * (x - mean)
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
double lmvn(double x[], double z[], double mean[], double cholSd[], int n) {
  int ii, jj, colI;
  double tmpSum;
  for(ii = 0; ii < n; ii++) z[ii] = x[ii] - mean[ii];
  // forward substitution
  for(ii = 0; ii < n; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) tmpSum += cholSd[colI + jj] * z[jj];
    z[ii] = (z[ii] - tmpSum)/cholSd[colI + ii];
  }
  tmpSum = 0.0;
  for(ii = 0; ii < n; ii++) tmpSum += z[ii] * z[ii];
  tmpSum *= 0.5;
  for(ii = 0; ii < n; ii++) tmpSum += log(cholSd[n*ii + ii]);
  return(-tmpSum);
}
