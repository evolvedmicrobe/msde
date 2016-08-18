#ifndef sdeCore_h
#define sdeCore_h 1

#include <Rcpp.h>
using namespace Rcpp;

extern const int nParams; //Number of model paramaters
extern const int nDims; //Number of model dimensions
extern const int nMaxPriorParams;
const int nMaxPriorParams = nParams + nDims;

class PropSim {
 private:
  double *mean;
  double *sd;
  double *z;
  int ncomp;
 public:
  PropSim(int pnComp) {
    ncomp = pnComp;
    mean = new double[ncomp*nDims];
    sd = new double[ncomp*nDims*nDims];
    z = new double[ncomp*nDims];
    for(int i = 0; i < ncomp*nDims; i++) {
      mean[i] = 0.0;
      z[i] = 0.0;
    }
    for(int i = 0; i < ncomp*nDims*nDims; i++) {
      sd[i] = 0.0;
    }
  }
  double *Mean() {
    return mean;
  }
  double *Sd() {
    return sd;
  }
  double *Z() {
    return z;
  }
  int nComp() {
    return ncomp;
  }
};


// eraker proposal mean and standard deviatiation
// NOTE: sde = upper triangular cholesky factor
void mvEraker(double mean[], double sd[], double x0[], double x2[], double b, double b2, double params[]);

// euler approximation mean and standard deviation
// NOTE: sde = upper triangular cholesky factor
void mvEuler(double mean[], double sd[], double x0[], double t, double t2, double params[]);

// componentwise data updates.
int missGibbsUpdate(int missInd[], int nMiss, double jumpSd[], Prior *prior, int gibbsAccept[],
		    int paramAccept[], PropSim *propSim, int *parIndex, int nMiss0, double *propAccept,
		    double *currData, double *currParams, double *propData,
		    double *dT, double *sqrtDT, double *B, double *sqrtB);

// loglikelihood evaluation.  this is always done on currData.
double loglik(double params[], PropSim *propSim, double *currData, double *dT, double *sqrtDT);

// componentwise vanilla MH parameter updates
void paramVnlUpdate(double jumpSD[], Prior *prior, int paramAccept[], PropSim *propSim, double *currData,
					double *currParams, double *propParams, double *dT, double *sqrtDT);

// multiresolution updates using eraker proposals
// flow is as follows: (1) drop down all low res samples
// (2) for each new missing data point, drop down + calculate acceptance rate
// this flow allows for easy parallelization of the relevant for-loops.
void multiUpdate(int multiInd[], int pastMultiInd[], int nMulti, int multiAccept[],
		 double logMultiAcc[], int nPastComp, int nPastSamples, int nMaxMultiTries,
		 Prior *prior, Prior *pastPrior, PropSim *propSim, double *pastData,
		 double *pastParams, double *currData, double *currParams, double *propData, double *propParams,
		 double *dT, double *sqrtDT, double *B, double *sqrtB, double *pastDT, double *pastSqrtDT);

// user functions
void sdeDr(double dr[], double x[], double t, double params[]);
void sdeDf(double df[], double x[], double sqrtT, double params[]);
int isValidData(double x[]);
int isValidParams(double params[]);

#endif
