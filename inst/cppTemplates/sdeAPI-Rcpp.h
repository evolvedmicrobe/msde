#ifndef sdeAPI_Rcpp_h
#define sdeAPI_Rcpp_h 1

#include <Rcpp.h>
using namespace Rcpp;

List sdeEulerSim(int nDataOut,
		 int N, int reps, int r, double delta, 
		 int MAXBAD, NumericVector initData, NumericVector params);
		 
// mcmc sampling of posterior distribution
List sdeEulerMCMC(int nParamsOut, int nDataOut, 
		  NumericVector initParams, NumericVector initData,
		  NumericVector deltaT, IntegerVector nObsDim, IntegerVector nCompData, 
		  int nSamples, int burn, 
		  IntegerVector dataOutRow, IntegerVector dataOutCol,
		  double updateParams, double updateData, double updateMulti,
		  NumericVector rwJumpSd, int priorType, List priorParams,
		  IntegerVector multiIndex, int nMultiInd, 
		  NumericVector past_Params, NumericVector past_Data,
		  int pastPriorType, List pastPriorParams, int nMaxMultiTries,
		  int nPast_Samples, int nPastCompData,
		  int updateLogMulti, int nLogMultiOut,
		  int updateLogLik, int nLogLikOut,
		  int updateLastMiss, int nLastMissOut);

// drift and diffusion functions, mainly for debugging
NumericVector sdeDrift(NumericVector xIn, NumericVector thetaIn, int nReps);
NumericVector sdeDiff(NumericVector xIn, NumericVector thetaIn, int nReps);

// SDE log-likelihood evaluation.
NumericVector sdeLogLik(NumericVector xIn, NumericVector thetaIn, NumericVector deltaT, 
			int nCompData, int nReps);
			



#endif
