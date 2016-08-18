//#include "sdeAPI-Rcpp.h"
//#include "ConvenienceFunctions.h"
//#include "UserDefinitionSpace.h"
//#include "sdeCore.h"

#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export("sde.model$sim")]]
List sdeEulerSim(int nDataOut,
		 int N, int reps, int r, double delta, 
		 int MAXBAD, NumericVector initData, NumericVector params) {
  RNGScope scope;
		
  double sqrtDelta = sqrt(delta);
  int ii, jj, kk;
  int bad = 0;
  // output
  NumericVector dataOut(nDataOut);
  int nBadDraws;
		
  double *X = new double[nDims];
  double *tmpX = new double[nDims];
  PropSim *propSim = new PropSim(1);
  // initialize
  for(ii = 0; ii < nDims; ii++) {
    X[ii] = 0.0;
    tmpX[ii] = 0.0;
  }
		
  for(ii = 0; ii < reps; ii++) {
    for(jj = 0; jj < N*r; jj++) {
      // initialize chains
      if(jj == 0) {
	for(kk = 0; kk < nDims; kk++) {
	  X[kk] = initData[ii*nDims + kk];
	}
      }
      else {
	mvEuler(propSim->Mean(), propSim->Sd(), X, delta, sqrtDelta, &params[ii*nParams]);
	for(kk = 0; kk < nDims; kk++) {
	  propSim->Z()[kk] = norm_rand();
	}
	xmvn(tmpX, propSim->Z(), propSim->Mean(), propSim->Sd(), nDims);
	// validate draw
	while(!isValidData(tmpX) && bad < MAXBAD) {
	  for(kk = 0; kk < nDims; kk++) {
	    propSim->Z()[kk] = norm_rand();
	  }
	  xmvn(tmpX, propSim->Z(), propSim->Mean(), propSim->Sd(), nDims);
	  bad++;
	}
	if (bad == MAXBAD) {
	  goto stop;
	}
	for(kk = 0; kk < nDims; kk++) {
	  X[kk] = tmpX[kk];
	}
      }
      // store
      if(jj % r == 0) {
	for(kk = 0; kk < nDims; kk++) {
	  dataOut[ii*N*nDims + (jj/r)*nDims + kk] = X[kk];
	}
      }
    }
  }
		
 stop:
  nBadDraws = bad;

  delete [] X;
  delete [] tmpX;
		
  return List::create(_["dataOut"] = dataOut, _["nBadDraws"] = nBadDraws);
}

// mcmc sampling of posterior distribution
//[[Rcpp::export("sde.model$post")]]
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
		  int updateLastMiss, int nLastMissOut) {
  RNGScope scope;
  int ii, jj, kk, updt, foundNaN;

  // input
  double *jumpSd = REAL(rwJumpSd);
  int *multiInd = INTEGER(multiIndex);
  int nMulti = nMultiInd;

  // initialize external variables
  double *dT = REAL(deltaT);
  int nComp = INTEGER(nCompData)[0];
  int nCompOut = INTEGER(nCompData)[1];
  int *parIndex = INTEGER(nObsDim);
  // unobserved states in first observation
  int nMiss0 = nDims-parIndex[0];
  // unobserved states in last observation
  int nMissN = nDims-parIndex[nComp-1];
  double *pastParams = REAL(past_Params);
  double *pastData = REAL(past_Data);
  int nPastSamples = nPast_Samples;
  int nPastComp = nPastCompData;

  // output variables
  NumericVector paramsOut(nParamsOut);
  NumericVector dataOut(nDataOut);
  NumericVector logMultiOut(nLogMultiOut);
  IntegerVector paramAcceptOut(nParams + nMiss0);
  IntegerVector gibbsAcceptOut(nComp);
  IntegerVector multiAcceptOut(1);
  NumericVector logLikOut(nLogLikOut);
  NumericVector lastMissOut(nLastMissOut);
  NumericVector lastIter(nParams + nComp*nDims);
  int *paramAccept = INTEGER(paramAcceptOut);
  int *gibbsAccept = INTEGER(gibbsAcceptOut);
  int *multiAccept = INTEGER(multiAcceptOut);

  // times
  double *sqrtDT = new double[nComp];
  double *B = new double[nComp];
  double *sqrtB = new double[nComp];
  for(ii = 0; ii < nComp-1; ii++) {
    sqrtDT[ii] = sqrt(dT[ii]);
    if(ii > 0) {
      B[ii] = dT[ii]/(dT[ii-1] + dT[ii]);
      sqrtB[ii] = sqrt((1-B[ii]) * dT[ii]);
    }
  }

  // data
  double *currData = new double[nComp*nDims];
  double *propData = new double[nComp*nDims];
  double *propAccept = new double[nComp];
  PropSim *propSim = new PropSim(nComp);
  for(ii = 0; ii < nComp; ii++) {
    propAccept[ii] = 0.0;
    for(jj = 0; jj < nDims; jj++) {
      currData[ii*nDims + jj] = initData[ii*nDims + jj];
    }
  }
  for(ii = 0; ii < nComp*nDims; ii++) {
    propData[ii] = 0.0;
  }
	
  // identify missing data indices, i.e. at least one component to update
  int nMiss = 0;
  for(ii = 0; ii < nComp; ii++) {
    nMiss += (parIndex[ii] < nDims);
  }
  int *missInd = new int[nMiss + (nMiss == 0)];
  jj = 0;
  for(ii = 0; ii < nComp; ii++) {
    if(parIndex[ii] < nDims) {
      missInd[jj++] = ii;
    }
  }
	
  // params
  double *currParams = new double[nParams];
  double *propParams = new double[nParams];
  for(ii = 0; ii < nParams; ii++) {
    currParams[ii] = initParams[ii];
    propParams[ii] = currParams[ii];
  }

  // multiresolution
  double *logMultiAcc = new double[nMulti];
  double *pastDT = new double[nMulti];
  double *pastSqrtDT = new double[nMulti];
  int *pastMultiInd = &multiInd[nMulti];
  for(ii = 0; ii < nMulti; ii++) {
    pastDT[ii] = dT[multiInd[ii]-1] + dT[multiInd[ii]];
    pastSqrtDT[ii] = sqrt(pastDT[ii]);
  }

  // parameter prior
  Prior *prior = NULL;
  Prior *pastPrior = NULL;
  if((Prior::Type)priorType == Prior::Flat) {
    prior = new FlatPrior();
  }
  else if((Prior::Type)priorType == Prior::Normal) {
    prior = new NormalPrior(priorParams, parIndex[0], nMiss0);
  }
  else if((Prior::Type)priorType == Prior::GCop) {
    prior = new GCopPrior(priorParams, parIndex[0], nMiss0);
  }
  else if((Prior::Type)priorType == Prior::Custom) {
    prior = new CustomPrior(priorParams, parIndex[0], nMiss0);
  }
  else {
    throw Rcpp::exception("ERROR: Unrecognized prior type\n");
  }

  // past parameter prior
  if(updateMulti > 0.0) {
    if((Prior::Type)pastPriorType == Prior::Flat) {
      pastPrior = new FlatPrior();
    }
    else if((Prior::Type)pastPriorType == Prior::Normal) {
      pastPrior = new NormalPrior(pastPriorParams, parIndex[0], nMiss0);
    }
    else if((Prior::Type)pastPriorType == Prior::GCop) {
      pastPrior = new GCopPrior(pastPriorParams, parIndex[0], nMiss0);
    }
    else if((Prior::Type)pastPriorType == Prior::Custom) {
      pastPrior = new CustomPrior(pastPriorParams, parIndex[0], nMiss0);
    }
    else {
      throw Rcpp::exception("ERROR: Unrecognized past prior type\n");
    }
  }
	
  // main MCMC loop
  jj = 0;
  for(int smp = -burn; smp < nSamples; smp++) {
    // missing data update
    updt = 0;
    if(updateData > 0.0) {
      if(updateData < 1.0) {
	updt = unif_rand() <= updateData;
      } 
      else {
	updt = (smp % (int) updateData) == 0;
      }
    }
		
    if(updt) {
      foundNaN = missGibbsUpdate(missInd, nMiss, jumpSd, prior, gibbsAccept, paramAccept, propSim, 
				 parIndex, nMiss0,
				 propAccept, currData, currParams, propData, dT, sqrtDT, B, sqrtB);
    }

    // parameter update
    updt = 0;
    if(updateParams > 0.0) {
      if(updateParams < 1.0) {
	updt = unif_rand() <= updateParams;
      }
      else {
	updt = (smp % (int) updateParams) == 0;
      }
    }
		
    if(updt) {	
      paramVnlUpdate(jumpSd, prior, paramAccept, propSim, currData, currParams, propParams, dT, sqrtDT);
    }
		
    // multiresolution update
    updt = 0;
    if(updateMulti > 0.0) {
      if(updateMulti < 1.0) {
	updt = unif_rand() <= updateMulti;
      } 
      else {
	updt = (smp % (int) updateMulti) == 0;
      }
    }
		
    if(updt) {
      multiUpdate(multiInd, pastMultiInd, nMulti, multiAccept, logMultiAcc, nPastComp, nPastSamples, 
		  nMaxMultiTries, prior, pastPrior, propSim,
		  pastData, pastParams, currData, currParams, propData, propParams, dT, sqrtDT, B, sqrtB, 
		  pastDT, pastSqrtDT);
    }

    // sufficient statistics
    if(smp >= 0) {
      // log-likelihood
      if(updateLogLik) logLikOut[smp] = loglik(currParams, propSim, currData, dT, sqrtDT);
    }
		
    // storage
    if(smp == dataOutRow[jj]) {
      if(updateData > 0.0) {
	for(ii = 0; ii < nCompOut; ii++) {
	  for(kk = 0; kk < nDims; kk++) {
	    dataOut[jj*nDims*nCompOut + ii*nDims + kk] = currData[dataOutCol[ii]*nDims + kk];
	  }
	}
      }
      if((updateMulti > 0.0) & updateLogMulti) {
	for(ii = 0; ii < nMulti; ii++) {
	  logMultiOut[jj*nMulti + ii] = logMultiAcc[ii];
	}
      }
      jj++;
    }
    if(updateParams > 0.0 && smp >= 0) {
      for(ii = 0; ii < nParams; ii++) paramsOut[smp*nParams + ii] = currParams[ii];
    }
    if(updateLastMiss && smp >= 0) {
      for(ii = 0; ii < nMissN; ii++) {
	lastMissOut[smp*nMissN + ii] = currData[(nComp-1)*nDims + parIndex[nComp-1] + ii];
      }
    }

    // NaN breakpoint
    if(foundNaN) {
      Rprintf("smp = %i\n", smp);
      goto stop;
    }
  }
	
	
  // store last iteration, to resume MCMC later if needed
  for(ii = 0; ii < nParams; ii++) {
    lastIter[ii] = currParams[ii];
  }
  for(ii = 0; ii < nDims*nComp; ii++) {
    lastIter[nParams+ii] = currData[ii];
  }
	
 stop:
  // delete dynamic variables
  delete [] sqrtDT;
  delete [] B;
  delete [] sqrtB;
  delete [] currData;
  delete [] propData;
  delete [] propAccept;
  delete [] currParams;
  delete [] propParams;
  delete [] missInd;
  delete [] pastDT;
  delete [] pastSqrtDT;
  delete [] logMultiAcc;
  //delete prior;
  //delete pastPrior;
	
  return List::create(_["paramsOut"] = paramsOut, _["dataOut"] = dataOut,
		      _["logMultiOut"] = logMultiOut, 
		      _["paramAccept"] = paramAcceptOut, 
		      _["gibbsAccept"] = gibbsAcceptOut,
		      _["multiAccept"] = multiAcceptOut, _["logLikOut"] = logLikOut,
		      _["lastMissOut"] = lastMissOut, _["lastIter"] = lastIter);
}

// drift and diffusion functions, mainly for debugging
//[[Rcpp::export("sde.model$drift")]]
NumericVector sdeDrift(NumericVector xIn, NumericVector thetaIn, int nReps) {
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector drOut(nReps*nDims);
  double *dr = REAL(drOut);
  for(int ii = 0; ii < nReps; ii++) {
    sdeDr(&dr[ii*nDims], &x[ii*nDims], 1.0, &theta[ii*nParams]);
  }
  return drOut;
}

//[[Rcpp::export("sde.model$diff")]]
NumericVector sdeDiff(NumericVector xIn, NumericVector thetaIn, int nReps) {
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector dfOut(nReps*nDims*nDims);
  double *df = REAL(dfOut);
  for(int ii = 0; ii < nReps; ii++) {
    sdeDf(&df[ii*nDims*nDims], &x[ii*nDims], 1.0, &theta[ii*nParams]);
  }
  return dfOut;
}

// SDE log-likelihood evaluation.
//[[Rcpp::export("sde.model$loglik")]]
NumericVector sdeLogLik(NumericVector xIn, NumericVector thetaIn, NumericVector deltaT, 
			int nCompData, int nReps) {
  int ii;
  double *currData;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector llOut(nReps);
  double *ll = REAL(llOut);
  int nComp = nCompData;
  double *dT = REAL(deltaT);
  double *sqrtDT = new double[nComp];
  PropSim *propSim = new PropSim(nComp);
  for(ii = 0; ii < nComp-1; ii++) {
    sqrtDT[ii] = sqrt(dT[ii]);
  }
  for(ii = 0; ii < nReps; ii++) {
    currData = &x[ii*nDims*nComp];
    ll[ii] = loglik(&theta[ii*nParams], propSim, currData, dT, sqrtDT);
  }
  delete [] sqrtDT;
	
  return llOut;
}
