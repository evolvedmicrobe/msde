//#include "sdeCore.h"
//#include "ConvenienceFunctions.h"
//#include "UserDefinitionSpace.h"


// eraker proposal mean and standard deviatiation
// NOTE: sde = upper triangular cholesky factor
void mvEraker(double mean[], double sd[], double x0[], double x2[], double b, double b2, double params[]) {
  for(int jj = 0; jj < nDims; jj++) {
    mean[jj] = x0[jj] * b + x2[jj] * (1-b);
  }
  sdeDf(sd, x0, b2, params);
  return;
}

// euler approximation mean and standard deviation
// NOTE: sde = upper triangular cholesky factor
void mvEuler(double mean[], double sd[], double x0[], double t, double t2, double params[]) {
  sdeDr(mean, x0, t, params);
  for(int jj = 0; jj < nDims; jj++) {
    mean[jj] += x0[jj];
  }
  sdeDf(sd, x0, t2, params);
  return;
}

// componentwise data updates.
int missGibbsUpdate(int missInd[], int nMiss, double jumpSd[], Prior *prior, int gibbsAccept[],
		    int paramAccept[], PropSim *propSim, int *parIndex, int nMiss0, double *propAccept,
		    double *currData, double *currParams, double *propData, double *dT, double *sqrtDT,
		    double *B, double *sqrtB) {
  int ii, II, jj, JJ;
  int startII, endII;
  int foundNaN = 0;
  // only elements in missInd are updated.
  // first and last components are handled separately.
  startII = (missInd[0] == 0) ? 1 : 0;
  endII = (missInd[nMiss-1] == propSim->nComp()-1) ? nMiss-1 : nMiss;
  // Markov chain elements are conditionally independent,
  // so every other can be updated in parallel
  for(JJ = 0; JJ < 2; JJ++) {
    for(II = startII+JJ; II < endII; II = II+2) { // *** PARALLELIZABLE FOR-LOOP ***
      ii = missInd[II];
      // intermediate data points
      if((0 < ii) && (ii < propSim->nComp()-1)) {
	mvEraker(&(propSim->Mean()[ii*nDims]), &(propSim->Sd()[ii*nDims*nDims]), &currData[(ii-1)*nDims],
		 &currData[(ii+1)*nDims], B[ii], sqrtB[ii], currParams);
	// partial observations
	if(parIndex[ii] == 0) {
	  for(jj = 0; jj < nDims; jj++) propSim->Z()[ii*nDims + jj] = norm_rand();
	}
	else {
	  zmvn(&(propSim->Z()[ii*nDims]), &currData[ii*nDims], &(propSim->Mean()[ii*nDims]),
	       &(propSim->Sd()[ii*nDims*nDims]), nDims, parIndex[ii]);
	  for(jj = parIndex[ii]; jj < nDims; jj++) {
	    propSim->Z()[ii*nDims + jj] = norm_rand();
	  }
	}
	// proposals
	xmvn(&propData[ii*nDims], &(propSim->Z()[ii*nDims]), &(propSim->Mean()[ii*nDims]),
	     &(propSim->Sd()[ii*nDims*nDims]), nDims);
	if(isNaN(&propData[ii*nDims], nDims)) {
	  Rprintf("found NaN in missGibbsUpdate, ii = %i.\n", ii);
	  foundNaN = 1;
	  goto stop;
	}
	// only calculate acceptance rate if proposal is valid
	if(isValidData(&propData[ii*nDims])) {
	  // acceptance rate
	  // proposal
	  propAccept[ii] = lmvn(&currData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
				&propSim->Sd()[ii*nDims*nDims], nDims);
	  propAccept[ii] -= lmvn(&propData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
				 &propSim->Sd()[ii*nDims*nDims], nDims);


	  // target 1
	  mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &currData[(ii-1)*nDims],
		  dT[ii-1], sqrtDT[ii-1], currParams);
	  propAccept[ii] += lmvn(&propData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
				 &propSim->Sd()[ii*nDims*nDims], nDims);
	  propAccept[ii] -= lmvn(&currData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
				 &propSim->Sd()[ii*nDims*nDims], nDims);


	  // target 2
	  mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &propData[ii*nDims], dT[ii],
		  sqrtDT[ii], currParams);
	  propAccept[ii] += lmvn(&currData[(ii+1)*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
				 &propSim->Sd()[ii*nDims*nDims], nDims);
	  mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &currData[ii*nDims], dT[ii],
		  sqrtDT[ii], currParams);
	  propAccept[ii] -= lmvn(&currData[(ii+1)*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
				 &propSim->Sd()[ii*nDims*nDims], nDims);


	  // evaluate mh ratio
	  if(exp(propAccept[ii]) >= unif_rand()) {
	    for(jj = 0; jj < nDims; jj++) {
	      currData[ii*nDims + jj] = propData[ii*nDims + jj];
	    }
	    gibbsAccept[ii]++;
	  }
	}
      }
    }
  }
  // special treatment for last datapoint
  if(missInd[nMiss-1] == propSim->nComp()-1) {
    ii = propSim->nComp()-1;
    mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &currData[(ii-1)*nDims],
	    dT[ii-1], sqrtDT[ii-1], currParams);
    // partial observations
    if(parIndex[ii] == 0) {
      for(jj = 0; jj < nDims; jj++) {
	propSim->Z()[ii*nDims + jj] = norm_rand();
      }
    }
    else {
      zmvn(&propSim->Z()[ii*nDims], &currData[ii*nDims], &propSim->Mean()[ii*nDims],
	   &propSim->Sd()[ii*nDims*nDims], nDims, parIndex[ii]);
      for(jj = parIndex[ii]; jj < nDims; jj++) {
	propSim->Z()[ii*nDims + jj] = norm_rand();
      }
    }
    // proposals
    xmvn(&propData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
	 &propSim->Sd()[ii*nDims*nDims], nDims);
    if(isNaN(&propData[ii*nDims], nDims)) {
      Rprintf("found NaN in missGibbsUpdate, ii = %i.\n", ii);
      foundNaN = 1;
      goto stop;
    }
    // acceptance is 100% as long as the proposal is valid
    if(isValidData(&propData[ii*nDims])) {
      for(jj = 0; jj < nDims; jj++) {
	currData[ii*nDims + jj] = propData[ii*nDims + jj];
      }
      gibbsAccept[ii]++;
    }
  }
  // special treatment for first datapoint
  if(missInd[0] == 0) {
    ii = 0;
    // initialize
    for(jj = 0; jj < nDims; jj++) {
      propData[jj] = currData[jj];
    }


    // random walk metropolis
    for(jj = 0; jj < nMiss0; jj++) {
      // proposal
      propData[parIndex[0]+jj] = currData[parIndex[0]+jj] + jumpSd[nParams+jj] * norm_rand();
      if(isValidData(&propData[ii*nDims])) {
	// acceptance rate.
	// target 1
	propAccept[ii] = prior->logPrior(currParams, propData);
	propAccept[ii] -= prior->logPrior(currParams, currData);


	// target 2
	mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &propData[ii*nDims], dT[ii],
		sqrtDT[ii], currParams);
	propAccept[ii] += lmvn(&currData[(ii+1)*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
			       &propSim->Sd()[ii*nDims*nDims], nDims);
	mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &currData[ii*nDims],
		dT[ii], sqrtDT[ii], currParams);
	propAccept[ii] -= lmvn(&currData[(ii+1)*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
			       &propSim->Sd()[ii*nDims*nDims], nDims);


	// evaluate mh ratio
	if(exp(propAccept[ii]) >= unif_rand()) {
	  currData[ii*nDims + parIndex[0]+jj] = propData[ii*nDims + parIndex[0]+jj];
	  paramAccept[nParams + jj]++;
	}
	else {
	  propData[ii*nDims + parIndex[0]+jj] = currData[ii*nDims + parIndex[0]+jj];
	}
      }
    }
  }
 stop:
  return foundNaN;
}

// loglikelihood evaluation.  this is always done on currData.
double loglik(double params[], PropSim *propSim, double *currData, double *dT, double *sqrtDT) {
  double ll = 0;
  for(int ii = 0; ii < propSim->nComp()-1; ii++) { // *** PARALLELIZABLE FOR-LOOP ***
    mvEuler(&(propSim->Mean()[ii*nDims]), &(propSim->Sd()[ii*nDims*nDims]), &(currData[ii*nDims]),
	    dT[ii], sqrtDT[ii], params);
    ll += lmvn(&(currData[(ii+1)*nDims]), &(propSim->Z()[ii*nDims]), &(propSim->Mean()[ii*nDims]),
	       &(propSim->Sd()[ii*nDims*nDims]), nDims);
  }
  return(ll);
}

// componentwise vanilla MH parameter updates
void paramVnlUpdate(double jumpSD[], Prior *prior, int paramAccept[], PropSim *propSim, double *currData,
		    double *currParams, double *propParams, double *dT, double *sqrtDT) {
  double acc, currLoglik, propLoglik;
  int ii;
  // initialize
  for(ii = 0; ii < nParams; ii++) {
    propParams[ii] = currParams[ii];
  }
  currLoglik = loglik(currParams, propSim, currData, dT, sqrtDT);
  // random walk metropolis updates
  for(ii = 0; ii < nParams; ii++) {
    // proposal
    propParams[ii] = currParams[ii] + jumpSD[ii] * norm_rand();
    // only calculate acceptance if valid
    if(isValidParams(propParams)) {
      propLoglik = loglik(propParams, propSim, currData, dT, sqrtDT);
      acc = propLoglik - currLoglik;

      acc += prior->logPrior(propParams, currData);
      acc -= prior->logPrior(currParams, currData);

      if(exp(acc) >= unif_rand()) {
	currParams[ii] = propParams[ii];
	currLoglik = propLoglik;
	paramAccept[ii]++;
      }
    }
    // propParams and currParams should only differ by one element
    propParams[ii] = currParams[ii];
  }
  return;
}

// multiresolution updates using eraker proposals
// flow is as follows: (1) drop down all low res samples
// (2) for each new missing data point, drop down + calculate acceptance rate
// this flow allows for easy parallelization of the relevant for-loops.
void multiUpdate(int multiInd[], int pastMultiInd[], int nMulti, int multiAccept[],
		 double logMultiAcc[], int nPastComp, int nPastSamples, int nMaxMultiTries,
		 Prior *prior, Prior *pastPrior, PropSim *propSim, double *pastData,
		 double *pastParams, double *currData, double *currParams, double *propData,
		 double *propParams, double *dT, double *sqrtDT, double *B, double *sqrtB,
		 double *pastDT, double *pastSqrtDT) {
  double acc = 0;
  double macc = 0;
  int ii, jj, kk, nTries;
  // pick an old sample
  int pastSmp = (int) floor(unif_rand()*nPastSamples);
  // params
  for(ii = 0; ii < nParams; ii++) {
    propParams[ii] = pastParams[pastSmp*nParams + ii];
  }
  // drop down all old data
  for(ii = 0; ii < nPastComp; ii++) {
    jj = pastMultiInd[ii];
    for(kk = 0; kk < nDims; kk++) {
      propData[jj*nDims + kk] = pastData[pastSmp*nPastComp*nDims + ii*nDims + kk];
    }
  }
  // sample new missing data
  for(jj = 0; jj < nMulti; jj++) { // *** PARALLELIZABLE FOR-LOOP ***
    ii = multiInd[jj];
    // simulate a draw
    mvEraker(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &propData[(ii-1)*nDims],
	     &propData[(ii+1)*nDims], B[ii], sqrtB[ii], propParams);
    nTries = 0;
    for(kk = 0; kk < nDims; kk++) {
      propSim->Z()[ii*nDims + kk] = norm_rand();
    }
    xmvn(&propData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
	 &propSim->Sd()[ii*nDims*nDims], nDims);
    while(!isValidData(&propData[ii*nDims]) && nTries < nMaxMultiTries) {
      for(kk = 0; kk < nDims; kk++) {
	propSim->Z()[ii*nDims + kk] = norm_rand();
      }
      xmvn(&propData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
	   &propSim->Sd()[ii*nDims*nDims], nDims);
      nTries++;
    }
    if(nTries == nMaxMultiTries) {
      goto stop;
    }
    // evaluate the likelihood
    // new proposal miss
    acc = -lmvn(&propData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
		&propSim->Sd()[ii*nDims*nDims], nDims);


    // new proposal obs
    mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &propData[(ii-1)*nDims],
	    pastDT[jj], pastSqrtDT[jj], propParams);
    acc -= lmvn(&propData[(ii+1)*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
		&propSim->Sd()[ii*nDims*nDims], nDims);


    // new target comp
    mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &propData[(ii-1)*nDims], dT[ii-1],
	    sqrtDT[ii-1], propParams);
    acc += lmvn(&propData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
		&propSim->Sd()[ii*nDims*nDims], nDims);

    mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &propData[ii*nDims],
	    dT[ii], sqrtDT[ii], propParams);
    acc += lmvn(&propData[(ii+1)*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
		&propSim->Sd()[ii*nDims*nDims], nDims);


    // old proposal miss
    mvEraker(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &currData[(ii-1)*nDims],
	     &currData[(ii+1)*nDims], B[ii], sqrtB[ii], currParams);
    acc += lmvn(&currData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
		&propSim->Sd()[ii*nDims*nDims], nDims);


    // old proposal obs
    mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &currData[(ii-1)*nDims],
	    pastDT[jj], pastSqrtDT[jj], currParams);
    acc += lmvn(&currData[(ii+1)*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
		&propSim->Sd()[ii*nDims*nDims], nDims);


    // old target comp
    mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &currData[(ii-1)*nDims],
	    dT[ii-1], sqrtDT[ii-1], currParams);
    acc -= lmvn(&currData[ii*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
		&propSim->Sd()[ii*nDims*nDims], nDims);

    mvEuler(&propSim->Mean()[ii*nDims], &propSim->Sd()[ii*nDims*nDims], &currData[ii*nDims],
	    dT[ii], sqrtDT[ii], currParams);
    acc -= lmvn(&currData[(ii+1)*nDims], &propSim->Z()[ii*nDims], &propSim->Mean()[ii*nDims],
		&propSim->Sd()[ii*nDims*nDims], nDims);


    logMultiAcc[jj] = acc;
    macc += acc;
  }

  // parameter priors
  macc += prior->logPrior(propParams, propData);
  macc -= prior->logPrior(currParams, currData);

  macc -= pastPrior->logPrior(propParams, propData);
  macc += pastPrior->logPrior(currParams, currData);


  if(exp(macc) >= unif_rand()) {
    for(ii = 0; ii < nParams; ii++) {
      currParams[ii] = propParams[ii];
    }
    for(ii = 0; ii < propSim->nComp()*nDims; ii++) {
      currData[ii] = propData[ii];
    }
    multiAccept[0]++;
  }
 stop:
  return;
}
