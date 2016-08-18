//#include "Priors.h"
//#include "sdeCore.h"
//#include "UserDefinitionSpace.h"

#include <Rcpp.h>
using namespace Rcpp;

FlatPrior::FlatPrior() {
  //do nothing
}

FlatPrior::~FlatPrior() {
  //do nothing
}

double FlatPrior::logPrior(double params[], double x[]) {
  return 0.0;
}


NormalPrior::NormalPrior(List priorParams, int index, int pnMiss0) {
  mean = REAL(priorParams["Mu"]);
  cholSd = REAL(priorParams["V"]);
  parindex0 = index;
  nMiss0 = pnMiss0;
}

NormalPrior::~NormalPrior() {
  //do nothing
}

double NormalPrior::logPrior(double params[], double x[]) {
  double tmpZ[nMaxPriorParams];
  double tmpX[nMaxPriorParams];
  double lp;
  int ii;
  for(ii = 0; ii < nParams; ii++) {
    tmpX[ii] = params[ii];
  }
  for(ii = 0; ii < nMiss0; ii++) {
    tmpX[nParams+ii] = x[parindex0+ii];
  }
  lp = lmvn(tmpX, tmpZ, mean, cholSd, nParams+nMiss0);
  return(lp);
}

GCopPrior::GCopPrior(List priorParams, int index, int pnMiss0) {
  nBreaks = INTEGER(priorParams["nbreaks"]);
  range = REAL(priorParams["rx"]);
  dx = REAL(priorParams["dx"]);
  pdf = REAL(priorParams["dens.y"]);
  logPdf = REAL(priorParams["ldens.y"]);
  cdf = REAL(priorParams["Dens.y"]);
  mean = REAL(priorParams["mean"]);
  sd = REAL(priorParams["sd"]);
  RhoCholSd = REAL(priorParams["Rho"]);
  parindex0 = index;
  nMiss0 = pnMiss0;
}

GCopPrior::~GCopPrior() {
  //do nothing
}

// parameter + unobserved first states gaussian copula prior.
double GCopPrior::logPrior(double params[], double x[]) {
  double qNorm[nMaxPriorParams];
  double tmpX[nMaxPriorParams];
  int n = nParams+nMiss0;
  int ii, jj, colI, densElt, start;
  double lp = 0.0;
  double tmpSum;
  for(ii = 0; ii < nParams; ii++)	{
    tmpX[ii] = params[ii];
  }
  for(ii = 0; ii < nMiss0; ii++) {
    tmpX[nParams+ii] = x[parindex0+ii];
  }
  // normal quantiles and marginal components
  start = 0;
  for(ii = 0; ii < n; ii ++) {
    densElt = (int) floor((tmpX[ii]-range[2*ii])/dx[ii]);
    if((densElt >= 0) & (densElt < nBreaks[ii])) {
      lp += logPdf[densElt + start];
      qNorm[ii] = tmpX[ii] - (range[2*ii] + densElt * dx[ii]);
      qNorm[ii] *= pdf[densElt + start];
      qNorm[ii] = Rf_qnorm5(cdf[densElt + start + ii] +  qNorm[ii], 0.0, 1.0, 1, 0);
    }
    else {
      lp += Rf_dnorm4(tmpX[ii], mean[ii], sd[ii], 1);
      qNorm[ii] = (tmpX[ii] - mean[ii])/sd[ii];
    }
    start += nBreaks[ii];
  }
  // copula components
  // iid standard normal densities
  tmpSum = 0.0;
  for(ii = 0; ii < n; ii++) {
    tmpSum += qNorm[ii] * qNorm[ii];
  }
  lp += 0.5 * tmpSum;
  // multivariate normal density
  for(ii = 0; ii < n; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) {
      tmpSum += RhoCholSd[colI + jj] * qNorm[jj];
    }
    qNorm[ii] = (qNorm[ii] - tmpSum)/RhoCholSd[colI + ii];
  }
  tmpSum = 0.0;
  for(ii = 0; ii < n; ii++) {
    tmpSum += qNorm[ii] * qNorm[ii];
  }
  tmpSum *= 0.5;
  for(ii = 0; ii < n; ii++) {
    tmpSum += log(RhoCholSd[n*ii + ii]);
  }
  lp -= tmpSum;
  return(lp);
}

// custom prior
CustomPrior::CustomPrior(List priorParams, int index, int pnMiss0) {
  customParams = new double[priorParams.size()];
  for(int ii = 0; ii < priorParams.size(); ii++) {
    customParams[ii] = *(REAL(priorParams[ii]));
  }
  parindex0 = index;
  nMiss0 = pnMiss0;
}

CustomPrior::~CustomPrior() {
  delete [] customParams;
}
