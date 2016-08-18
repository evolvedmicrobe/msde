#ifndef Priors_h
#define Priors_h 1

#include <Rcpp.h>
using namespace Rcpp;

//base prior class
class Prior {
 public:
  enum Type {
    Flat = 1,
    Normal = 2,
    GCop = 3,
    Custom = 4
  };

  //pure virtual function must be defined in each derived class
  virtual double logPrior(double params[], double x[])=0;
};


class FlatPrior : public Prior
{
 public:
  FlatPrior();
  ~FlatPrior();
  double logPrior(double params[], double x[]);
};



class NormalPrior : public Prior
{
  double *mean;
  double *cholSd;
  int parindex0, nMiss0;
 public:
  NormalPrior(List priorParams, int index, int pnMiss0);
  ~NormalPrior();
  double logPrior(double params[], double x[]);
};



class GCopPrior : public Prior
{
  // gaussian copula parameters
  int *nBreaks;
  double *range, *dx, *pdf, *logPdf, *cdf;
  double *mean, *sd, *RhoCholSd;
  int parindex0, nMiss0;
 public:
  GCopPrior(List priorParams, int index, int pnMiss0);
  ~GCopPrior();
  // parameter + unobserved first states gaussian copula prior.
  double logPrior(double params[], double x[]);
};



class CustomPrior : public Prior {
  // custom prior parameters
  double *customParams;
  int parindex0, nMiss0;
 public:
  CustomPrior(List priorParams, int index, int pnMiss0);
  ~CustomPrior();
  double logPrior(double params[], double x[]);
};
#endif
