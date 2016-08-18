//Convenience functions that can be used in user defined functions located in the UserDefinitionSpace
#ifndef ConvenienceFunctions_h
#define ConvenienceFunctions_h 1

int inRange(double x[], double lBound[], double uBound[], int n);
int isNaN(double x[], int n);
void chol_decomp(double U[], double A[], int n);
void xmvn(double x[], double z[], double mean[], double cholSd[], int n);
void zmvn(double z[], double x[], double mean[], double cholSd[], int n, int nMax);
double lmvn(double x[], double z[], double mean[], double cholSd[], int n);

#endif
