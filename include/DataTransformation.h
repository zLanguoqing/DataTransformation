#ifndef _DataTransformation_H_
#define _DataTransformation_H_
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#define MIN (1e-30)
void FFT(double pr[],double pi[],int n, int k,double fr[],double fi[],int l,int il);
void  Walsh(double p[],int n,int k,double x[]);
#endif