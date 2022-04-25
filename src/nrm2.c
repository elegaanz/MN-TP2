#include <stdlib.h>
#include "mnblas.h"
#include "complexe.h"
#include <math.h> 

float  mnblas_snrm2(const int N, const float *X, const int incX){
	float res = 0.0;
	for (int i = 0; i < N; i += incX) {
		res += X[i] * X[i];
	}
	return sqrt(res);
}

double mnblas_dnrm2(const int N, const double *X, const int incX){
	double res = 0.0;
	for (int i = 0; i < N; i += incX) {
		res += X[i] * X[i];
	}
	return sqrt(res);
}

float  mnblas_scnrm2(const int N, const void *X, const int incX){
	float res = 0.0;
  complexe_float_t *PX = (complexe_float_t *)X;
	for (int i = 0; i < N; i += incX) {
		res += PX[i].real * PX[i].real + PX[i].imaginary * PX[i].imaginary;
	}
	return sqrt(res);
}

double mnblas_dznrm2(const int N, const void *X, const int incX){
  double res = 0.0;
  complexe_double_t *PX = (complexe_double_t *)X;
	for (int i = 0; i < N; i += incX) {
		res += PX[i].real * PX[i].real + PX[i].imaginary * PX[i].imaginary;
	}
	return sqrt(res);
}
