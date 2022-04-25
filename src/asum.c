#include "complexe.h"
#include "mnblas.h"
#include <omp.h>

float mnblas_sasum(const int N, const float *X, const int incX) {
  float res = 0.0;

#pragma omp parallel for reduction(+ : res)
  for (register unsigned int i = 0; i < N; i += incX) {
    if (X[i] < 0)
      res += (-1) * X[i];
    else
      res += X[i];
  }

  return res;
}

float mnblas_scasum(const int N, const void *X, const int incX) {
  complexe_float_t *PX = (complexe_float_t *)X;
  float res = 0.0;

#pragma omp parallel for reduction(+ : res)
  for (register unsigned int i = 0; i < N; i += incX) {
    if (PX[i].imaginary < 0)
      res += (-1) * PX[i].imaginary;
    else
      res += PX[i].imaginary;

    if (PX[i].real < 0)
      res += (-1) * PX[i].real;
    else
      res += PX[i].real;
  }

  return res;
}

double mnblas_dasum(const int N, const double *X, const int incX) {
  float res = 0.0;

#pragma omp parallel for reduction(+ : res)
  for (register unsigned int i = 0; i < N; i += incX) {
    if (X[i] < 0)
      res += (-1) * X[i];
    else
      res += X[i];
  }

  return res;
}

double mnblas_dzasum(const int N, const void *X, const int incX) {
  complexe_float_t *PX = (complexe_float_t *)X;
  double res = 0.0;

#pragma omp parallel for reduction(+ : res)
  for (register unsigned int i = 0; i < N; i += incX) {
    if (PX[i].imaginary < 0)
      res += (-1) * PX[i].imaginary;
    else
      res += PX[i].imaginary;

    if (PX[i].real < 0)
      res += (-1) * PX[i].real;
    else
      res += PX[i].real;
  }

  return res;
}