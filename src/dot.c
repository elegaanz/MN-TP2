#include "complexe.h"
#include "mnblas.h"
#include <stdio.h>

/*
float mncblas_sdot(const int N, const float *X, const int incX,
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}
*/

float mncblas_sdot(const int N, const float *X, const int incX, const float *Y,
                   const int incY) {
  register unsigned int i = 0;
  register unsigned int j = 0;
  float dot = 0.0;

  for (i = 0; i < N; i += incX) {
    dot += X[i] * Y[j];
    j += incY;
  }

  return dot;
}

double mncblas_ddot(const int N, const double *X, const int incX,
                    const double *Y, const int incY) {
  register unsigned int i = 0;
  register unsigned int j = 0;
  double dot = 0.0;

  for (i = 0; i < N; i += incX) {
    dot += X[i] * Y[j];
    j += incY;
  }

  return dot;
}

void mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu) {
  complexe_float_t *PX = (complexe_float_t *)X;
  complexe_float_t *PY = (complexe_float_t *)Y;
  register unsigned int i = 0;
  register unsigned int j = 0;
  complexe_float_t *dot = (complexe_float_t *)dotu;

  for (i = 0; i < N; i += incX) {
    *dot = add_complexe_float(*dot, mult_complexe_float(PX[i], PY[j]));
    j += incY;
  }
}

void mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc) {
  complexe_float_t *PX = (complexe_float_t *)X;
  PX->imaginary = -PX->imaginary;
  complexe_float_t *PY = (complexe_float_t *)Y;
  register unsigned int i = 0;
  register unsigned int j = 0;
  complexe_float_t *dot = (complexe_float_t *)dotc;

  for (i = 0; i < N; i += incX) {
    *dot = add_complexe_float(*dot, mult_complexe_float(PX[i], PY[j]));
    j += incY;
  }
}

void mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu) {
  complexe_double_t *PX = (complexe_double_t *)X;
  complexe_double_t *PY = (complexe_double_t *)Y;
  register unsigned int i = 0;
  register unsigned int j = 0;
  complexe_double_t *dot = (complexe_double_t *)dotu;

  for (i = 0; i < N; i += incX) {
    *dot = add_complexe_double(*dot, mult_complexe_double(PX[i], PY[j]));
    j += incY;
  }
}

void mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc) {
  complexe_double_t *PX = (complexe_double_t *)X;
  complexe_double_t *PY = (complexe_double_t *)Y;
  PX->imaginary = -PX->imaginary;
  register unsigned int i = 0;
  register unsigned int j = 0;
  complexe_double_t *dot = (complexe_double_t *)dotc;

  for (i = 0; i < N; i += incX) {
    *dot = add_complexe_double(*dot, mult_complexe_double(PX[i], PY[j]));
    j += incY;
  }
}
