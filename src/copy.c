#include "complexe.h"
#include "mnblas.h"

void mncblas_scopy(const int N, const float *X, const int incX, float *Y,
                   const int incY) {
  register unsigned int i = 0;
  register unsigned int j = 0;

  for (; ((i < N) && (j < N)); i += incX, j += incY) {
    Y[j] = X[i];
  }

  return;
}

void mncblas_dcopy(const int N, const double *X, const int incX, double *Y,
                   const int incY) {
  register unsigned int i = 0;
  register unsigned int j = 0;

  for (; ((i < N) && (j < N)); i += incX, j += incY) {
    Y[j] = X[i];
  }

  return;
}

void mncblas_ccopy(const int N, const void *X_, const int incX, void *Y_,
                   const int incY) {

  complexe_float_t *X = (complexe_float_t *)X_;
  complexe_float_t *Y = (complexe_float_t *)Y_;
  register unsigned int i = 0;
  register unsigned int j = 0;

  for (; ((i < N) && (j < N)); i += incX, j += incY) {
    Y[j] = X[i];
  }

  return;
}

void mncblas_zcopy(const int N, const void *X_, const int incX, void *Y_,
                   const int incY) {

  complexe_double_t *X = (complexe_double_t *)X_;
  complexe_double_t *Y = (complexe_double_t *)Y_;
  register unsigned int i = 0;
  register unsigned int j = 0;

  for (; ((i < N) && (j < N)); i += incX, j += incY) {
    Y[j] = X[i];
  }

  return;
}
