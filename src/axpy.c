#include "complexe.h"
#include "mnblas.h"

void mnblas_saxpy(const int N, const float alpha, const float *X,
                  const int incX, float *Y, const int incY) {

  register unsigned int j = 0;
#pragma omp parallel for
  for (register unsigned int i = 0; i < N; i += incX) {
    Y[j] = alpha * X[i] + Y[j];
    j += incY;
  }

  return;
}

void mnblas_daxpy(const int N, const double alpha, const double *X,
                  const int incX, double *Y, const int incY) {
  register unsigned int j = 0;

#pragma omp parallel for
  for (register unsigned int i = 0; i < N; i += incX) {
    Y[j] = alpha * X[i] + Y[j];
    j += incY;
  }

  return;
}

void mnblas_caxpy(const int N, const void *alpha_, const void *X_,
                  const int incX, void *Y_, const int incY) {

  complexe_float_t *alpha = (complexe_float_t *)alpha_;
  complexe_float_t *X = (complexe_float_t *)X_;
  complexe_float_t *Y = (complexe_float_t *)Y_;
  int j = 0;
#pragma omp parallel for

  for (int i = 0; i < N; i += incX) {
    *Y = add_complexe_float(*Y, mult_complexe_float(*alpha, *X));
    j += incY;
  }
}

void mnblas_zaxpy(const int N, const void *alpha_, const void *X_,
                  const int incX, void *Y_, const int incY) {
  complexe_double_t *alpha = (complexe_double_t *)alpha;
  complexe_double_t *X = (complexe_double_t *)X_;
  complexe_double_t *Y = (complexe_double_t *)Y_;
  int j = 0;
#pragma omp parallel for

  for (int i = 0; i < N; i += incX) {
    *Y = add_complexe_double(*Y, mult_complexe_double(*alpha, *X));
    j += incY;
  }
}
