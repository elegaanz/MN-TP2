#include "../include/complexe.h"
#include "../include/mnblas.h"
#include <stdio.h>

void mncblas_sgemv(const MNCBLAS_LAYOUT layout, const MNCBLAS_TRANSPOSE TransA,
                   const int M, const int N, const float alpha, const float *A,
                   const int lda, const float *X, const int incX,
                   const float beta, float *Y, const int incY) {
  float tmp[M];
  int tmpj = 0;
  float res;
#pragma omp parallel for
  for (int tmpi = 0; tmpi < M; tmpi++) {
    res = 0;
    tmpj = 0;
    while (tmpj < N) {
      res = A[tmpi * M + tmpj] * X[tmpj];
      tmpj++;
    }
    tmp[tmpi] = res;
    tmpi++;
  }

#pragma omp parallel for
  for (int tmpi = 0; tmpi < M; tmpi++) {
    Y[tmpi] = alpha * tmp[tmpi] + beta * Y[tmpi];
    tmpi++;
  }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M,
                   const int N, const double alpha, const double *A,
                   const int lda, const double *X, const int incX,
                   const double beta, double *Y, const int incY) {

  double tmp[M];
  int tmpj = 0;
  float res;
#pragma omp parallel for
  for (int tmpi = 0; tmpi < M; tmpi++) {
    res = 0;
    tmpj = 0;
    while (tmpj < N) {
      res = A[tmpi * M + tmpj] * X[tmpj];
      tmpj++;
    }
    tmp[tmpi] = res;
    tmpi++;
  }
#pragma omp parallel for
  for (int tmpi = 0; tmpi < M; tmpi++) {
    Y[tmpi] = alpha * tmp[tmpi] + beta * Y[tmpi];
    tmpi++;
  }
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M,
                   const int N, const void *alpha_, const void *A,
                   const int lda, const void *X, const int incX,
                   const void *beta_, void *Y, const int incY) {

  complexe_float_t *PX = (complexe_float_t *)X;
  complexe_float_t *PY = (complexe_float_t *)Y;
  complexe_float_t *PA = (complexe_float_t *)A;
  complexe_float_t *alpha = (complexe_float_t *)alpha_;
  complexe_float_t *beta = (complexe_float_t *)beta_;

  complexe_float_t tmp[M];
  int tmpj = 0;
  complexe_float_t res;
#pragma omp parallel for
  for (int tmpi = 0; tmpi < M; tmpi++) {
    res.real = 0;
    res.imaginary = 0;
    tmpj = 0;
    while (tmpj < N) {
      res = add_complexe_float(
          res, mult_complexe_float(PA[tmpi * M + tmpj], PX[tmpj]));
      tmpj++;
    }
    tmp[tmpi] = res;
    tmpi++;
  }
#pragma omp parallel for
  for (int tmpi = 0; tmpi < M; tmpi++) {
    PY[tmpi] = add_complexe_float(mult_complexe_float(*alpha, tmp[tmpi]),
                                  mult_complexe_float(*beta, PY[tmpi]));
    tmpi++;
  }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA, const int M,
                   const int N, const void *alpha_, const void *A,
                   const int lda, const void *X, const int incX,
                   const void *beta_, void *Y, const int incY) {

  complexe_double_t *PX = (complexe_double_t *)X;
  complexe_double_t *PY = (complexe_double_t *)Y;
  complexe_double_t *PA = (complexe_double_t *)A;

  complexe_double_t *alpha = (complexe_double_t *)alpha_;
  complexe_double_t *beta = (complexe_double_t *)beta_;

  complexe_double_t tmp[M];
  int tmpj = 0;
  complexe_double_t res;
#pragma omp parallel for
  for (int tmpi = 0; tmpi < M; tmpi++) {
    res.real = 0;
    res.imaginary = 0;
    tmpj = 0;
    while (tmpj < N) {
      res = add_complexe_double(
          res, mult_complexe_double(PA[tmpi * M + tmpj], PX[tmpj]));
      tmpj++;
    }
    tmp[tmpi] = res;
    tmpi++;
  }
#pragma omp parallel for
  for (int tmpi = 0; tmpi < M; tmpi++) {
    PY[tmpi] = add_complexe_double(mult_complexe_double(*alpha, tmp[tmpi]),
                                   mult_complexe_double(*beta, PY[tmpi]));
    tmpi++;
  }
}