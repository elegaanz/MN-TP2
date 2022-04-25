#include "../include/complexe.h"
#include "../include/mnblas.h"

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const float alpha, const float *A,
                   const int lda, const float *B, const int ldb,
                   const float beta, float *C, const int ldc) {

#pragma omp parallel for
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      float AB = 0.0;
      for (int k = 0; k < K; k++) {
        AB += A[i * K + k] * B[k * N + j];
      }
      C[i * M + j] = alpha * AB + beta * C[i * M + j];
    }
  }
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const double alpha, const double *A,
                   const int lda, const double *B, const int ldb,
                   const double beta, double *C, const int ldc) {
#pragma omp parallel for
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      double AB = 0.0;
      for (int k = 0; k < K; k++) {
        AB += A[i * K + k] * B[k * N + j];
      }
      C[i * M + j] = alpha * AB + beta * C[i * M + j];
    }
  }
}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const void *alpha_, const void *A_,
                   const int lda, const void *B_, const int ldb,
                   const void *beta_, void *C_, const int ldc) {
  complexe_float_t *A = (complexe_float_t *)A_;
  complexe_float_t *B = (complexe_float_t *)B_;
  complexe_float_t *C = (complexe_float_t *)C_;
  complexe_float_t *alpha = (complexe_float_t *)alpha_;
  complexe_float_t *beta = (complexe_float_t *)beta_;
#pragma omp parallel for
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      complexe_float_t AB = {0.0, 0.0};
      for (int k = 0; k < K; k++) {
        AB = add_complexe_float(mult_complexe_float(A[i * K + k], B[k * N + j]),
                                AB);
      }
      C[i * M + j] =
          add_complexe_float(mult_complexe_float(*alpha, AB),
                             mult_complexe_float(*beta, C[i * M + j]));
    }
  }
}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const void *alpha_, const void *A_,
                   const int lda, const void *B_, const int ldb,
                   const void *beta_, void *C_, const int ldc) {
  complexe_double_t *A = (complexe_double_t *)A_;
  complexe_double_t *B = (complexe_double_t *)B_;
  complexe_double_t *C = (complexe_double_t *)C_;
  complexe_double_t *alpha = (complexe_double_t *)alpha_;
  complexe_double_t *beta = (complexe_double_t *)beta_;
#pragma omp parallel for
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      complexe_double_t AB = {0.0, 0.0};
      for (int k = 0; k < K; k++) {
        AB = add_complexe_double(
            mult_complexe_double(A[i * K + k], B[k * N + j]), AB);
      }
      C[i * M + j] =
          add_complexe_double(mult_complexe_double(*alpha, AB),
                              mult_complexe_double(*beta, C[i * M + j]));
    }
  }
}