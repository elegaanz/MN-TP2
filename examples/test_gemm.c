#include "flop.h"
#include "mnblas.h"
#include <stdio.h>

#define ITER 10000
#define M 10
#define N 20
#define K 5

float *zero_matrix(int m, int n) {
  float *matrix = malloc(sizeof(float) * m * n);
  for (int i = 0; i < m * n; i++) {
    matrix[i] = 0.0;
  }
  return matrix;
}

int main() {
  struct timeval start, end;

  float *A = zero_matrix(K, M);
  float *B = zero_matrix(N, K);
  float *C = zero_matrix(M, N);

  init_flop_micro();
  for (int iter = 0; iter < ITER; iter++) {
    TOP_MICRO(start);

    mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K,
                  12.0, A, M, B, K, 42.0, C, N);
    TOP_MICRO(end);
  }
  calcul_flop_micro("gemm micro", M * N * K * ITER, tdiff_micro(&start, &end));
}