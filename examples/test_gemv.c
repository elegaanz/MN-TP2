#include <stdio.h>

#include "complexe.h"
#include "mnblas.h"

#include "flop.h"

#define VECSIZE 65536

#define NB_FOIS 10

typedef float vfloat[VECSIZE];

vfloat vec1, vec2, vec3;

typedef double vdouble[VECSIZE];

vdouble vec1_double, vec2_double, vec3_double;

void vector_init(vfloat V, float x) {
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = x;

  return;
}

void vector_print(vfloat V) {
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("%f ", V[i]);
  printf("\n");

  return;
}

void vector_init_double(vdouble V, double x) {
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    V[i] = x;

  return;
}

void vector_print_double(vdouble V) {
  register unsigned int i;

  for (i = 0; i < VECSIZE; i++)
    printf("%f ", V[i]);
  printf("\n");

  return;
}

int main(int argc, char **argv) {

  MNCBLAS_LAYOUT layout = MNCblasRowMajor;
  MNCBLAS_TRANSPOSE transpose = MNCblasNoTrans;
  struct timeval start, end;
  unsigned long long int start_tsc, end_tsc;

  int i;

  complexe_float_t c1[VECSIZE];
  complexe_float_t c2[VECSIZE];
  complexe_float_t c3[VECSIZE];
  complexe_double_t cd1[VECSIZE];
  complexe_double_t cd2[VECSIZE];
  complexe_double_t cd3[VECSIZE];
  for (int i = 0; i < VECSIZE; i++) {
    complexe_float_t tmp1 = {0.1 * i, 0.2 * i};
    complexe_float_t tmp2 = {0.2 * i, 0.3 * i};
    complexe_float_t tmp3 = {0.3 * i, 0.4 * i};
    c1[i] = tmp1;
    c2[i] = tmp2;
    c3[i] = tmp3;
  }
  for (int i = 0; i < VECSIZE; i++) {
    complexe_double_t tmp1 = {0.1 * i, 0.2 * i};
    complexe_double_t tmp2 = {0.2 * i, 0.3 * i};
    complexe_double_t tmp3 = {0.3 * i, 0.4 * i};
    cd1[i] = tmp1;
    cd2[i] = tmp2;
    cd3[i] = tmp3;
  }
  void *C1;
  C1 = &c1[0];
  void *C2;
  C2 = &c2[0];
  void *C3;
  C3 = &c3[0];
  void *CD1;
  CD1 = &cd1[0];
  void *CD2;
  CD2 = &cd2[0];
  void *CD3;
  CD3 = &cd3[0];

  complexe_float_t alphac = {3.3, 4.4};
  complexe_float_t betac = {2.3, 2.4};
  complexe_double_t alphaz = {3.3, 4.4};
  complexe_double_t betaz = {2.3, 2.4};

  void *ac;
  ac = &alphac;
  void *bc;
  bc = &betac;
  void *az;
  az = &alphaz;
  void *bz;
  bz = &betaz;

  for (i = 0; i < NB_FOIS; i++) {
    vector_init(vec1, 1.0);
    vector_init(vec2, 2.0);
    vector_init(vec3, 3.0);

    start_tsc = _rdtsc();
    mncblas_sgemv(layout, transpose, 1, VECSIZE, 2, vec1, 1, vec2, 1, 3, vec3,
                  1);
    end_tsc = _rdtsc();

    calcul_flop_tsc("sgemv nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("==========================================================\n");

  init_flop_micro();

  for (i = 0; i < NB_FOIS; i++) {
    vector_init(vec1, 1.0);
    vector_init(vec2, 2.0);
    vector_init(vec3, 3.0);

    TOP_MICRO(start);
    mncblas_sgemv(layout, transpose, 1, VECSIZE, 2, vec1, 1, vec2, 1, 3, vec3,
                  1);
    TOP_MICRO(end);

    calcul_flop_micro("sgemv micro", 2 * VECSIZE, tdiff_micro(&start, &end));
  }

  printf("==========================================================\n");

  init_flop_tsc();

  for (i = 0; i < NB_FOIS; i++) {
    vector_init(vec1, 1.0);
    vector_init(vec2, 2.0);
    vector_init(vec3, 3.0);

    start_tsc = _rdtsc();
    mncblas_sgemv(layout, transpose, 1, VECSIZE, 2, vec1, 1, vec2, 1, 3, vec3,
                  1);
    end_tsc = _rdtsc();

    calcul_flop_nano("sgemv nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("==========================================================\n");

  printf("\n \n===========================Tests "
         "dgemv===============================\n \n");

  for (i = 0; i < NB_FOIS; i++) {
    vector_init_double(vec1_double, 1.0);
    vector_init_double(vec2_double, 2.0);
    vector_init_double(vec3_double, 3.0);

    start_tsc = _rdtsc();
    mncblas_dgemv(layout, transpose, 1, VECSIZE, 2, vec1_double, 1, vec2_double,
                  1, 3, vec3_double, 1);
    end_tsc = _rdtsc();

    calcul_flop_tsc("dgemv nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("==========================================================\n");

  init_flop_micro();

  for (i = 0; i < NB_FOIS; i++) {
    vector_init_double(vec1_double, 1.0);
    vector_init_double(vec2_double, 2.0);
    vector_init_double(vec3_double, 3.0);

    TOP_MICRO(start);
    mncblas_dgemv(layout, transpose, 1, VECSIZE, 2, vec1_double, 1, vec2_double,
                  1, 3, vec3_double, 1);
    TOP_MICRO(end);

    calcul_flop_micro("dgemv micro", 2 * VECSIZE, tdiff_micro(&start, &end));
  }

  printf("==========================================================\n");

  init_flop_tsc();

  for (i = 0; i < NB_FOIS; i++) {
    vector_init_double(vec1_double, 1.0);
    vector_init_double(vec2_double, 2.0);
    vector_init_double(vec3_double, 3.0);

    start_tsc = _rdtsc();
    mncblas_dgemv(layout, transpose, 1, VECSIZE, 2, vec1_double, 1, vec2_double,
                  1, 3, vec3_double, 1);
    end_tsc = _rdtsc();

    calcul_flop_nano("dgemv nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("==========================================================\n");

  printf("\n \n===========================Tests "
         "cgemv===============================\n \n");

  for (i = 0; i < NB_FOIS; i++) {
    start_tsc = _rdtsc();
    mncblas_cgemv(layout, transpose, 1, VECSIZE, ac, C1, 1, C2, 1, bc, C3, 1);
    end_tsc = _rdtsc();

    calcul_flop_tsc("cgemv nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("==========================================================\n");

  init_flop_micro();

  for (i = 0; i < NB_FOIS; i++) {

    TOP_MICRO(start);
    mncblas_cgemv(layout, transpose, 1, VECSIZE, ac, C1, 1, C2, 1, bc, C3, 1);
    TOP_MICRO(end);

    calcul_flop_micro("cgemv micro", 2 * VECSIZE, tdiff_micro(&start, &end));
  }

  printf("==========================================================\n");

  init_flop_tsc();

  for (i = 0; i < NB_FOIS; i++) {

    start_tsc = _rdtsc();
    mncblas_cgemv(layout, transpose, 1, VECSIZE, ac, C1, 1, C2, 1, bc, C3, 1);
    end_tsc = _rdtsc();

    calcul_flop_nano("cgemv nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("==========================================================\n");

  printf("\n \n===========================Tests "
         "zgemv===============================\n \n");

  for (i = 0; i < NB_FOIS; i++) {
    start_tsc = _rdtsc();
    mncblas_zgemv(layout, transpose, 1, VECSIZE, az, CD1, 1, CD2, 1, bz, CD3,
                  1);
    end_tsc = _rdtsc();

    calcul_flop_tsc("zgemv nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("==========================================================\n");

  init_flop_micro();

  for (i = 0; i < NB_FOIS; i++) {

    TOP_MICRO(start);
    mncblas_zgemv(layout, transpose, 1, VECSIZE, az, CD1, 1, CD2, 1, bz, CD3,
                  1);
    TOP_MICRO(end);

    calcul_flop_micro("zgemv micro", 2 * VECSIZE, tdiff_micro(&start, &end));
  }

  printf("==========================================================\n");

  init_flop_tsc();

  for (i = 0; i < NB_FOIS; i++) {

    start_tsc = _rdtsc();
    mncblas_zgemv(layout, transpose, 1, VECSIZE, az, CD1, 1, CD2, 1, bz, CD3,
                  1);
    end_tsc = _rdtsc();

    calcul_flop_nano("zgemv nano ", 2 * VECSIZE, end_tsc - start_tsc);
  }

  printf("==========================================================\n");
}