#include <stdio.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;

vfloat vec1, vec2 ;

typedef double vdouble [VECSIZE] ;

vdouble vec1_double, vec2_double ;

void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

void vector_init_double (vdouble V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print_double (vdouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

int main (int argc, char **argv)
{
  
 struct timeval start, end ;
 unsigned long long int start_tsc, end_tsc ;

  
  complexe_float_t c[VECSIZE];
  complexe_double_t c2[VECSIZE];
  for(int i=0; i<VECSIZE; i++){
    complexe_float_t tmp={0.1*i,0.2*i};
    c[i]=tmp;
  }
  for(int i=0; i<VECSIZE; i++){
    complexe_double_t tmp={0.1*i,0.2*i};
    c2[i]=tmp;
  }
  void *C;
  C=&c[0];
  void *C2;
  C2=&c2[0];
 complexe_double_t dot ={2,4};
 void *DOT;
 DOT = &dot;
 float res ;
 int i ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        mncblas_sdot (VECSIZE, vec1, 1, vec2, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("sdot nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     res = 0.0 ;
     
     TOP_MICRO(start) ;
        mncblas_sdot (VECSIZE, vec1, 1, vec2, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("sdot micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        mncblas_sdot (VECSIZE, vec1, 1, vec2, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("sdot nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

 printf ("\n \n===========================Tests ddot===============================\n \n") ;
  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     vector_init_double (vec2_double, 2.0) ;
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        mncblas_ddot (VECSIZE, vec1_double, 1, vec2_double, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("ddot nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     vector_init_double (vec2_double, 2.0) ;
     res = 0.0 ;
     
     TOP_MICRO(start) ;
        mncblas_ddot (VECSIZE, vec1_double, 1, vec2_double, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("ddot micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     vector_init_double (vec2_double, 2.0) ;
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        mncblas_ddot (VECSIZE, vec1_double, 1, vec2_double, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("ddot nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

   printf ("\n \n===========================Tests cdot===============================\n \n") ;
  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        mncblas_cdotu_sub (VECSIZE,C,1,C2,1,DOT) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("cdot nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     TOP_MICRO(start) ;
        mncblas_cdotu_sub (VECSIZE,C,1,C2,1,DOT) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("cdot micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        mncblas_cdotu_sub (VECSIZE,C,1,C2,1,DOT) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("cdot nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

   printf ("\n \n===========================Tests zdot===============================\n \n") ;
  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        mncblas_zdotu_sub (VECSIZE,C,1,C2,1,DOT) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("zdot nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     TOP_MICRO(start) ;
        mncblas_zdotu_sub (VECSIZE,C,1,C2,1,DOT) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("zdot micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        mncblas_zdotu_sub (VECSIZE,C,1,C2,1,DOT) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("zdot nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
}
