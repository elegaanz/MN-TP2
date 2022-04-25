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
 int i ;
  complexe_float_t alphaF= {1.0, 2.0} ;
  void *af;
 af=&alphaF;
  complexe_double_t alphaD= {1.0, 2.0} ;
  void *ad;
 ad=&alphaD;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     
     start_tsc = _rdtsc () ;
        mnblas_saxpy (VECSIZE, 12.0, vec1, 1, vec2, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("saxpy tsc ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     
     TOP_MICRO(start) ;
        mnblas_saxpy (VECSIZE, 12.0, vec1, 1, vec2, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("saxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     
     start_tsc = _rdtsc () ;
        mnblas_saxpy (VECSIZE, 12.0, vec1, 1, vec2, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("saxpy nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("==========================================================\n") ;

     printf ("\n \n===========================Tests daxpy===============================\n \n") ;

   init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     vector_init_double (vec2_double, 2.0) ;
     
     start_tsc = _rdtsc () ;
        mnblas_daxpy (VECSIZE, 12.0, vec1_double, 1, vec2_double, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("daxpy tsc ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     vector_init_double (vec2_double, 2.0) ;
     
     TOP_MICRO(start) ;
        mnblas_daxpy (VECSIZE, 12.0, vec1_double, 1, vec2_double, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("daxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     vector_init_double (vec2_double, 2.0) ;
     
     start_tsc = _rdtsc () ;
        mnblas_daxpy (VECSIZE, 12.0, vec1_double, 1, vec2_double, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("daxpy nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("==========================================================\n") ;

     printf ("\n \n===========================Tests caxpy===============================\n \n") ;

   init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {  
     start_tsc = _rdtsc () ;
        mnblas_caxpy (VECSIZE, af, C, 1, C, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("caxpy tsc ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     TOP_MICRO(start) ;
        mnblas_caxpy (VECSIZE, af, C, 1, C, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("caxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     start_tsc = _rdtsc () ;
        mnblas_caxpy (VECSIZE, af, C, 1, C, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("caxpy nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("==========================================================\n") ;

     printf ("\n \n===========================Tests zaxpy===============================\n \n") ;

   init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   { 
     start_tsc = _rdtsc () ;
        mnblas_zaxpy (VECSIZE, ad, C2, 1, C2, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("zaxpy tsc ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     TOP_MICRO(start) ;
        mnblas_zaxpy (VECSIZE, ad, C2, 1, C2, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("zaxpy micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   { 
     start_tsc = _rdtsc () ;
        mnblas_zaxpy (VECSIZE, ad, C2, 1, C2, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("zaxpy nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("==========================================================\n") ;
}
