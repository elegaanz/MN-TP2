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
 
 float res ;
 int i ;

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
  
 printf ("=============================Tests srnm2=============================\n") ;
  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        res = mnblas_snrm2 (VECSIZE, vec1, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("snrm2 nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     res = 0.0 ;
     
     TOP_MICRO(start) ;
        res = mnblas_snrm2 (VECSIZE, vec1, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("snrm2 micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        res = mnblas_snrm2 (VECSIZE, vec1, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("snrm2 nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

   printf ("=============================Tests dnrm2=============================\n") ;
  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        res = mnblas_dnrm2 (VECSIZE, vec1_double, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("dnrm2 nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     res = 0.0 ;
     
     TOP_MICRO(start) ;
        res = mnblas_dnrm2 (VECSIZE, vec1_double, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("dnrm2 micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        res = mnblas_dnrm2 (VECSIZE, vec1_double, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("dnrm2 nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

     printf ("=============================Tests cnrm2=============================\n") ;
  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        res = mnblas_scnrm2 (VECSIZE, C, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("scnrm2 nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     TOP_MICRO(start) ;
        res = mnblas_scnrm2 (VECSIZE, C, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("scnrm2 micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {

     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        res = mnblas_scnrm2 (VECSIZE, C, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("scnrm2 nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

  printf ("=============================Tests znrm2=============================\n") ;
  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        res = mnblas_dznrm2 (VECSIZE, C2, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("dznrm2 nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0.0 ;
     
     TOP_MICRO(start) ;
        res = mnblas_dznrm2 (VECSIZE, C2, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("dznrm2 micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {

     res = 0.0 ;
     
     start_tsc = _rdtsc () ;
        res = mnblas_dznrm2 (VECSIZE, C2, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("dznrm2 nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %f\n", res) ;
 printf ("==========================================================\n") ;
}