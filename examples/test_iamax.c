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
 
 int res ;
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

   printf ("===========================Tests Isamax===============================\n \n") ;
  
 printf ("===========================Isamax tsc===============================\n") ;
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     res = 0;
     
     start_tsc = _rdtsc () ;
        res = mnblas_isamax (VECSIZE, vec1, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("isamax nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %d\n", res) ;
 printf ("===========================Isamax micro===============================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     res = 0 ;
     
     TOP_MICRO(start) ;
        res = mnblas_isamax (VECSIZE, vec1, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("isamax micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %d\n", res) ;
 printf ("===========================Isamax nano===============================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     res = 0;
     
     start_tsc = _rdtsc () ;
        res = mnblas_idamax (VECSIZE, vec1_double, 1);
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("idamax nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %d\n", res) ;
   printf ("\n \n===========================Tests Idamax===============================\n \n") ;
  
 printf ("===========================Idamax tsc===============================\n") ;

  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     res = 0;
     
     start_tsc = _rdtsc () ;
        res = mnblas_idamax (VECSIZE, vec1_double, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("idamax tsc ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %d\n", res) ;
 printf ("===========================Idamax micro===============================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     res = 0 ;
     
     TOP_MICRO(start) ;
        res = mnblas_idamax (VECSIZE, vec1_double, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("idamax micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %d\n", res) ;
 printf ("===========================Idamax nano===============================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     res = 0;
     
     start_tsc = _rdtsc () ;
        res = mnblas_idamax (VECSIZE, vec1_double, 1);
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("idamax nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %d\n", res) ;
 printf ("==========================================================\n") ;

    printf ("\n \n===========================Tests Icamax===============================\n \n") ;
  
 printf ("===========================Icamax tsc===============================\n") ;

  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0;
     
     start_tsc = _rdtsc () ;
        res = mnblas_icamax (VECSIZE, C, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("icamax nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %d\n", res) ;
 printf ("===========================Icamax micro===============================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0 ;
     
     TOP_MICRO(start) ;
        res = mnblas_icamax (VECSIZE, C, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("icamax micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %d\n", res) ;
 printf ("===========================Icamax nano===============================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0;
     
     start_tsc = _rdtsc () ;
        res = mnblas_icamax (VECSIZE, C, 1);
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("icamax nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %d\n", res) ;
 printf ("==========================================================\n") ;

   printf ("\n \n===========================Tests Izamax===============================\n \n") ;
  
 printf ("===========================Izamax tsc===============================\n") ;

  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0;
     
     start_tsc = _rdtsc () ;
        res = mnblas_izamax (VECSIZE, C2, 1) ;
     end_tsc = _rdtsc () ;
     
     calcul_flop_tsc ("izamax nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %d\n", res) ;
 printf ("===========================Izamax micro===============================\n") ;
 
 init_flop_micro () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0 ;
     
     TOP_MICRO(start) ;
        res = mnblas_izamax (VECSIZE, C2, 1) ;
     TOP_MICRO(end) ;
     
     calcul_flop_micro ("izamax micro", 2 * VECSIZE, tdiff_micro (&start, &end)) ;
   }

 printf ("res = %d\n", res) ;
 printf ("===========================Izamax nano===============================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     res = 0;
     
     start_tsc = _rdtsc () ;
        res = mnblas_izamax (VECSIZE, C2, 1);
     end_tsc = _rdtsc () ;
     
     calcul_flop_nano ("izamax nano ", 2 * VECSIZE, end_tsc-start_tsc) ;
   }

 printf ("res = %d\n", res) ;
 printf ("==========================================================\n") ;

}
