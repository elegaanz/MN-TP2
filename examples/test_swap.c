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
  
 printf ("=============================Tests sspwap=============================\n") ;
  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     
     start_tsc = _rdtsc () ;
        mncblas_sswap (VECSIZE, vec1, 1, vec2, 1) ;
     end_tsc = _rdtsc () ;

     long  data=4*(long)VECSIZE*NB_FOIS;//4:float
     printf("%lld Mo/s\n",((end_tsc-start_tsc)*data)/1600000000000000);
     //160000... car durée cycle cpu*10^-6 pour avoir en Mo
   }

 printf ("==========================================================\n") ;

 printf ("=============================Tests dspwap=============================\n") ;
  
 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init_double (vec1_double, 1.0) ;
     vector_init_double (vec2_double, 2.0) ;
     
     start_tsc = _rdtsc () ;
       mncblas_dswap (VECSIZE, vec1_double, 1, vec2_double, 1) ;
     end_tsc = _rdtsc () ;

     long  data=8*(long)VECSIZE*NB_FOIS;//8:double
     printf("%lld Mo/s\n",((end_tsc-start_tsc)*data)/1600000000000000);
     //160000... car durée cycle cpu*10^-6 pour avoir en Mo
   }

 printf ("=============================Tests cspwap=============================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     
     start_tsc = _rdtsc () ;
        mncblas_cswap (VECSIZE, C, 1, C, 1) ;
     end_tsc = _rdtsc () ;

     long  data=8*2*(long)VECSIZE*NB_FOIS;//8*2:deux complexes
     printf("%lld Mo/s\n",((end_tsc-start_tsc)*data)/1600000000000000);
     //160000... car durée cycle cpu*10^-6 pour avoir en Mo
     
 printf ("=============================Tests zspwap=============================\n") ;

 init_flop_tsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     
     start_tsc = _rdtsc () ;
        mncblas_zswap (VECSIZE, C2, 1, C2, 1) ;
     end_tsc = _rdtsc () ;

     long  data=8*2*(long)VECSIZE*NB_FOIS;//8*2:deux complexes
     printf("%lld Mo/s\n",((end_tsc-start_tsc)*data)/1600000000000000);
     //160000... car durée cycle cpu*10^-6 pour avoir en Mo
}