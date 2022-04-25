#include <stdio.h>
#include <stdlib.h>

#include "../include/complexe2.h"


#define    NB_FOIS        512

#include "flop.h"

int main (int argc, char **argv)
{
 complexe_float_t c1= {1.0, 2.0} ;
 complexe_float_t c2= {3.0, 6.0} ;
 complexe_float_t c3 ;

 complexe_double_t cd1 ;
 complexe_double_t cd2 ;
 complexe_double_t cd3 ;

 unsigned long long int start, end  ;

 int i ;

 init_flop_tsc () ;

 printf ("******************Tests sur les fonctions pour les additions***********************\n") ;

 printf ("c1.r %f c&.i %f\n", c1.real, c1.imaginary) ;

 start = _rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     c3 = add_complexe_float (c1, c2) ;

 end = _rdtsc () ;

 printf ("apres boucle c3.real %f c3.imaginary %f %Ld cycles \n", c3.real, c3.imaginary, end-start) ;

 calcul_flop_micro ("calcul complexe nano ", NB_FOIS*2, end-start) ;

 cd1 = (complexe_double_t) {10.0, 7.0} ;
 cd2 = (complexe_double_t) {25.0, 32.0} ;

 init_flop_tsc () ;

 start = _rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     cd3 = add_complexe_double (cd1, cd2) ;

 end = _rdtsc () ;

 printf ("apres boucle cd3.real %f cd3.imaginary %f %Ld cycles \n", cd3.real, cd3.imaginary, end-start) ;

 calcul_flop_micro ("calcul complexe nano ", NB_FOIS*2, end-start) ;


 printf ("******************Tests sur les fonctions pour les multiplications***********************\n") ;

 printf ("c1.r %f c&.i %f\n", c1.real, c1.imaginary) ;

 init_flop_tsc () ;

 start = _rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     c3 = mult_complexe_float (c1, c2) ;

 end = _rdtsc () ;

 printf ("apres boucle c3.real %f c3.imaginary %f %Ld cycles \n", c3.real, c3.imaginary, end-start) ;

 calcul_flop_micro ("calcul complexe nano ", NB_FOIS*2, end-start) ;

 init_flop_tsc () ;

 start = _rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     cd3 = mult_complexe_double (cd1, cd2) ;

 end = _rdtsc () ;

 printf ("apres boucle cd1.real %f cd1.imaginary %f %Ld cycles \n", cd1.real, cd1.imaginary, end-start) ;

 calcul_flop_micro ("calcul complexe nano ", NB_FOIS*2, end-start) ;


 printf ("******************Tests sur les fonctions pour les divisions***********************\n") ;

 printf ("c1.r %f c&.i %f\n", c1.real, c1.imaginary) ;

 init_flop_tsc () ;

 start = _rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     c3 = div_complexe_float (c1, c2) ;

 end = _rdtsc () ;

 printf ("apres boucle c3.real %f c3.imaginary %f %Ld cycles \n", c3.real, c3.imaginary, end-start) ;

 calcul_flop_micro ("calcul complexe nano ", NB_FOIS*2, end-start) ;

 init_flop_tsc () ;

 start = _rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     cd3 = div_complexe_double (cd1, cd2) ;

 end = _rdtsc () ;

 printf ("apres boucle cd3.real %f cd3.imaginary %f %Ld cycles \n", cd3.real, cd3.imaginary, end-start) ;

 calcul_flop_micro ("calcul complexe nano ", NB_FOIS*2, end-start) ;

 exit (0) ;
 
}


