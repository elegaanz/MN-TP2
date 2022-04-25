#include <stdio.h>
#include <stdlib.h>

#include "../include/complexe.h"

#include "flop.h"

#define    NB_FOIS        512

int main (int argc, char **argv)
{
 complexe_float_t c1 = {1.0, 2.0} ;
 complexe_float_t c2 = {3.0, 6.0} ;
 complexe_float_t c3 ;

 complexe_double_t cd1 ;
 complexe_double_t cd2 ;
 complexe_double_t cd3 ;

 struct timeval start, end ;
 
 int i ;

 init_flop_micro () ;
 
 printf ("******************Tests sur les fonctions pour les additions***********************\n") ;

 printf ("c1.r %f c&.i %f\n", c1.real, c1.imaginary) ;

 TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     c3 = add_complexe_float (c1, c2) ;

 TOP_MICRO(end) ;

 printf ("apres boucle c3.real %f c3.imaginary %f duree pour les floats : %f \n", c3.real, c3.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;

 cd1 = (complexe_double_t) {10.0, 7.0} ;
 cd2 = (complexe_double_t) {25.0, 32.0} ;


 init_flop_micro () ;

 TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     cd1 = add_complexe_double (cd1, cd2) ;

 TOP_MICRO(end) ;

 printf ("apres boucle cd3.real %f cd3.imaginary %f duree pour les doubles : %f \n", cd1.real, cd1.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;

 printf ("******************Tests sur les fonctions pour les multiplications***********************\n") ;

 printf ("c1.r %f c1.i %f\n", c1.real, c1.imaginary) ;

 init_flop_micro () ;

 TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     c3 = mult_complexe_float (c1, c2) ;

 TOP_MICRO(end) ;

 printf ("apres boucle c3.real %f c3.imaginary %f duree pour les floats : %f \n", c3.real, c3.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;

 init_flop_micro () ;

 cd1 = (complexe_double_t) {10.0, 7.0} ;
  
 TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     cd3 = mult_complexe_double (cd1, cd2) ;

 TOP_MICRO(end) ;

 printf ("apres boucle cd3.real %f cd3.imaginary %f duree pour les doubles : %f \n", cd1.real, cd1.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;

 printf ("******************Tests sur les fonctions pour les divisions***********************\n") ;

 printf ("c1.r %f c1.i %f\n", c1.real, c1.imaginary) ;

 init_flop_micro () ;

 TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     c3 = div_complexe_float (c1, c2) ;

 TOP_MICRO(end) ;

 printf ("apres boucle c3.real %f c3.imaginary %f duree pour les floats : %f \n", c3.real, c3.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;

 cd1 = (complexe_double_t) {10.0, 7.0} ;
 cd2 = (complexe_double_t) {25.0, 32.0} ;

 init_flop_micro () ;
 
 TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
     cd3 = div_complexe_double (cd1, cd2) ;

 TOP_MICRO(end) ;

 printf ("apres boucle cd3.real %f cd3.imaginary %f duree pour les doubles : %f \n", cd3.real, cd3.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;
 

 exit (0) ;
}
