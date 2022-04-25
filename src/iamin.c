
#include <stdlib.h>
#include "mnblas.h"
#include "complexe.h"

#define abs(x) ((x < 0) ? -x : x)

CBLAS_INDEX mnblas_isamin(const int N, const float  *X, const int incX){
  register unsigned int i = 0 ;
  CBLAS_INDEX res = 0; 
  float minValue = abs(X[0]);
  
  for (; i < N ; i += incX)
    {	
      if (minValue > abs(X[i])){
        minValue = abs(X[i]);
        res = i;
        }
    }

  return res;
}

CBLAS_INDEX mnblas_idamin(const int N, const double *X, const int incX){
  register unsigned int i = 0 ;
  CBLAS_INDEX res = 0; 
  double minValue = abs(X[0]);
  
  for (; i < N ; i += incX)
    {	
      if (minValue > abs(X[i])){
        minValue = abs(X[i]);
        res = i;
        }
    }

  return res;
}

CBLAS_INDEX mnblas_icamin(const int N, const void   *X_, const int incX){
	register unsigned int i = 0 ;
  CBLAS_INDEX res = 0; 
	complexe_float_t *X = (complexe_float_t *)X_;
  float minValue = abs(X[0].real) + abs(X[0].imaginary);
  
  for (; i < N ; i += incX)
    {
      if (minValue > (abs(X[i].real) + abs(X[i].imaginary))){
        minValue = (abs(X[i].real) + abs(X[i].imaginary));
        res = i;
        }
    }

  return res;
  
}

CBLAS_INDEX mnblas_izamin(const int N, const void   *X_, const int incX){
  register unsigned int i = 0 ;
  CBLAS_INDEX res = 0; 
	complexe_double_t *X = (complexe_double_t *)X_;
  double minValue = abs(X[0].real) + abs(X[0].imaginary);
  
  for (; i < N ; i += incX)
    {
      if (minValue > (abs(X[i].real) + abs(X[i].imaginary))){
        minValue = (abs(X[i].real) + abs(X[i].imaginary));
        res = i;
        }
    }

  return res;
};