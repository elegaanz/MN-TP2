#include "mnblas.h"
#include "complexe.h"

void mncblas_sswap(const int N, float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX, 
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_cswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
  complexe_float_t *PX = (complexe_float_t *)X;
  complexe_float_t *PY = (complexe_float_t *)Y;

  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register complexe_float_t  save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = PY [j] ;
      PY [j] = PX [i] ;
      PX [i] = save ;
    }

  return ;
}

void mncblas_zswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
  complexe_double_t *PX = (complexe_double_t *)X;
  complexe_double_t *PY = (complexe_double_t *)Y;

  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register complexe_double_t  save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = PY [j] ;
      PY [j] = PX [i] ;
      PX [i] = save ;
    }

  return ;
}

