#include "complexe.h"

complexe_float_t add_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;
  
  return r ;
}

complexe_double_t add_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;
  
  return r ;
}

complexe_float_t mult_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  r.real = (c1.real * c2.real) - (c1.imaginary * c2.imaginary);
  r.imaginary = (c1.real * c2.imaginary) + (c1.imaginary * c2.real) ;
  
  return r ;
}

complexe_double_t mult_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
  {
  complexe_double_t r ;

  r.real = (c1.real * c2.real) - (c1.imaginary * c2.imaginary);
  r.imaginary = (c1.real * c2.imaginary) + (c1.imaginary * c2.real) ;
  
  return r ;
}
  

complexe_float_t div_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r, tmp1, tmp2, conjugc2 ;

  if(c2.imaginary == 0){
    r.real = c1.real/c2.real;
    r.imaginary = c1.imaginary/c2.imaginary;
    return r;
  }

  conjugc2.real = c2.real;
  conjugc2.imaginary = -c2.imaginary;

  tmp1 = mult_complexe_float(c1, conjugc2) ;
  tmp2 = mult_complexe_float(c2, conjugc2) ;
  
  r.real = tmp1.real / tmp2.real;
  r.imaginary = tmp1.imaginary / tmp2.real;

  return r ;
}

complexe_double_t div_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r, tmp1, tmp2, conjugc2 ;

  if(c2.imaginary == 0){
    r.real = c1.real/c2.real;
    r.imaginary = c1.imaginary/c2.imaginary;
    return r;
  }

  conjugc2.real = c2.real;
  conjugc2.imaginary = -c2.imaginary;

  tmp1 = mult_complexe_double(c1, conjugc2) ;
  tmp2 = mult_complexe_double(c2, conjugc2) ;
  
  r.real = tmp1.real / tmp2.real;
  r.imaginary = tmp1.imaginary / tmp2.real;
  
  return r ;
}
