//
//
// (c) Takahiro Hashimoto
//
//

#include <iostream>
#include <gsl/gsl_sf_bessel.h>

int
main (void)
{
  double x = 5.0;
  double y = gsl_sf_bessel_J0 (x);
  std::cout << "J0(" << x << ") = " << y << "\n";
  return 0;
}