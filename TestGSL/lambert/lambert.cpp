//
// (c) Takahiro Hashimoto
//

#include <iostream>
#include <gsl/gsl_sf_lambert.h>

int
main (void)
{

	const double e = 2.71828182846;

	const int ArrSize = 100;
	const double Start = -1;
	const double Space = 0.1;
	double x[ArrSize], y[ArrSize];

	std::cout << "Lambert W function: " << "\n";

	for (int idx = 0; idx < ArrSize; ++idx)
	{
		x[idx] = Start + Space*idx;
		gsl_sf_result result;
		int status = gsl_sf_lambert_W0_e(x[idx], &result);
		if ( status != 0 )
		{
			y[idx] = -1.;
		} 
		else
		{
			y[idx] = result.val;
		}
		std::cout << x[idx] << " " << y[idx] << "\n";
	}

	return 0;
}