//
// (c) 2016 Takahiro Hashimoto
//

//
// Finite difference solution of Lotka-Volterra Equations
//

#include <iostream>
#include <math.h>
#include <gsl/gsl_sf_lambert.h>

const int TimeStepMax = 10000;

const double dt = 0.01;

const double popInit1 = 0.9;
const double popInit2 = 0.9;

const double a = 2./3;
const double b = 4./3;
const double c = 1.;
const double d = 1.;

const string ResultFilename = "result.dat";

void PrintUsageAndExit(const string& arg0);
double W0Extended(const double x);

int main(int argc, char const *argv[])
{

	bool isUseConservedQuantity = false;

	// process command line arguments
	for (int argIdx = 1; argIdx < argc; argIdx++)
	{
		if ( string(argv[argIdx]) == "--help" || string(argv[argIdx]) == "-h")
		{
			PrintUsageAndExit(argv[0]);
		}
		else if ( string(argv[argIdx]) == "-c" )
		{
			cout << "simulate with conserved quantity" << "\n";
			isUseConservedQuantity = true;
		}
		else
		{
			PrintUsageAndExit(argv[0]);
		}
	}

	double t[TimeStepMax];
	double pop1[TimeStepMax], pop2[TimeStepMax];

	// initialize t, pop1, pop2
	t[0] = 0;
	pop1[0] = popInit1;
	pop2[0] = popInit2;

	double v0 = -d*popInit1 + c*log(popInit1)
				-b*popInit2 + a*log(popInit2);

	// update population
	for ( int tIdx = 0; tIdx < TimeStepMax-1; tIdx++)
	{
		t[tIdx + 1] = dt + t[tIdx];
		pop1[tIdx+1] =
			(a*pop1[tIdx] - b*pop1[tIdx]*pop2[tIdx])*dt + pop1[tIdx];

		if ( isUseConservedQuantity )
		{
			const double p1 = pop1[tIdx+1];
			double p2 = W0Extended(-b/a*exp((v0-d*log(p1)+c*p1)/a));
			p2 = -a/b*p2;
			pop2[tIdx+1] = p2;
		} 
		else
		{
			pop2[tIdx+1] =
				(-c*pop2[tIdx] + d*pop1[tIdx]*pop2[tIdx])*dt + pop2[tIdx];	
		}
	}

	double v[TimeStepMax];
	for ( int tIdx = 0; tIdx < TimeStepMax; tIdx++)
	{
		v[tIdx] = -d*pop1[tIdx] + c*log(pop1[tIdx])
				  -b*pop2[tIdx] + a*log(pop2[tIdx]);
	}

	// save result

	cout << "save result as " << ResultFilename << "\n";

	cout << "# population evolution simulated by Lotka-Volterra equations" << "\n";
	cout << "# t pop1 pop2 v" << "\n";
	for ( int tIdx = 0; tIdx < TimeStepMax; tIdx++)
	{
		std::cout << t[tIdx] << " "
				  << pop1[tIdx] << " "
				  << pop2[tIdx] << " "
				  << v[tIdx] << "\n";
	}

	return 0;
}

void PrintUsageAndExit(const string& arg0)
{
	cout << "Usage: " << arg0 << "[options]" << "\n";
	cout << "Options: -h | --help       print this usage" << "\n";
	cout << "Options: -c                use conserved quantity" << "\n";
}

double W0Extended(const double x)
{

	gsl_sf_result result;
	int status = gsl_sf_lambert_W0_e(x, &result);

	double y;
	if ( status != 0 )
	{
		y = -1.;
	} 
	else
	{
		y = result.val;
	}

	return y;
}