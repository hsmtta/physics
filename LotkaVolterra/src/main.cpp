//
// (c) 2016 Takahiro Hashimoto
//

//
// Finite difference solution of Lotka-Volterra Equations
//

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_sf_lambert.h>

using namespace std;

const int TimeStepMax = 10000;
const double dt = 0.01;
const string ResultFilename = "result.dat";

void PrintUsageAndExit(const string& arg0);
double Wm1Extended(const double x, bool& isInsideRange);
double W0Extended(const double x, bool& isInsideRange);

int main(int argc, char const *argv[])
{
	double popInit1 = 1.5;
	double popInit2 = 1.0;

	double a = 2./3;
	double b = 4./3;
	double c = 1.;
	double d = 1.;

	bool isUseInvariant = false;

	// process command line arguments
	for (int argIdx = 1; argIdx < argc; argIdx++)
	{
		if ( string(argv[argIdx]) == "--help" || string(argv[argIdx]) == "-h")
		{
			PrintUsageAndExit(argv[0]);
		}
		else if ( string(argv[argIdx]) == "-i" )
		{
			cout << "simulate with invariant" << "\n";
			isUseInvariant = true;
		}
		else if ( string(argv[argIdx]) == "-p" )
		{
			cout << "set parameters" << "\n";
			if ( argIdx + 4 < argc )
			{
				a = atof(argv[argIdx+1]);
				b = atof(argv[argIdx+2]);
				c = atof(argv[argIdx+3]);
				d = atof(argv[argIdx+4]);
				argIdx = argIdx + 4;
			} 
			else
			{
				PrintUsageAndExit(argv[0]);
			}
		}
		else if ( string(argv[argIdx]) == "-n" ) 
		{
			if ( argIdx + 2 < argc )
			{
				popInit1 = atof(argv[argIdx+1]);
				popInit2 = atof(argv[argIdx+2]);
				argIdx = argIdx + 2;
			}
			else
			{
				PrintUsageAndExit(argv[0]);
			}
		}
		else
		{
			PrintUsageAndExit(argv[0]);
		}
	}

	cout << "simulation parameters: a = " << a 
		 << ", b = " << b
		 << ", c = " << c
		 << ", d = " << d << "\n";
	cout << "initial populations: pop1 = " << popInit1
		 << ", pop2 = " << popInit2 << "\n";

	// init glut

	double t[TimeStepMax];
	double pop1[TimeStepMax], pop2[TimeStepMax];

	// initialize t, pop1, pop2
	t[0] = 0;
	pop1[0] = popInit1;
	pop2[0] = popInit2;

	double v0 = -a*log(popInit2) + b*popInit2
				-c*log(popInit1) + d*popInit1;

	// update population
	for ( int tIdx = 0; tIdx < TimeStepMax-1; tIdx++)
	{
		t[tIdx + 1] = dt + t[tIdx];
		pop1[tIdx+1] =
			(a*pop1[tIdx] - b*pop1[tIdx]*pop2[tIdx])*dt + pop1[tIdx];

		pop2[tIdx+1] =
			(-c*pop2[tIdx] + d*pop1[tIdx]*pop2[tIdx])*dt + pop2[tIdx];

		if ( isUseInvariant )
		{
			const double p1 = pop1[tIdx+1];
			const double p2 = pop2[tIdx+1];
			double p2_;
			bool isInsideRange;

			if ( b/a*p2 > 1 )
			{
				//cout << "-1 branch" << "\n";
				p2_ = -a/b*Wm1Extended(
					-b/a*(exp(1./a*(-c*log(p1)+d*p1-v0))), 
					isInsideRange );
			}
			else
			{
				//cout << "0 branch" << "\n";
				p2_ = -a/b*W0Extended(
					-b/a*(exp(1./a*(-c*log(p1)+d*p1-v0))), 
					isInsideRange );
			}

			if ( isInsideRange ) pop2[tIdx+1] = p2_;
		}
	}

	double v[TimeStepMax];
	for ( int tIdx = 0; tIdx < TimeStepMax; tIdx++)
	{
		v[tIdx] = -a*log(pop2[tIdx]) + b*pop2[tIdx]
				  -c*log(pop1[tIdx]) + d*pop1[tIdx];
	}

	// save result
	cout << "save result as " << ResultFilename << "\n";

	ofstream ofs;
	ofs.open(ResultFilename.c_str());

	ofs << "# population evolution simulated by Lotka-Volterra equations" << "\n";
	ofs << "# t pop1 pop2 v" << "\n";
	for ( int tIdx = 0; tIdx < TimeStepMax; tIdx++)
	{
		ofs << t[tIdx] << " "
			<< pop1[tIdx] << " "
			<< pop2[tIdx] << " "
			<< v[tIdx] << "\n";
	}

	ofs.close();

	return 0;
}

void PrintUsageAndExit(const string& arg0)
{

	cout << "Usage: " << arg0 << " <options>" << "\n";
	cout << "Options:" << "\n";
	cout << "-h | --help          print this usage" << "\n";
	cout << "-i                   use invariant" << "\n";
	cout << "-p <a> <b> <c> <d>   set parameters a, b, c and d" << "\n";
	cout << "-n <p1> <p2>         set initial population p1, p2" << "\n";

	exit(0);
}

double Wm1Extended(const double x, bool& isInsideRange)
{

	gsl_sf_result result;
	int status = gsl_sf_lambert_Wm1_e(x, &result);

	double y;
	if ( status != 0 )
	{
		y = -1.;
		isInsideRange = false;
	} 
	else
	{
		y = result.val;
		isInsideRange = true;
	}

	return y;
}

double W0Extended(const double x, bool& isInsideRange)
{

	gsl_sf_result result;
	int status = gsl_sf_lambert_W0_e(x, &result);

	double y;
	if ( status != 0 )
	{
		y = -1.;
		isInsideRange = false;
	} 
	else
	{
		y = result.val;
		isInsideRange = true;
	}

	return y;
}






