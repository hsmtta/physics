//
// (c) 2016 Takahiro Hashimoto
//

//
// Finite difference solution of Lotka-Volterra Equations
//

#include <iostream>
#include <fstream>
#include <math.h>
#include <gsl/gsl_sf_lambert.h>
#include <GL/glut.h>

using namespace std;

const int TimeBuff = 10000;
const double dt = 0.01;
const string ResultFilename = "result.dat";

double a, b, c, d;
double v0;

bool isUseInvariant = false;

double time[TimeBuff];
double pop1[TimeBuff], pop2[TimeBuff];

void PrintUsageAndExit(const string& arg0);
double Wm1Extended(const double x, bool& isInsideRange);
double W0Extended(const double x, bool& isInsideRange);
double getInvariant(const double p1, const double p2);
void update(const double t, const double p1, const double p2,  
	        double &tn, double &pn1, double &pn2);

void display(void);
void initDisplay(void);

int main(int argc, char *argv[])
{
	double popInit1 = 1.5;
	double popInit2 = 1.0;

	bool isUseGlut = false;

	a = 2./3;
	b = 4./3;
	c = 1.;
	d = 1.;

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
		else if ( string(argv[argIdx]) == "-d")
		{
			isUseGlut = true;
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


	v0 = getInvariant(popInit1, popInit2);

	time[0] = 0;
	pop1[0] = popInit1;
	pop2[0] = popInit2;

	if ( isUseGlut )
	{
		// init glut
		glutInit(&argc, argv);
		glutCreateWindow(argv[0]);
		glutDisplayFunc(display);

		initDisplay();

		glutMainLoop();
	}
	else
	{
		// simulate and save result
		double v[TimeBuff];
		for ( int tIdx = 0; tIdx < TimeBuff-1; tIdx++)
		{
			// update population
			const double t = time[tIdx];
			const double p1 = pop1[tIdx];
			const double p2 = pop2[tIdx];
			double tn, pn1, pn2;
			update(t, p1, p2, tn, pn1, pn2);
			time[tIdx+1] = tn;
			pop1[tIdx+1] = pn1;
			pop2[tIdx+1] = pn2;

			v[tIdx+1] = getInvariant(p1, p2);

		}

		// save result
		cout << "save result as " << ResultFilename << "\n";

		ofstream ofs;
		ofs.open(ResultFilename.c_str());

		ofs << "# population evolution simulated by Lotka-Volterra equations" << "\n";
		ofs << "# t pop1 pop2 v" << "\n";
		for ( int tIdx = 0; tIdx < TimeBuff; tIdx++)
		{
			ofs << time[tIdx] << " "
				<< pop1[tIdx] << " "
				<< pop2[tIdx] << " "
				<< v[tIdx] << "\n";
		}

		ofs.close();
	}

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

double getInvariant(const double p1, const double p2)
{
	return -a*log(p2) + b*p2 -c*log(p1) + d*p1;
}

void update(const double t, const double p1, const double p2, 
	        double &tn, double &pn1, double &pn2)
{
	tn = dt + t;
	pn1 = (a*p1 - b*p1*p2)*dt + p1;
	pn2 = (-c*pn2 + d*p1*p2)*dt + p2;

	if ( isUseInvariant )
	{
		double pn2_;
		bool isInsideRange;

		if ( b/a*p2 > 1 )
		{
			pn2_ = -a/b*Wm1Extended(
				-b/a*(exp(1./a*(-c*log(p1)+d*p1-v0))), 
				isInsideRange );
		}
		else
		{
			pn2_ = -a/b*W0Extended(
				-b/a*(exp(1./a*(-c*log(p1)+d*p1-v0))), 
				isInsideRange );
		}

		if ( isInsideRange ) pn2 = pn2_;
	}
}

void display(void)
{
	static int count = 0;

	const double t = time[count];
	const double p1 = pop1[count];
	const double p2 = pop2[count];

	double tn, pn1, pn2;
	update(p1, p2, t, pn1, pn2, tn);

	if ( ++count == TimeBuff -1) count = 0;

	time[count] = tn;
	pop1[count] = pn1;
	pop2[count] = pn2;

	glClear(GL_COLOR_BUFFER_BIT);

	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(15, 0, 0);
	glEnd();

	glFlush();
}

void initDisplay(void)
{
	glClearColor(0.0, 0.0, 1.0, 1.0);
}