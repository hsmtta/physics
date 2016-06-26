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
#include <GLUT/glut.h>

using namespace std;

const int TimeBuff = 200;
const double dt = 0.1;
const string ResultFilename = "result.dat";

double a, b, c, d;
double v0;

bool isUseInvariant = false;

double tArr[TimeBuff];
double pop1[TimeBuff], pop2[TimeBuff];

double *xa;

void PrintUsageAndExit(const string& arg0);
double Wm1Extended(const double x, bool& isInsideRange);
double W0Extended(const double x, bool& isInsideRange);
double getInvariant(const double p1, const double p2);
void update(const double t, const double p1, const double p2,  
	        double &tn, double &pn1, double &pn2);

void idle(void);
void timer(int value);
void display(void);
void init(void);
void deinit(void);

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

	tArr[0] = 0;
	pop1[0] = popInit1;
	pop2[0] = popInit2;

	if ( isUseGlut )
	{
		// init glut
		glutInit(&argc, argv);
		glutInitWindowSize(300,300);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA );
		glutCreateWindow("Evolution of population 1 and 2");
		glutDisplayFunc(display);
		// glutIdleFunc(idle);
		// glEnable(GL_BLEND);
		// glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		// glEnable(GL_LINE_SMOOTH);
		// glHint(GL_LINE_SMOOTH, GL_NICEST);

		// glEnable(GL_POINT_SMOOTH);
		// glHint(GL_POINT_SMOOTH, GL_NICEST);

		glutTimerFunc(50, timer, 0);

		init();
		glutMainLoop();
		deinit();
	}
	else
	{
		// simulate and save result
		double v[TimeBuff];
		v[0] = v0;

		for ( int tIdx = 0; tIdx < TimeBuff-1; tIdx++)
		{
			// update population
			const double t = tArr[tIdx];
			const double p1 = pop1[tIdx];
			const double p2 = pop2[tIdx];
			double tn, pn1, pn2;
			update(t, p1, p2, tn, pn1, pn2);
			tArr[tIdx+1] = tn;
			pop1[tIdx+1] = pn1;
			pop2[tIdx+1] = pn2;

			v[tIdx+1] = getInvariant(pn1, pn2);

		}

		// save result
		cout << "save result as " << ResultFilename << "\n";

		ofstream ofs;
		ofs.open(ResultFilename.c_str());

		ofs << "# population evolution simulated by Lotka-Volterra equations" << "\n";
		ofs << "# t pop1 pop2 v" << "\n";
		for ( int tIdx = 0; tIdx < TimeBuff; tIdx++)
		{
			ofs << tArr[tIdx] << " "
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
	cout << "-d                   display result\n";

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
	pn2 = (-c*p2 + d*p1*p2)*dt + p2;

	if ( isUseInvariant )
	{
		double pn2_;
		bool isInsideRange;

		if ( b/a*p2 > 1 )
		{
			pn2_ = -a/b*Wm1Extended(
				-b/a*(exp(1./a*(-c*log(pn1)+d*pn1-v0))), 
				isInsideRange );
		}
		else
		{
			pn2_ = -a/b*W0Extended(
				-b/a*(exp(1./a*(-c*log(pn1)+d*pn1-v0))), 
				isInsideRange );
		}

		if ( isInsideRange ) pn2 = pn2_;
	}
}

void idle(void)
{
  glutPostRedisplay();
}

void timer(int value) 
{
	glutPostRedisplay();
	glutTimerFunc(16 , timer , 0);
}

void display(void)
{
	static int count = 0;

	const double t = tArr[count];
	const double p1 = pop1[count];
	const double p2 = pop2[count];

	double tn, pn1, pn2;
	update(t, p1, p2, tn, pn1, pn2);

	if ( ++count == TimeBuff) count = 0;

	tArr[count] = tn;
	pop1[count] = pn1;
	pop2[count] = pn2;

	glClear(GL_COLOR_BUFFER_BIT);

	glLineWidth(1.0);
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINE_STRIP);
	for (int i=count+1, j=0; j<TimeBuff; ++i, ++j)
	{
		if ( i > TimeBuff -1) i = 0;
		glVertex2d(xa[j], pop1[i]/7.);
	}
	glEnd();

	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINE_STRIP);
	for (int i=count+1, j=0; j<TimeBuff; ++i, ++j)
	{
		if ( i > TimeBuff -1) i = 0;
		glVertex2d(xa[j], pop2[i]/7);
	}
	glEnd();

	// exit(1);

	// glVertex2d(0,0);
	// glVertex2d(0,0.5);
	// glVertex2d(0.5,0.5);
	// glVertex2d(0.5,0);

	// glEnd();

	glFlush();
	glutSwapBuffers();
}

void init(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	
	const double xStart = -1;
	const double xEnd = 1;
	const double delta = (xEnd - xStart)/(TimeBuff -1);
	xa = new double[TimeBuff];

	for (int i = 0; i < TimeBuff; ++i)
	{
		xa[i] = xStart + i*delta;
	}
}

void deinit(void)
{
	delete[] xa;
}





