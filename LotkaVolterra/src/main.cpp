//
// (c) 2016 Takahiro Hashimoto
//

// 
// Finete difference solution of Lotka-Voltera Equations
//

#include <iostream>
#include <math.h>

const int TimeStepMax = 10000;

const float dt = 0.01;

const float popInit1 = 0.9;
const float popInit2 = 0.9;

const float a = 2./3;
const float b = 4./3;
const float c = 1.;
const float d = 1.;


int main(int argc, char const *argv[])
{

	float t[TimeStepMax];
	float pop1[TimeStepMax], pop2[TimeStepMax];

	// initialize t, pop1, pop2
	t[0] = 0;
	pop1[0] = popInit1;
	pop2[0] = popInit2;

	// update population
	for ( int tIdx = 0; tIdx < TimeStepMax-1; tIdx++)
	{
		t[tIdx + 1] = dt + t[tIdx]; 
		pop1[tIdx+1] = 
			(a*pop1[tIdx] - b*pop1[tIdx]*pop2[tIdx])*dt + pop1[tIdx];
		pop2[tIdx+1] =
			(-c*pop2[tIdx] + d*pop1[tIdx]*pop2[tIdx])*dt + pop2[tIdx];
	}

	float v[TimeStepMax];
	for ( int tIdx = 0; tIdx < TimeStepMax; tIdx++)
	{
		v[tIdx] = -d*pop1[tIdx] + c*log(pop1[tIdx]) -
				   b*pop2[tIdx] + a*log(pop2[tIdx]);
	}

	// save result
	std::cout << "# population evolution simulated by Lotka-Voltera equations" << "\n";
	std::cout << "# t pop1 pop2 v" << "\n";
	for ( int tIdx = 0; tIdx < TimeStepMax; tIdx++)
	{
		std::cout << t[tIdx] << " " 
				  << pop1[tIdx] << " " 
				  << pop2[tIdx] << " " 
				  << v[tIdx] << "\n";
	}

	return 0;
}