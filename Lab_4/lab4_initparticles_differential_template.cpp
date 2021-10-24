#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

double m = 1.67e-24; //g
double T = 1.0e5; //K
double k = 1.38e-16; //boltsmann
double vth = sqrt(2.0*k*T/m);
double vmin = -5.0*vth;
double vmax = 5.0*vth;
double Pi = 3.141592653;
int N = 1000; //number of particles

//differential Maxwell distribution
//we do not need norm here!
double fmaxwell(double v)
{
	double f = exp(-m*pow(v,2)/(2*k*T));
	return f;
}

//return random value from 0 to 1
double GetRand()
{
	double x = rand() / (RAND_MAX + 1.0);
	return (x);
}

//return value, distributed by Gaussian
double RandMaxwell()
{
	double* v0 = new double;
	double* f0 = new double;
	*v0 = vmin + GetRand() * (vmax - vmin);
	*f0 = GetRand();
	if (*f0 <= fmaxwell(*v0))
	{
		return *v0;
	}
	else
	{
		delete v0, f0;
		return RandMaxwell();
	}
}

int main(int argc, char **argv)
{
	//std::cout << vth << std::endl;
	FILE *output = fopen("Particles_diff.dat", "wt");
	double f, vx, vy, vz;
	for (int i = 0; i < N; i++)
	{
		vx = RandMaxwell();
		vy = RandMaxwell();
		vz = RandMaxwell();
		fprintf (output, "%i %le %le %le %le\n", i, vx, vy, vz, sqrt(vx*vx + vy*vy + vz*vz));
	}
	fclose (output);
}