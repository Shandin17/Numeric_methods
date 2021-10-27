#include <stdlib.h>
#include <cstdio>
#include <math.h>
#include <fstream>
#include <iostream>
using namespace std;

int Nt = 100;  //number of time steps
int Nx = 100; //number of points
double xmin = -1.0; //large intervals are taken in order to solution stay in modelig box
double xmax = 3.0;
double tmax = 1.0;


double q0 (double x)
{
	if ((x < 0) || (x >= 1))
		return 0.0;
	return (sin(M_PI*x));
}
//flux
double f(double qu)
{
	return qu*qu*0.5;
}

double* Upwind(double c1, double* Q, double* dQ)
{
	int i;
	int n;
	for(n = 1; n <= Nt; n++)
	{
		for(i = 1; i < Nx-1; i++)
		{
			if(tan(f(Q[i])) > 0.0)
			dQ[i] = - c1*(f(Q[i]) - f(Q[i-1]));
			else if(tan(f(Q[i])) < 0.0)
			dQ[i] = - c1*(f(Q[i+1]) - f(Q[i]));
		}
		for(i = 1; i < Nx-1; i++)
			Q[i] = Q[i] + dQ[i];
	}
	return Q;
}

double* Leapfrog(double c1, double* Qodd, double* Qeven,
				double* dQodd, double* dQeven)
{
	unsigned i;
	for (i = 1; i < Nx - 1; i++)
	{
		if (tan(f(Qeven[i])) > 0.0)
			dQeven[i] = -c1 * (f(Qeven[i]) - f(Qeven[i - 1]));
		else if (tan(f(Qeven[i])) < 0.0)
			dQeven[i] = -c1 * (f(Qeven[i + 1]) - f(Qeven[i]));
	}
	for (i = 1; i < Nx - 1; i++)
		Qeven[i] = Qeven[i] + dQeven[i];

	for(int n = 3; n <= Nt; n++)
	{
		if(n%2 == 0)
		{
			for(i = 1; i < Nx-1; i++)
			dQeven[i] = - c1*(f(Qodd[i+1]) - f(Qodd[i-1]));
			for(i = 1; i < Nx-1; i++)
			Qeven[i] = Qeven[i] + dQeven[i]; 
		}
		else
		{
			for(i = 1; i < Nx-1; i++)
			dQodd[i] = - c1*(f(Qeven[i+1]) - f(Qeven[i-1]));
			for(i = 1; i < Nx-1; i++)
			Qodd[i] = Qodd[i] + dQodd[i]; 
		}
	}
	if(Nt%2 == 0)
	return Qeven;
	else
	return Qodd;
}

int main (int argc, char **argv)
{
	double** q = new double *[4]; //4 arrays: 1 for Upwind; 2 for odd and even
	// steps of leapfrog; 1 for lax-wendroff
	double** dq = new double *[4]; //4 arrays for incremets at each step
	double t, c1;
	double *x = new double[Nx];
	for (int i = 0; i < 4; i++)
	{
		q[i] = new double [Nx];
		dq[i] = new double [Nx];
	}
	double dx = (xmax - xmin) / (Nx - 1);  
	double dt = tmax / (Nt - 1);
	c1 = dt/dx;
	for (int i = 0; i < Nx; i++) //initial conditions
	{
		x[i] = xmin + i*dx;
		q[0][i] = q[1][i] = q[2][i] = q[3][i] = q0(x[i]);
		dq[0][i] = dq[1][i] = dq[2][i] = dq[3][i] = 0.0;
	}
	
	double* temp;
	temp = q[0];
	q[0] = Upwind(c1, q[0], dq[0]);
	delete [] temp;
	temp = q[2];
	q[2] = Leapfrog(c1, q[1], q[2], dq[1], dq[2]);
	delete [] temp;
	
	FILE *output = fopen ("q.dat", "wt");
	for (int i = 0; i < Nx; i++)
		fprintf (output, "%f %f %f %f\n", x[i], q[0][i], q[2][i], q[3][i]);
	fclose (output);
	cout << 5%2 << endl;
	return 0;
}