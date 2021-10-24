#include <stdlib.h>
#include <cstdio>
#include <math.h>
#include <fstream>
#include <iostream>
using namespace std;

int Nc = 30; //number of characteristics
int Nx = 10; //number of points along each characteristic
double xmin = 0;
double xmax = 2.0;
double x0max = 1.0;

double q0 (double x)
{
	return (sin(M_PI*x));
}

int main (int argc, char **argv)
{
	double x0, x, t;
	double dx = (xmax - xmin) / (Nx - 1);  //interval between characteristics points
	double dx0 = (x0max - xmin) / (Nc - 1); //interval between beginnings of characterisics
	FILE *output = fopen ("Characteristics.dat", "wt");
	FILE *output1 = fopen ("X.dat", "wt");
	//got Nc rows and Nx cols

	//fprintf (output, "%i %le %le %le %le\n", dx0, dx );
	x0 = 0;
	for (int j = 0; j <Nx; j++)
	{
		fprintf(output, "%le\t", 0.0 ); // 0.0 point case
		fprintf(output1, "%le\t", 0.0 );
	}
    fprintf(output, "\n");
	fprintf(output1, "\n");

	x0 = x0 + dx0;
	for (int i = 0; i <Nc-2; i++)
	{
		x = x0;
		for (int j = 0; j <Nx; j++)
		{
			t = (x - x0)/q0(x0);
			fprintf(output, "%le\t", t );
			fprintf(output1, "%le\t", x );
			x = x + dx;
		}
		x0 = x0 + dx0;
		fprintf(output, "\n" );
		fprintf(output1, "\n" );
	}

	for (int j = 0; j <Nx; j++)
	{
		fprintf(output, "%le\t", 0.0 ); // 1.0 point case
		fprintf(output1, "%le\t", 0.0 );
	}
	


	fclose(output);
	fclose(output1);
}