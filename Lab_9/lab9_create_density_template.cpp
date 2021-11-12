#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
int NC = 100;	//number of cells
int NP = 10000; //number of particles
double L_cell = 1.0;
double L = L_cell * NC; //box length in normalized units

//density profile
double rho(double x)
{
	double r = 1.0 + 0.5 * sin(2 * M_PI * x / L);
	return r;
}

//random value from xmin to xmax
double GetRand(double xmin, double xmax)
{
	double x = xmin + (xmax - xmin) * rand() / (RAND_MAX + 1.0);
	return x;
}

int main(int argc, char **argv)
{
	double x = 0.0;
	double y = 0.0;
	double w = 0.0;
	int* temp = new int[NC];
	for (int i = 0; i < NC; i++)
		temp[i] = 0;
	
	FILE *output = fopen("Coordinates.dat", "wt");
	fprintf(output, "%s %s\n", "[i]", "[x]");
	for (int i = 0; i < NP; i++)
	{
		do
		{
			x = GetRand(0.0, L);
			y = GetRand(0.0, 1.5);
		} while (y > rho(x));
		fprintf(output, "%i %.1f\n", i + 1, x);
		if(x < L && x > 0.0)
		temp[int(x)]++;
	}
	for (int i = 0; i < NC; i++)
		w = w + temp[i];
	w = NC / w;
	fprintf(output, "\n");
	fprintf(output, "%s %f\n","Weight: ", w);
	fclose(output);
	delete[] temp;
	return 0;
}