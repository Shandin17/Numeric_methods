#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>




typedef struct tagPARTICLE
{
	double r[3]; //coordinate
	double v[3]; //velocity
	struct tagPARTICLE *p;
}PARTICLE;

typedef struct tagCELL  
{
	int r[3]; //coordinate
	double E[3]; //electric field
	double B[3]; //magnetic field
	double rho_c; //charge density
	double j[3]; //current density
	struct tagCELL* neighbours[3][3][3]; //neighbours
	PARTICLE *p; //particles
}CELL;


int main (int argc, char **argv)
{
	const double c = 3.0e10;
	const double me = 0.911e-27;
	const double q = 4.8e-10;
	const double pi = 3.141592653589793;
	const double n = 1.0;


	double le = c*sqrt(me)/sqrt( 4.0*pi * n * pow(q,2) );
	double R = 1.0e18;
	double a = 1.0;
	double N_cell = pow(R,3) / pow(a*le,3); // amount of macroparticles
	double N_p = 10 * N_cell; //full amount of particles
	double Comp = ( N_cell*sizeof(CELL) + N_p*sizeof(PARTICLE) )/(pow(1024,3));
	double Cluster_mem = 158976.0 * 32.0;
	std::cout << Comp / Cluster_mem << " - amount of clusters required" << std::endl;

	std::cout << sizeof(CELL) << std::endl;

}