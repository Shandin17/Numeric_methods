#include "stdio.h"
#include "math.h"
#include <iostream>
#include <fstream>

using namespace std;

#define NP 1024 //number of particles
#define NC 8 //number of cells in each direction

typedef struct tagPART
{
    double x, y; //coordinates
    double vx, vy; //velocity
    double vxi, vyi; //initial velocity
    double phi; //angle
    double vpar; //velocity component along initial
    double vperp; //velocity component perpendicular to initial
    double h; //0.5*(v(t)^2 - v(0)^2)
} PART;

 
double Ex[NC][NC];  //field components
double Ey[NC][NC]; 
PART electrons[NP];  // 
double H = 1.0; //cell size in lambda_d
double k = 0.5;  //timestep in wpe^-1
double Vte = k/H; //initial electrons thermal speed, it is equal to lambda_d*wpe, and in normilized units is k/H
double L = double(NC);  //area size in cells
double A = 0.1; //amplitude of electric field fluctuations, in normilized units of k^2*e/(m*H)
//Note! If we normalize everything to dt and dx, then there is no dx and dt in discrete equations!

//random value from min to max
double GetRand(double min, double max)
{
	double rand_val = min + (max - min) * rand() / (RAND_MAX + 1.0);
	return rand_val;
}

//Maxwell distributin
double fmaxwell(double v)
{
	double f = (1/(M_PI*pow(Vte,2))) * exp(-pow(v,2)); // normilized function (2d maxwell)
	return f;										   // v = v/Vte 
}

//random value distributed by Maxwell, any method you prefer
double RandMaxwell(double Vt)
{
	double* v0 = new double;
	double* f0 = new double;
	*v0 = GetRand(-Vte*5.0, Vte*5.0);
	*f0 = GetRand(0.0,1.0);
	if (*f0 <= fmaxwell(*v0))
	{
		return *v0;
	}
	else
	{
		delete v0, f0;
		return RandMaxwell(Vte);
	}
}

//initialize particles
void Init()
{
	int i, j;
	double v;
	for (i = 0; i < NP; i++)
	{
		electrons[i].x = GetRand(0,L);
		electrons[i].y = GetRand(0,L);
		electrons[i].vx = RandMaxwell(Vte);
		electrons[i].vy = RandMaxwell(Vte);
		electrons[i].vxi = electrons[i].vx;
		electrons[i].vyi = electrons[i].vy;
		electrons[i].phi = 0.0;
		v = sqrt (electrons[i].vx*electrons[i].vx +
		  + electrons[i].vy*electrons[i].vy);
		electrons[i].vpar = v;
		electrons[i].vperp = 0.0;
	}
}

void MoveParticles()
{
	double v, vi;
	int temp_x, temp_y;
	for (int i = 0; i < NP; i++)
	{
		temp_x = int(electrons[i].x) % NC;
		temp_y = int(electrons[i].y) % NC;
		if (temp_x < 0)
			temp_x = temp_x + NC;
		if (temp_y < 0)
			temp_y = temp_y + NC;
		electrons[i].x = electrons[i].x + electrons[i].vx;
		electrons[i].y = electrons[i].y + electrons[i].vy;
		electrons[i].vx = electrons[i].vx + Ex[temp_x][temp_y];
		electrons[i].vy = electrons[i].vy + Ey[temp_x][temp_y];		
	}
}

void WriteOutput(double time, FILE *fout)
{
	double h_av, vi, v;
	for (int i = 0; i < NP; i++)
	{
		v = sqrt(electrons[i].vx * electrons[i].vx +
				+electrons[i].vy * electrons[i].vy);
		vi = sqrt(electrons[i].vxi * electrons[i].vxi +
				+electrons[i].vyi * electrons[i].vyi);
		electrons[i].h = 0.5 * (pow(v, 2) - pow(vi, 2));
		h_av = h_av + electrons[i].h/NP; 
	}
	fprintf(fout, "%f %f\n", time, h_av);
}

int main(int argc, char **argv)
{
	int step;
	int T = 1000;
	Init();
	FILE *fout = fopen("output.dat", "wt");
	for (step = 0; step < T; step++)
	{
		//put random fields into cells
		for (int i = 0; i < NC; i++)
			for (int j = 0; j < NC; j++)
			{
				Ex[i][j] = GetRand(-A, A);
				Ey[i][j] = GetRand(-A, A);
			}
		MoveParticles();
		WriteOutput(step, fout);
	}
	fclose(fout);
	return 0;
}