#include <ctime>
#include "math.h"
#include "stdio.h"

int NP = 10000000;
double c  = 3e10;
double mp = 1.67e-24;
double me = 9.1e-28; 
double q = 4.8e-10;
double* v = new double [NP]; //velocity
double *a = new double [NP]; //acceleration
double *E = new double [NP]; //energy
double *p = new double [NP]; //momentum
double Vt = 1e8;

double fmaxwell (double v)
{
	return (exp(-v*v));
}

double GetRand (double min, double max)
{
	return (min + (max - min)*rand()/(RAND_MAX + 1.0));
}

double RandMaxwell (double Vt)
{
	double v, f;
	for (;;)
	{
		v = GetRand (-5.0, 5.0);
		f = GetRand (0,1);
		if (f < fmaxwell(v))
			return v*Vt;
	}
}

/* 
float Q_rsqrt( float number ) //Quake alg
{	
	const float x2 = number * 0.5F;
	const float threehalfs = 1.5F;

	union {
		float f;
		int i;
	} conv = {number}; // member 'f' set to value of 'number'.
	conv.i = 0x5f3759df - ( conv.i >> 1 );
	conv.f *= threehalfs - x2 * conv.f * conv.f;
	return  conv.f;
}
*/

int main (int argc, char **argv)
{
	clock_t t1 = 0, t2 = 0;
	double worktime = 0.0;
	int cyc = 10;
	double mp2 = mp*0.5, me2 = me*0.5;
	double ap_c = q*mp/c, ae_c = q*me/c;
	
//initialization
	for (int i = 0; i < NP; i++)
	{
		v[i] = RandMaxwell (Vt);
	}	
//not optimized version
	for (int j = 0; j < cyc; j++)
	{
		t1 = t1 + clock();
		for (int i = 0; i < NP; i++)
		{
			if (i % 2 == 0)
			{
				a[i] = q*v[i]*mp/c;
				E[i] = mp*pow(v[i],2)/2;
				p[i] = sqrt(2*mp*E[i]);
			}
			else
			{
				a[i] = -q*v[i]*me/c;
				E[i] = me*pow(v[i],2)/2;
				p[i] = sqrt(2*me*E[i]);
			}
		}
		t2 = t2 + clock();
	}
	worktime = double(t2 - t1)/CLOCKS_PER_SEC;
	printf ("av total work time before opt = %f seconds\n", worktime/cyc);
	t1 = t2 = 0;
// optimized version
	for (int j = 0; j < cyc; j++)
	{
		t1 = t1 + clock();
		if (NP % 2 == 0)
		{
			for (int i = 0; i < NP; i+=2)
			{
				a[i+1] = -ae_c*v[i+1];
				E[i+1] = me2*v[i+1]*v[i+1];
				p[i+1] = me*abs(v[i+1]);

				a[i] = ap_c*v[i];
				E[i] = mp2*v[i]*v[i];
				p[i] = mp*abs(v[i]);
			}
		}
		else if (NP % 2 != 0)
		{
			for (int i = 0; i < NP-1; i+=2)
			{
				a[i+1] = -ae_c*v[i+1];
				E[i+1] = me2*v[i+1]*v[i+1];
				p[i+1] = me*abs(v[i+1]);

				a[i] = ap_c*v[i];
				E[i] = mp2*v[i]*v[i];
				p[i] = mp*abs(v[i]);
			}
			a[NP-1] = ap_c*v[NP-1];
			E[NP-1] = mp2*v[NP-1]*v[NP-1];
			p[NP-1] = mp*abs(v[NP-1]);
		}
		t2 = t2 + clock();
	}
	worktime = double(t2 - t1)/CLOCKS_PER_SEC;
	printf ("av total work time after opt = %f seconds\n", worktime/cyc);
	return 0;
}