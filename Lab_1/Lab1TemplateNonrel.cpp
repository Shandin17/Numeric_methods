#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>

double B0 = 1e-6;//1e-6; //Gs, magnetic field
double E0 = 0.0;//Gs
double v0 = 3e8; //cm/s
double B[3]; //magnetic field
double E[3]; //electric field
double q = 4.8e-10; //CGS units, electric charge of proton/electron
double m = 1.67e-27; //g, proton mass
double qm = q/m;
double v1[3]; //particle velocities for the first pusher
double r1[3]; //particle positions for the first pusher 
double v2[3]; //particle velocities for the second pusher
double r2[3]; //particle positions for the second pusher
double c = 3.0e10; //velocity of light
double dt = 0.001; //seconds, time step


double F(double r[3], double v[3], unsigned i)
{
    switch (i)
        {
        case 0:
            return (  qm * (E[0] + (1.0/c) * (v[1] * B[2] - v[2] * B[1]))  ) ;
        break;
        case 1:
            return (  qm * (E[1] - (1.0/c) * (v[0] * B[2] - v[2] * B[0]) ) ) ;
        break;
        case 2:
            return (  qm * (E[2] + (1.0/c) * (v[0] * B[1] - v[1] * B[0]) ) ) ;
        break;
        }
    return 1;


}


void pusher1 (double dt)
{
    //Euler method

    double Lorenz[3];

    for (unsigned i = 0; i<3; i++ )
    {
       Lorenz[i] = F(r1,v1,i); 
    }

    for (unsigned i = 0; i<3; i++ )
    {
        r1[i] = r1[i] + v1[i] * dt;
       v1[i] = v1[i] + Lorenz[i] * dt;
    }

}

void pusher2 (double dt)
{
    double ar[3];
    double av[3];
    double br[3];
    double bv[3];
    double tmp1[3];
    double tmp2[3];

    for (unsigned i = 0; i<3; i++ )
    {
        ar[i] = av[i] = br[i] = bv[i] = tmp1[i] = tmp2[i] = 0.0;
    }

    for (unsigned i = 0; i<3; i++ )
    {
        ar[i] = dt * v2[i] ;
        av[i] = dt * F(r2,v2,i);
        br[i] = dt * (v2[i] + av[i]) ;
        tmp1[i] = r2[i] + ar[i];
        tmp2[i] = v2[i] + av[i];
    }

    for (unsigned i = 0; i<3; i++ )
    {
      bv[i] = dt * F(tmp1, tmp2, i);
      r2[i] = r2[i] + 0.5 * (ar[i] + br[i] );
      v2[i] = v2[i] + 0.5 * (av[i] + bv[i] );
    }


}

int main (int argc, char ** argv)
{
	double t = 0.0; //current time
	double maxt = 10.0;
	B[2] = B0;  //Bz
	B[0] = B[1] = 0.0;
    E[0] = E[1] = E[2] = 0.0;
	for (int i = 0; i < 3; i++)
	{
        r1[i] = r2[i] = 0.0;
		if (i == 0) 
			v1[i] = v2[i] = v0; //vx
		else
            v1[i] = v2[i] = 0.0;
	}
	char FileName[100];
	sprintf (FileName, "Compare_particle_movers_dt_%1.0le.dat", dt);
	FILE *output = fopen(FileName, "wt");
	for (t = 0.0; t < maxt; t += dt)
	{
		fprintf (output, "%1.4le ", t);
        fprintf (output, "%1.4le %1.4le %1.4le %1.4le %1.4le %1.4le ", r1[0], r1[1], r1[2], v1[0], v1[1], v1[2]);
        fprintf (output, "%1.4le %1.4le %1.4le %1.4le %1.4le %1.4le ", r2[0], r2[1], r2[2], v2[0], v2[1], v2[2]);
		fprintf (output, "\n");
        pusher1(dt);
        pusher2(dt);
	}
	fclose (output);


    return 0;
}
	
