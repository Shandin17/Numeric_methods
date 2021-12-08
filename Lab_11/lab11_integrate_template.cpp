#include <stdio.h>
#include <math.h>
#include <time.h>

double a = 0.0;
double b = 0.5*M_PI;
int *TEST;
//number of iterations during intergrate
int num_iter = 10;
//number of integrations, you can change it
int num_integr = 10000;

//integrand
double f (double x)
{
	return (cos(x));
}

//probability distribution
double p (double x)
{
	return (4.0/M_PI*(1.0 - 2.0*x/M_PI));
}

//simple random from a to b
double rand_simple(double a, double b)
{
	double rand_val = a + (b - a) * rand() / (RAND_MAX + 1.0);
	return rand_val;
}

//return value, distributed by p on (a,b)
double rand_distr(double a, double b)
{
	double* x_rand = new double;
	double* y_rand = new double;
	*x_rand = rand_simple(a,b);
	*y_rand = rand_simple(p(b), p(a));
	if( *y_rand <= p(*x_rand))
	{
		double tmp = *y_rand;
		delete x_rand, y_rand;
		return tmp;
	}
	else
	{
		delete x_rand, y_rand;
		return rand_distr(a,b);
	}
}


double Monte_Carlo_simple(double N)
{
	int n_fit = 0;
	int n = 0;
	double tmp_x, tmp_y;
	double S_result = (1.0 - 0.0)*(b - a);
	for(int i = 0; i < N ; i++)
	{
		tmp_x = rand_simple(a,b);
		tmp_y = rand_simple(0.0, 1.0);
		n++;
		if(tmp_y <= f(tmp_x))
			n_fit++;
	}
	S_result = S_result * n_fit/n;
	return S_result;
}

double Monte_Carlo_distr(double N)
{
	int n_fit = 0;
	int n = 0;
	double tmp_x, tmp_y;
	double S_result = (1.0 - 0.0)*(b - a);
	for(int i = 0; i < N ; i++)
	{
		tmp_x = rand_simple(a,b);
		tmp_y = rand_distr(a,b);
		n++;
		if(tmp_y <= f(tmp_x))
			n_fit++;
	}
	S_result = S_result * n_fit/n;
	return S_result;
}



int main (int argc, char **argv)
{
	srand(time(NULL));
	//two copies of the same integral for the two methods
	double I_simple = 0.0, I_distr = 0.0;
	double disp_simple = 0.0, disp_distr = 0.0;
	double aver_simple = 0.0, aver_square_simple = 0.0;
	double aver_distr = 0.0, aver_square_distr = 0.0;
	//num_integr times calculate I_simple and I_distr by 10-iterations Monte-Carlo and find dispersion of the two methods
	for (int j = 0; j < num_integr; j++)
	{
		I_simple = Monte_Carlo_simple(num_iter);
		I_distr = Monte_Carlo_distr(num_iter);
		aver_square_simple = aver_square_simple + pow(I_simple,2);
		aver_simple = aver_simple + I_simple;
		aver_square_distr = aver_square_distr + pow(I_distr,2);
		aver_distr = aver_distr + I_distr;
	}
	aver_square_simple = aver_square_simple/num_integr;
	aver_simple = aver_simple/num_integr;
	aver_square_distr = aver_square_distr/num_integr;
	aver_distr = aver_distr/num_integr;
	disp_simple = aver_square_simple - pow(aver_simple,2);
	disp_distr = aver_square_distr - pow(aver_distr,2);
	printf ("dispersion for simple method = %f dispersion for distribution method = %f", disp_simple, disp_distr);
	return 0;
}
