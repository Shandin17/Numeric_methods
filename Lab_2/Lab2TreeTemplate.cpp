#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

typedef struct tagCELL
{
	double x_c;  //coordinate of centre
	double x_m; //coordinate of mass center
    double x_p;
	double size; //size
	double M; //summary mass

	struct tagCELL* daughters[2]; //left and right
} CELL;

int np = 3; //number of particles
double L = 1.0; //size of root cell
double x[3] = {0.7, 0.27, 0.3}; //coordinates of particles
double m = 1.0; //mass of 1 particle

//create left and right daughter cells
void DivideCell (CELL *curc)
{
	CELL* left = new CELL;
	CELL* right = new CELL;
	
	curc->daughters[0] = left;
	curc->daughters[0]->size = curc->size * 0.5;
	curc->daughters[0]->M = 0.0;
	curc->daughters[0]->x_m= 0.0;
    curc->daughters[0]->x_p= 0.0;
	curc->daughters[0]->x_c = curc->x_c - 0.25*curc->size;
	
	curc->daughters[1] = right;
	curc->daughters[1]->size = curc->size * 0.5;
	curc->daughters[1]->M = 0.0;
	curc->daughters[1]->x_m= 0.0;
    curc->daughters[1]->x_p= 0.0;
	curc->daughters[1]->x_c = curc->x_c + 0.25*curc->size;
	


	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			curc->daughters[i]->daughters[j] = NULL;
}

//add particle with coordinate x to the tree
void UpdateTree (double x, CELL *c)
{   

    if (c->M == 0.0 )
    {
        c->M = m;
        c->x_m = x;
    }
    else
    { 
        c->x_m =  (c->x_m*c->M + x*m)/(c->M+m);
        c->M = c->M + m;
        if(x < c->x_c)
        {
            if (c->daughters[0] == NULL)
                {
                DivideCell(c);
                }
              UpdateTree(x,c->daughters[0]);
        }
        if(x >= c->x_c)
        {
            if (c->daughters[1] == NULL)
            {
             DivideCell(c);
             }
            UpdateTree(x,c->daughters[1]);
        }
        np = np-1;
    }
}

//walk tree and print cells parameters
//COOL METHOD OF TREE MOVING
void WalkTree (CELL *curc, FILE* f)
{
	fprintf (f, "Cell size %f x_c %f x_m %f M %f\n", curc->size, curc->x_c, curc->x_m, curc->M);
	if (curc->daughters[0] != NULL)
		WalkTree (curc->daughters[0], f);
	if (curc->daughters[1] != NULL)
		WalkTree (curc->daughters[1], f);
}

int main ( int argc, char *argv[] )
{
    if (m < 0)
    return -1;

	int i;
	CELL *root = new CELL;
	root->size = L;
	root->M = 0.0;
	root->x_c = 0.5*L;
	root->x_m = 0.0;
	root->daughters[0] = root->daughters[1] = NULL;
	//Build tree
    float k = 0;
	for (i = 0; i < 3; i++)
	{
        UpdateTree(x[i], root);

	}

	FILE *output = fopen ("Tree.dat", "wt");
	WalkTree (root, output);
	fclose (output);

   return 0;
}


/*
if(x < c->x_c)
{
    if (c->daughters[0] == NULL)
       {
        DivideCell(c);
        }
    UpdateTree(x,c->daughters[0]);
}
else
{
    if (c->daughters[1] == NULL)
    {
     DivideCell(c);

     }
    UpdateTree(x,c->daughters[1]);
}
*/
