#include "pluto.h"

void Init (double *v, double x1, double x2, double x3)
{
  // Will: this is where you can use g_inputParam[BONDI] in place of 4.0, so you don't have to re-compile the code
  // every time you want to sample a new value for the Bondi radius
  v[RHO] = exp(4.0/x1); 
  v[VX1] = 0.0;
  g_isoSoundSpeed = 1.0;

}

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  int   i, j, k;
  double  *x1 = grid->x[IDIR];
  
  // Will: so you are prescribing the density at the outer radial boundary
  // but not the velocity, so the gas may well be flowing through the domain,
  // especially if you don't enforce a solid boundary condition at X1_BEG!
  // 1. impose a boundary condition on VX1 as well,
  // 2. make sure the X1_BEG boundary prevents mass to pass through
  // (you can use the preset 'reflective', or (my advice) encode it yourself here
  if (side == X1_END){
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = exp(4.0/x1[i]);
    }
  }
  if (side == X1_BEG){
      BOX_LOOP(box,k,j,i){
        d->Vc[VX1][k][j][i] = 0.0;
    }
  }
}

double BodyForcePotential(double x1, double x2, double x3)
{
  return -4.0/x1; // Will: same advice about using user_def_parameters
}

void InitDomain (Data *d, Grid *grid)
{
}

void Analysis (const Data *d, Grid *grid)
{
  double  *x1 = grid->x[IDIR];
  int i,j,k;
  char fname[512];
  FILE *fp;
  
  sprintf (fname, "%s/density.txt",RuntimeGet()->output_dir);
  // Will: open with "a" to add new rows at the end of the file, 
  // if you store the data row-wise (one space instead of a new line between every element)
  // a bit annoying that gnuplot cannot read rows, but it will make it easier to have lot's of
  // data in a single file, capturing the time evolution of the system. 
  fp = fopen(fname,"w"); 
  DOM_LOOP(k,j,i){
    fprintf(fp,"%.4e\n",d->Vc[RHO][k][j][i]);
  }
  fclose(fp);
  
  sprintf (fname, "%s/error.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"w"); 
  DOM_LOOP(k,j,i){
    fprintf(fp,"%.4e\n",exp(4.0/x1[i])-d->Vc[RHO][k][j][i]);
  }
  fclose(fp);
  printf("lol");
}
