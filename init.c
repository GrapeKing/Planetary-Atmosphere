#include "pluto.h"
#include <stdlib.h>

double *mass_gas = NULL;
double *tot_Mg = NULL;
double *pos_grid = NULL;
double _x1;
int cnt=2;
int num_div=0;

void Init (double *v, double x1, double x2, double x3)
{
  // Will: this is where you can use g_inputParam[BONDI] in place of 4.0, so you don't have to re-compile the code
  // every time you want to sample a new value for the Bondi radius
  v[RHO] = g_inputParam[RHOINF]*exp(g_inputParam[BONDI]/x1);
  v[VX1] = 0.0;
  g_isoSoundSpeed = g_inputParam[CS];

  // Count number of points in domain (uses fact this is 1D)
  num_div += 1;
}

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
{
  int   i, j, k;
  double  *x1 = grid->x[IDIR];
  double  *dx = grid->dx[IDIR];
  double  *dy = grid->dx[JDIR];
  double  *dz = grid->dx[KDIR];

  // Will: so you are prescribing the density at the outer radial boundary
  // but not the velocity, so the gas may well be flowing through the domain,
  // especially if you don't enforce a solid boundary condition at X1_BEG!
  // 1. impose a boundary condition on VX1 as well in the X1_END loop
  // 2. make sure the X1_BEG boundary prevents mass to pass through
  // (you can use the preset 'reflective', or (my advice) encode it yourself here
  if (side == X1_END){
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = g_inputParam[RHOINF]*exp(g_inputParam[BONDI]/x1[i]);
    }
  }
  // Will: same here, you should also impose a condition on RHO in the X1_BEG loop
  // to make sure the boundary conditions are under control
  if (side == X1_BEG){
      BOX_LOOP(box,k,j,i){
        d->Vc[VX1][k][j][i] = 0.0;
    }
  }

  pos_grid = x1;
  if (side == 0) {
    TOT_LOOP(k,j,i) {
      mass_gas[i] = d->Vc[RHO][k][j][i]*dx[i]*dy[j]*dz[k];
      //printf("%f\n", mass_gas[i]);
    }
    /*tot_Mg[0] = mass_gas[0];
    TOT_LOOP(k,j,i) {
      if (i==0) continue;
      tot_Mg[i] = mass_gas[i] + tot_Mg[i-1];
    }*/
  }
}

/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3) {
  double Mc = g_inputParam[BONDI]*g_isoSoundSpeed*g_isoSoundSpeed/g_inputParam[G];
  double Mg = 0;
  int i,j,k,pos;

  TOT_LOOP(k,j,i) {
    if (pos_grid[i] == x1) {
      pos = i;
    }
  }
  for(i=0; i<=pos; ++i) {
    Mg += mass_gas[i];
  }
  //printf("%d, %f\n", pos, Mg);

  g[IDIR] = -(Mc+Mg)*g_inputParam[G]/(x1*x1);
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}

double BodyForcePotential(double x1, double x2, double x3)
{
  double Mc = g_inputParam[BONDI]*g_isoSoundSpeed*g_isoSoundSpeed/g_inputParam[G];
  double Mg = 0;
  int i,j,k;

  TOT_LOOP(k,j,i) {
    if (pos_grid[i] == x1) {
      Mg = tot_Mg[i];
    }
  }

  return 0; // Will: same advice about using user_def_parameters
}

void InitDomain (Data *d, Grid *grid)
{
  pos_grid = grid->x[IDIR];

  // Create a global array to store the mass of the gas cloud
  mass_gas = (double *) malloc(num_div*sizeof(double));
  tot_Mg = (double *) malloc(num_div*sizeof(double));
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
    fprintf(fp,"%.4e %.4e %.4e\n", x1[i], d->Vc[RHO][k][j][i], g_inputParam[RHOINF]*exp(g_inputParam[BONDI]/x1[i]));
  }
  fclose(fp);

  sprintf (fname, "%s/massProf.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"w");
  DOM_LOOP(k,j,i){
    fprintf(fp,"%.4e %.4e\n", x1[i], mass_gas[i]);
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
