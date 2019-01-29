#include "pluto.h"
#include <stdlib.h>

double *mass_gas = NULL;
double *tot_Mg = NULL;
int num_div=0;

void Init (double *v, double x1, double x2, double x3)
{
  // Set initial conditions
  v[RHO] = g_inputParam[RHOINF]*exp(g_inputParam[BONDI]/x1);
  v[VX1] = 0.0;

  // Set global sound speed
  g_isoSoundSpeed = g_inputParam[CS];
}

void InitDomain (Data *d, Grid *grid) {
  // Allocate space for mass profile of atmosphere (only space for active domain)
  static int first_call = 1;
  if (first_call == 1) {
    first_call = 0;
    num_div = grid->np_int[IDIR];
    mass_gas = (double *) malloc(num_div*sizeof(double));
    tot_Mg = (double *) malloc(num_div*sizeof(double));
  }
}


void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {
  int   i, j, k;
  double  *x1 = grid->x[IDIR];

  // Impose no flow at outer boundary and fixed density at infinity
  if (side == X1_END){
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = g_inputParam[RHOINF]*exp(g_inputParam[BONDI]/x1[i]);
        d->Vc[VX1][k][j][i] = 0.0;
    }
  }

  // Boundary condition is reflective, all particles falling to the core will "bounce back"
  if (side == X1_BEG) {
      BOX_LOOP(box,k,j,i){
        //d->Vc[VX1][k][j][i] = 0.0;
    }
  }

  // Integrate the density profile over active cells
  if (side == 0) {
    DOM_LOOP(k,j,i) {
      // Get volume cell in spherical coordinates
      double dV = 4*M_PI*x1[i]*x1[i] * (grid->xr[IDIR][i] - grid->xl[IDIR][i]);

      mass_gas[i-2] = d->Vc[RHO][k][j][i]*dV;

      // Mass  profile goes from 0 to nDiv (only active domain, ignore boundary cells)
      if (i == 2) {
        tot_Mg[0] = mass_gas[0];
      } else {
        tot_Mg[i-2] = mass_gas[i-2] + tot_Mg[i-3];
      }
    }
  }
}

void BodyForceVector(double *v, double *g, double x1, double x2, double x3) {
  double Mg = 0;
  int i,j,k;
  static int pos=0;

  // Keep track of current cell position
  if (pos >= num_div) {
    pos = 0;
  }
  Mg = tot_Mg[pos];
  ++pos;

  // Get acceleration from gravity of the envelope
  g[IDIR] = -Mg*g_inputParam[G]/(x1*x1);
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}

double BodyForcePotential(double x1, double x2, double x3)
{
  // Get gravitational potential from core
  return -g_inputParam[BONDI]*g_isoSoundSpeed*g_isoSoundSpeed/x1;
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
