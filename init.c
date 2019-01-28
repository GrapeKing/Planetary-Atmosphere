#include "pluto.h"
#include <stdlib.h>

double *mass_gas = NULL;
double *tot_Mg = NULL;
double *test;
int num_div=0; // Will: see comment in Init()

void Init (double *v, double x1, double x2, double x3)
{
  // Count number of points in domain (uses fact this is 1D)
  //num_div += 1;
  // Will: this is not going to work, the Init function might be called several times per cell.
  // The number of points in the domain is contained in grid->np_int[IDIR]
  // you can find the attributes of the grid structure at
  // http://plutocode.ph.unito.it/Doxygen/API-Reference_Guide/struct_grid.html

  v[RHO] = g_inputParam[RHOINF]*exp(g_inputParam[BONDI]/x1);
  v[VX1] = 0.0;
  g_isoSoundSpeed = g_inputParam[CS];
}

void InitDomain (Data *d, Grid *grid) {
  // Will: try & understand this:
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
  test=x1;

  // Will: you are prescribing the density at the outer radial boundary
  // but not the velocity, so the gas may well be flowing through the domain,
  // especially if you don't enforce a solid boundary condition at X1_BEG.
  // 1. impose a boundary condition on RHO and VX1 at X1_BEG and X1_END
  // 2. make sure the X1_BEG boundary prevents mass to pass through (not 'outflow' !)
  if (side == X1_END){
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = g_inputParam[RHOINF]*exp(g_inputParam[BONDI]/x1[i]);
        d->Vc[VX1][k][j][i] = 0.0;
    }
  }
  // Will: same here, you should also impose a condition on RHO in the X1_BEG loop
  // to make sure the boundary conditions are under control
  if (side == X1_BEG) {
      BOX_LOOP(box,k,j,i){
        //d->Vc[VX1][k][j][i] = 0.0;
    }
  }

  // Integrate the density profile
  if (side == 0) {
    // /!\ TOT_LOOP will loop over the ghost+active cells
    // --> you should allocate mass_gas large enough to include the ghost cells (grid->np_tot[IDIR])
    DOM_LOOP(k,j,i) {
      // Will: two things:
      // 1. since you're doing 1D, the Y and Z coordinates are not really used,
      // so you can/should remove the dy[j] and dz[k] parts to avoid generating
      // weird results if you ever change the coordinate boundaries in pluto.ini
      // 2. in spherical coordinates, the volume element is not dx*dy*dz;
      // since you are doing 1D, you have already integrated over the angles and all you need is
      // dV = 4*M_PI*x[i]*x[i] * (grid->xr[i] - grid->xl[i])

      double dV = 4*M_PI*x1[i]*x1[i] * (grid->xr[IDIR][i] - grid->xl[IDIR][i]);
      //printf("%f %f\n", dV, d->Vc[RHO][k][j][i]);
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

  if (pos >= num_div) {
    pos = 0;
  }
  Mg = tot_Mg[pos];
  //printf("%d, %f %f\n", pos, Mg, test[pos+2]);
  ++pos;


  // Find the mass of the envelope
  // Will: so for every coordinate x1, you are loopîng over every cell index i
  // i.e. a full N^2 complexity... you can optimize this in several ways
  // Will: same thing, for every coordinate you are looping over every cell index,
  // giving you N(N-1)/2 number of operations for N cells;
  // why not compute the integral just once in UserDefBoundary:side==0 ?

  g[IDIR] = -Mg*g_inputParam[G]/(x1*x1);
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}

double BodyForcePotential(double x1, double x2, double x3)
{
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
