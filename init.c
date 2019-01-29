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
  // Will: try & understand this (the 'static' property and how it can be used):
  static int first_call = 1;
  if (first_call == 1) {
    first_call = 0;
    num_div = grid->np_int[IDIR];
    mass_gas = (double *) malloc(num_div*sizeof(double));
    tot_Mg = (double *) malloc(num_div*sizeof(double));
  }
}


void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {
  int i, j, k;
  double *x1 = grid->x[IDIR];
  double dV;
  test=x1; // 

  if (side == X1_END){
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = g_inputParam[RHOINF]*exp(g_inputParam[BONDI]/x1[i]);
        //d->Vc[VX1][k][j][i] = 0.0; // Lorenzo's
	// Will: VX1=0 at the outer boundary is one possibility, since indeed you want the system to reach a steady state
	// however if you set VX1=0 at both boundary (with the reflective condition at X1_BEG), then
	// you are allowing sound waves to bounce back and forth through the domain, so you might obtain a relatively slow convergence
	// + if the system wants to collapse, the X1_END boundary will somehow slow the infall.
	// Suggestion: extrapolate the velocity field (try to understand the following):
	d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND] + (x1[i] - x1[IEND]) * (d->Vc[VX1][k][j][IEND] - d->Vc[VX1][k][j][IEND-1])/(x1[IEND] - x1[IEND-1])
    }
  }
  // Will: you can keep this block even if you use 'reflective' in pluto.ini, 
  if (side == X1_BEG) {
      BOX_LOOP(box,k,j,i){
	// Will's suggestions: symmetrizing the density (even parity), anti-symmetrizing the velocity (odd parity)
	// (try to understand it! and of course, switch to 'userdef' in pluto.ini if you want to use it)
        d->Vc[RHO][k][j][i] = +d->Vc[RHO][k][j][2*IBEG-1-i];
	d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][2*IBEG-1-i];
    }
  }

  // Integrate the density profile
  if (side == 0) {
    // Will: you should still allocate mass_gas large enough to include the ghost cells (grid->np_tot[IDIR] instead of np_int)
    // since the DOM_LOOP will start at i=IBEG!=0, and finish at IEND=IBEG+np_int>np_int
    DOM_LOOP(k,j,i) {
      // volume element
      dV = 4*M_PI*x1[i]*x1[i] * (grid->xr[IDIR][i] - grid->xl[IDIR][i]);
      //printf("%f %f\n", dV, d->Vc[RHO][k][j][i]);
      mass_gas[i-2] = d->Vc[RHO][k][j][i]*dV; // Will: allocate a size np_tot*sizeof(double) and start at i, that will avoid confusions later

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
  // Will: so pos=0 initially, then you increment pos+1 but you're not sure how many steps will be incremented, 
  // and then you never reset it to zero (that's the 'static' property), so this will eventually crash;
  // the cleanest way is to identify the cell index 'i' given its location 'x1', there are various ways to do this.

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
