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

  int i, j, k;
  double *x1 = grid->x[IDIR];
  double dV;
  test=x1; // 

  // Impose no flow at outer boundary and fixed density at infinity

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

  // Integrate the density profile over active cells
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

  // Keep track of current cell position
  if (pos >= num_div) {
    pos = 0;
  }
  Mg = tot_Mg[pos];
  ++pos;

  // Will: the problem is that you don't know how many times BodyForceVector is actually called (if it includes the ghost cells,
  // if it is called on the N+1 boundaries of the N cells, if it is called 2 times per timestep if you use RK2 (definitions.h)...
  // The cleanest way is to identify the cell index 'i' given its location 'x1'; there are various ways to do this
  // (go through every cell (complexity N^2), use a binary tree (N.logN, optimal), look for neighbors of the previous cell (N but unsafe)

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
  // Will: you should open with "a" to add new rows at the end of the file,
  // so you can put lots of data in a single file, capturing the time evolution of the system.
  // you could export the radial grid first (boundaries best) and then the current time at the beginning of every row)
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
