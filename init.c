#include "pluto.h"
#include <stdlib.h>

double *mass_gas = NULL;
double *tot_Mg = NULL;

void Init (double *v, double x1, double x2, double x3)
{
  v[RHO] = g_inputParam[RHOINF]*exp(-g_inputParam[BONDI]/g_domEnd[IDIR])*exp(g_inputParam[BONDI]/x1);
  v[VX1] = 0.0;

  g_isoSoundSpeed = g_inputParam[CS];

  char fname[512];
  FILE *fp;

  sprintf (fname, "%s/density.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"w");
  fclose(fp);

  sprintf (fname, "%s/coordinate.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"w");
  fclose(fp);
}

void InitDomain (Data *d, Grid *grid) {
  static int first_call = 1;
  if (first_call == 1) {
    first_call = 0;
    mass_gas = (double *) malloc(grid->np_tot[IDIR]*sizeof(double));
    tot_Mg = (double *) malloc(grid->np_tot[IDIR]*sizeof(double));
  }
}


void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {

  int i, j, k;
  double *x1 = grid->x[IDIR];
  double dV;

  if (side == X1_END){
    BOX_LOOP(box,k,j,i){
      d->Vc[RHO][k][j][i] = g_inputParam[RHOINF];
      d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND] + (x1[i] - x1[IEND]) * (d->Vc[VX1][k][j][IEND] - d->Vc[VX1][k][j][IEND-1]) / (x1[IEND] - x1[IEND-1]);
    }
  }

  if (side == X1_BEG) {
      BOX_LOOP(box,k,j,i){
	d->Vc[RHO][k][j][i] = +d->Vc[RHO][k][j][2*IBEG-1-i];
	d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][2*IBEG-1-i];
    }
  }

  if (side == 0) {
    DOM_LOOP(k,j,i) {

      dV = 4*M_PI*x1[i]*x1[i] * (grid->xr[IDIR][i] - grid->xl[IDIR][i]);
      mass_gas[i] = d->Vc[RHO][k][j][i]*dV; 

      if (i == IBEG) {
        tot_Mg[i] = mass_gas[i]; 
      } else {
        tot_Mg[i] = mass_gas[i] + tot_Mg[i-1];
      }
    }
  }
}

void BodyForceVector(double *v, double *g, double x1, double x2, double x3) {

  /* Will: problems with this implementation:
  / 1. 15.0 and 128.0 are problem dependent = danger!
  / --> Use g_domEnd[IDIR]-g_domBeg[IDIR], and extract the 128 from another function (global variable?)
  / 2. You are using a logarithmic grid ('l+' in pluto.ini), so (x1-x0)/dx is not at the appropriate cell!
  / 3. optional: this version won't be parallelized well...
  /
  / My advice: design a function with prototype that takes x1 as a parameter and returns the integer index of the corresponding cell; 
  / you will need the location of the cells on a global context, an then either: 
  / loop over all indices (complexity N), use a binary tree (logN, optimal), or check the neighbors of the previously used cell (complexity 1, but sometimes unsafe)
  */
  int i, j, k;
  double Mg = 0;
  int temp = 0;
  double dx = (g_domEnd[IDIR]-g_domBeg[IDIR])/256.0;
  double x0 = 1.0;

  temp = round((x1-x0)/dx);

  Mg = tot_Mg[IBEG+temp-1];

  g[IDIR] = -Mg/(x1*x1);
}

// Gravity potential for the core alone
double BodyForcePotential(double x1, double x2, double x3) {
  return -(g_inputParam[BONDI])/x1; 
}


void Analysis (const Data *d, Grid *grid) {

  int i,j,k;
  FILE *fp;
  static int first_call = 1; 
  static int call = 0; 

  if (first_call == 1) {
    first_call = 0;
    fp = fopen("coordinate.txt","a");

    //fprintf(fp,"%.4e ",g_domBeg[0]); // coordinate of the inner boundary
    DOM_LOOP(k,j,i){
      fprintf(fp,"%.4e ",grid->x[IDIR][i]); // coordinate of every cell boundary
    }
    fprintf(fp,"\n");
    fclose(fp);
  }


  fp = fopen("density.txt","a");
  // current time
  //fprintf(fp,"%.4e ",g_time); 
  DOM_LOOP(k,j,i){
    fprintf(fp,"%.4e ",d->Vc[RHO][k][j][i]); // density row
  }
  fprintf(fp,"\n");
  fclose(fp);
}

