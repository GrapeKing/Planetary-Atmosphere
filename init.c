#include "pluto.h"
#include <stdlib.h>

double *mass_gas = NULL;
double *tot_Mg = NULL;
double *cor = NULL;

void Init (double *v, double x1, double x2, double x3) {
  
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
  int i, j, k;
  if (first_call == 1) {
    first_call = 0;
    mass_gas = (double *) malloc(grid->np_tot[IDIR]*sizeof(double));
    tot_Mg = (double *) malloc(grid->np_tot[IDIR]*sizeof(double));
    cor = (double *) malloc(grid->np_tot[IDIR]*sizeof(double));
    TOT_LOOP(k,j,i) { //Will: replaced the DOM loop by a TOT loop (include ghost zones)
      cor[i] = grid->x[IDIR][i]; //Will: my advice: save the cell boundaries grid->xr[IDIR][i] instead of their center
    }
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

  int i, j, k;
  double Mg = 0;
  int temp0 = IBEG;
  int temp1 = IEND;
  int try = round((temp0+temp1)/2.0);
  double xmax = g_domEnd[IDIR];
  double xmin = g_domBeg[IDIR];

  //Will: it is a tree-based search ;-)
  // the 0.001 comparison is not elegant, but if it works...
  //(make sure it works with a basic python code if you can!)
  while (abs(x1 - cor[try])>0.001){
    if (x1 > cor[try]){
      temp0 = try;
      try = round((temp0+temp1)/2.0);
    } else{
      temp1 = try;
      try = round((temp0+temp1)/2.0);
    }
  }

  Mg = tot_Mg[try];

  g[IDIR] = -Mg/(x1*x1);
}

// Gravity potential for the core alone
double BodyForcePotential(double x1, double x2, double x3) {
  //Will: one of the reasons why you get oscillations is that the system is initialized a bit abruptly;
  // for example, if you start with a constant density and progressively increase the mass of the core
  // up to its nominal value, then the system is initialized smoothly and the oscillations will be much weaker;
  // another good trick would be to damp the velocity fluctuations manually near one boundary, inside a "buffer". 
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

