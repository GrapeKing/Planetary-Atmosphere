#include "pluto.h"
#include <stdlib.h>

double *mass_gas = NULL;
double *tot_Mg = NULL;
double *cor = NULL;

void Init (double *v, double x1, double x2, double x3) {
  
  v[RHO] = g_inputParam[RHOINF];
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
    TOT_LOOP(k,j,i) {
      cor[i] = grid->x[IDIR][i];
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
      d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][2*IEND+1-i];
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
      dV = 4.0/3.0*M_PI*(pow(grid->xr[IDIR][i],3) - pow(grid->xl[IDIR][i],3));
      mass_gas[i] = d->Vc[RHO][k][j][i]*dV; 
      if (i == IBEG) {
        tot_Mg[i] = mass_gas[i]; 
      } else {
        tot_Mg[i] = mass_gas[i] + tot_Mg[i-1];
      }
      //Will:
      //Outer buffer to help killing wave reflections
      if (g_time<4*g_domEnd[0]/g_isoSoundSpeed) {
	if (x1[i]>g_domEnd[0]-1.0) d->Vc[VX1][k][j][i] *= exp(-0.5*g_dt*g_isoSoundSpeed);
      }
    }
  }
}


// Tree-searching for an element in an ordered array
/* ************************************************* */
int binary_search(double* A, double x, int imin, int imax) {
  int imid;
  while (imin <= imax) {
    imid = (imin+imax)/2;
    if (x>=A[imid] && x<A[imid+1]) {
      return imid;
    } else if (x>A[imid]) {
      imin = imid + 1;
    } else {
      imax = imid - 1;
    }
  }
  return 0; //necessary when searching in an empty array
}



void BodyForceVector(double *v, double *g, double x1, double x2, double x3) {

  double Mg = 0;
  int i = binary_search(cor, x1, 0, NX1_TOT-1);
  Mg = tot_Mg[i]; //Will: removed the danger of tot_Mg[i-1] being out of bounds
  g[IDIR] = -Mg/(x1*x1);
}


// Gravity potential for the core alone
double BodyForcePotential(double x1, double x2, double x3) {
  double dbdt = 1.0; // d[Bondi]/dt
  double bondi = dbdt * g_time;
  if (bondi>g_inputParam[BONDI]) bondi = g_inputParam[BONDI];
  return -bondi/x1; 
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

