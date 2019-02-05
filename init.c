#include "pluto.h"
#include <stdlib.h>

double *mass_gas = NULL;
double *tot_Mg = NULL;
int num_div = 0;

void Init (double *v, double x1, double x2, double x3)
{
  // Set initial conditions
  v[RHO] = g_inputParam[RHOINF];
  v[VX1] = 0.0;

  // Set global sound speed
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
    num_div = grid->np_tot[IDIR];
    mass_gas = (double *) malloc(num_div*sizeof(double));
    tot_Mg = (double *) malloc(num_div*sizeof(double));
  }
}


void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {

  int i, j, k;
  double *x1 = grid->x[IDIR];
  double dV;

  if (side == X1_END){
    BOX_LOOP(box,k,j,i){
      d->Vc[RHO][k][j][i] = g_inputParam[RHOINF];
      //If the following is added I get oscillations:
      //d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND];
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
      // volume element
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
}

double BodyForcePotential(double x1, double x2, double x3)
{
  double Mg = 0;
  int temp = 0;
  double dx = 15.0/128.0;
  double x0 = 1.0;

  temp = round((x1-x0)/dx);

  Mg = tot_Mg[IBEG+temp-1];
//Here it should definitely be g_inputParam[BONDI]+Mg but the following produces the expected curve:P
  return -(g_inputParam[BONDI]-Mg)/x1;
}


void Analysis (const Data *d, Grid *grid)
{

  int i,j,k;
  FILE *fp;
  static int first_call = 1; 
  static int call = 0; 

  if (first_call == 1) {
    // Initializing the output files
    first_call = 0; // now zero for every subsequent call
    fp = fopen("coordinate.txt","a");

    //fprintf(fp,"%.4e ",g_domBeg[0]); // coordinate of the inner boundary
    DOM_LOOP(k,j,i){
      fprintf(fp,"%.4e ",grid->x[IDIR][i]); // coordinate of every cell boundary
    }
    fprintf(fp,"\n"); // new line
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

