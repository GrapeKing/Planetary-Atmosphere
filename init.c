#include "pluto.h"
#include <stdlib.h>

double *mass_gas = NULL;
double *tot_Mg = NULL;
int num_div = 0;

void Init (double *v, double x1, double x2, double x3)
{
  // Set initial conditions
  v[RHO] = g_inputParam[RHOINF]; //Will: you may want to start with the analytic hydrostatic solution
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
    mass_gas = (double *) malloc(grid->np_tot[IDIR]*sizeof(double)); //Will: removed the 'num_div' thing
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
      //If the following is added I get oscillations:
      //Will: you still need to impose something; here you are extrapolating VX1 with a constant value,
      //but if the gas wants to fall, then there is a good chance that the velocity profile will not be flat,
      //for example with rho=constant, the velocity will tend to vary as 1/r^2 (div(rho.v)=0, no sound waves);
      //maybe try a linear extrapolation? Otherwise, are the oscillations sustained/amplified, or damped over time?
      //d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND];
    }
  }

  if (side == X1_BEG) {
      BOX_LOOP(box,k,j,i){
	d->Vc[RHO][k][j][i] = +d->Vc[RHO][k][j][2*IBEG-1-i]; //Will: symmetrizing the density, 
	d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][2*IBEG-1-i]; // anti-symmetrizing the velocity --> no mass flux
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
  double dx = 15.0/128.0; //Will: ouch, completely problem dependent! Use g_domBeg[IDIR], g_domEnd[IDIR], 
  double x0 = 1.0;

  temp = round((x1-x0)/dx);

  Mg = tot_Mg[IBEG+temp-1];
  //Here it should definitely be g_inputParam[BONDI]+Mg but the following produces the expected curve:P
  return -(g_inputParam[BONDI]-Mg)/x1; //Will: this is not the solution of Poisson's equation!
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

