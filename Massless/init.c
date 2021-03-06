#include "pluto.h"

void Init (double *v, double x1, double x2, double x3)
{
  
  //v[RHO] = 0.5*exp(4.0/x1);
  v[RHO] = exp(1.0/4.0);
  v[VX1] = 0.0;
  g_isoSoundSpeed = 1.0;

  char fname[512];
  FILE *fp;

  sprintf (fname, "%s/density.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"w");
  fclose(fp);

  sprintf (fname, "%s/coordinate.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"w");
  fclose(fp);
}

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  int   i, j, k;
  double  *x1 = grid->x[IDIR];

  if (side == X1_END){
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = exp(4.0/x1[i]);
      }
  }
}

double BodyForcePotential(double x1, double x2, double x3)
{
  return -4.0/x1;
}

void InitDomain (Data *d, Grid *grid)
{
}

void Analysis (const Data *d, Grid *grid)
{
  // Will@Dao:
  int i,j,k;
  FILE *fp;
  static int first_call = 1; // static variable, keeps its value between calls of the function
  
  if (first_call == 1) {
    // Initializing the output files
    first_call = 0; // now zero for every subsequent call
    fp = fopen("density.txt","a");
    fprintf(fp,"%.4e ",g_domBeg[0]); // coordinate of the inner boundary
    DOM_LOOP(k,j,i){
      fprintf(fp,"%.4e ",grid->xr[IDIR][i]); // coordinate of every cell boundary
    }
    fprintf(fp,"\n"); // new line
    fclose(fp);
  }

  // exporting data (you can repeat for the velocity field)
  fp = fopen("density.txt","a");
  fprintf(fp,"%.4e ",g_time); // current time
  DOM_LOOP(k,j,i){
    fprintf(fp,"%.4e ",d->Vc[RHO][k][j][i]); // density row
  }
  fprintf(fp,"\n");
  fclose(fp);

  /*
  //I don't know how to execute it only once so I let it run anyway:P
  //Also I did not do this in init because it seems to loop over more cells than the DOM_LOOP command
  sprintf (fname, "%s/coordinate.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"w");
  DOM_LOOP(k,j,i){
    fprintf(fp,"%.4e ",x1[i]);
  }
  fclose(fp);
  */
  
}
