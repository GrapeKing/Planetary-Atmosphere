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
  double  *x1 = grid->x[IDIR];
  int i,j,k;
  char fname[512];
  FILE *fp;
  sprintf (fname, "%s/density.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"a");
  DOM_LOOP(k,j,i){
    fprintf(fp,"%.4e ",d->Vc[RHO][k][j][i]);
  }
  fprintf(fp,"\n");
  fclose(fp);
  
  sprintf (fname, "%s/coordinate.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"w");
  DOM_LOOP(k,j,i){
    fprintf(fp,"%.4e ",x1[i]);
  }
  fclose(fp);

  printf("lol");
}
