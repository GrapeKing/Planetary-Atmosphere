#include "pluto.h"

void Init (double *v, double x1, double x2, double x3)
{
  
  v[RHO] = exp(4.0/x1);
  v[VX1] = 0.0;
  g_isoSoundSpeed = 1.0;

}

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  int   i, j, k;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == X1_END){
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = exp(1.0/4.0);
        
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
  int i,j,k;
  char fname[512];
  FILE *fp;
  sprintf (fname, "%s/density.txt",RuntimeGet()->output_dir);
  fp = fopen(fname,"w");
  DOM_LOOP(k,j,i){
    
    fprintf(fp,"%.1e\n",k);
  }
  fclose(fp);
  printf("lol");
}