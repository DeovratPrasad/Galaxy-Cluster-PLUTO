#include "pluto.h"
/* ***************************************************************** */
double YCool (double T)
/*!
 *   returns Y[T] in cgs units.
 * 
 ******************************************************************* */
{
  int    klo, khi, kmid;
  static int ntab;
  double  Tmid, scrh, dT;
  static double *L_tab, *T_tab, *Y_tab;

  FILE *fcool;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

  if (T_tab == NULL){
    print (" > Reading table from disk...\n");
    fcool = fopen("ct_townsend.dat","r");
    if (fcool == NULL){
      print ("! YCool: ct_townsend.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);
    Y_tab = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf %lf\n", T_tab + ntab,
                                       L_tab + ntab, Y_tab + ntab)!=EOF) {
      ntab++;
    }
  }

/* ---------------------------------------------
            Make sure that T is well-defined 
   --------------------------------------------- */

  if (T != T){
    printf (" ! Nan found in YCool \n");
    printf (" ! T = %12.6e\n", T);
    QUIT_PLUTO(1);
  }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  klo = 0;
  khi = ntab - 1;

  if (T < T_tab[klo]){
    //print (" T<Tlo    %12.6e\n",T);
    //QUIT_PLUTO(1);
    return Y_tab[klo];
  }

  if (T > T_tab[khi]){
    //print (" T>Thi   %12.6e\n",T);
    //QUIT_PLUTO(1);
    return Y_tab[khi];
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Tmid = T_tab[kmid];
    if (T <= Tmid){
      khi = kmid;
    }else if (T > Tmid){
      klo = kmid;
    }
  }

/* -----------------------------------------------
    Compute and return Y 
   ----------------------------------------------- */

  dT       = T_tab[khi] - T_tab[klo];
  scrh     = Y_tab[klo]*(T_tab[khi] - T)/dT + Y_tab[khi]*(T - T_tab[klo])/dT;
  return scrh;
}


/* ***************************************************************** */
double YinvCool (double Y)
/*!
 *   returns T in cgs units.
 *
 ******************************************************************* */
{
  int    klo, khi, kmid;
  static int ntab;
  double  Ymid, scrh, dY;
  static double *L_tab, *T_tab, *Y_tab;

  FILE *fcool;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

  if (T_tab == NULL){
    print (" > Reading table from disk...\n");
    fcool = fopen("ct_townsend.dat","r");
    if (fcool == NULL){
      print ("! YCool: ct_townsend.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);
    Y_tab = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf %lf\n", T_tab + ntab,
                                       L_tab + ntab, Y_tab + ntab)!=EOF) {
      ntab++;
    }
  }
  /* ---------------------------------------------
            Make sure that T is well-defined
   --------------------------------------------- */

  if (Y != Y){
    printf (" ! Nan found in YinvCool \n");
    printf (" ! Y = %12.6e\n", Y);
    QUIT_PLUTO(1);
  }

  /* ----------------------------------------------
        Table lookup by binary search
   ---------------------------------------------- */

  klo = 0;
  khi = ntab - 1;

  if (Y < Y_tab[klo]){
    //print (" Y<Ylo    %12.6e\n",Y);
    //QUIT_PLUTO(1);
    return T_tab[klo];
  }

  if (Y > Y_tab[khi]){
    //print (" Y>Yhi   %12.6e\n",Y);
    //QUIT_PLUTO(1);
    return T_tab[khi];
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Ymid = Y_tab[kmid];
    if (Y <= Ymid){
      khi = kmid;
    }else if (Y > Ymid){
      klo = kmid;
    }
  }

/* -----------------------------------------------
    Compute and return Lambda
   ----------------------------------------------- */

  dY       = Y_tab[khi] - Y_tab[klo];
  scrh     = T_tab[klo]*(Y_tab[khi] - Y)/dY + T_tab[khi]*(Y - Y_tab[klo])/dY;
  return scrh;
}
