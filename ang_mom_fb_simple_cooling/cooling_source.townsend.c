/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Implement Townsend's "exact cooling" and heating 

  On output, a time-step estimate for the next time level is computed using
  the relative or absolute variation obtained during the integration of the ODE 
  (from t(n) --> t(n+1))  
  \f[
     \Delta t_c = \min_{ijk}\left[\frac{\Delta t^n M_R}{\epsilon}\right]
     \,\quad\rm{where}\quad
     \epsilon = \max\left(\left|\frac{p^{n+1}}{p^n} - 1\right|,\,
                |X^{n+1}-X^n|\right)
  \f]
  where \f$ M_R \f$ is the maximum cooling rate (defined by the global variable  
  ::g_maxCoolingRate) and X are the chemical species.
  
  \b References
     - "AN EXACT INTEGRATION SCHEME FOR RADIATIVE COOLING IN HYDRODYNAMICAL 
     SIMULATIONS" \n
     Townsend, ApJS (2009), 181, 391 

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Sharma 
  \date    November 12, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CoolingSource (const Data *d, double dt, Time_Step *Dts, Grid *GXYZ)
/*!
 * Integrate cooling and reaction source terms.
 *
 * \param [in,out]  d   pointer to Data structure
 * \param [in]     dt   the time step to be taken
 * \param [out]    Dts  pointer to the Time_Step structure
 * \param [in]    GXYZ  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int  k, j, i, ncool, cnt;
  double scrh, dtsub, unit_t, unit_q;
  double T0, n_e, n_i, n_tot;
  double Ytn, Yinv, const_mui = 1./(1./CONST_mu - 1./CONST_mue);
  double LamCool(double);
  double YCool(double);
  double YinvCool(double); 
  double xloc, xglob;

/*  ----------------------------------------------------------- 
                   first calculate the volume avg n_e n_i Lambda 
    -----------------------------------------------------------  */

  unit_t = UNIT_LENGTH/UNIT_VELOCITY;
  unit_q = UNIT_DENSITY*pow(UNIT_VELOCITY,3.0); unit_q /= UNIT_LENGTH;

  // first calculate the cooling time
  DOM_LOOP(k,j,i){
    T0  = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*CONST_mu; //in cgs units
    scrh = LamCool(T0);
    //print ("%12.6e %12.6e\n", T0, scrh);
    n_e = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(CONST_mue*CONST_amu);
    n_i = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(const_mui*CONST_amu);//cgs number densities 

    scrh = n_e*n_i*scrh/unit_q; //code units
    scrh = d->Vc[PRS][k][j][i]/(scrh*(g_gamma-1.0)); //code units
    Dts->dt_cool = MIN(Dts->dt_cool, scrh);
  } /* -- end loop on points -- */
  #ifdef PARALLEL
  xloc = Dts->dt_cool;
  MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  Dts->dt_cool = xglob;
  #endif
  ncool = ceil(dt/Dts->dt_cool); ncool = MAX(ncool, 50); // imposed maximum subcycling of 50
  dtsub = dt/ncool;

  for (cnt=0; cnt<ncool; ++cnt) {
    DOM_LOOP(k,j,i){  /* -- span the computational domain -- */
      #if INTERNAL_BOUNDARY == YES
       if (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY) continue;
      #endif
      if (d->flag[k][j][i] & FLAG_SPLIT_CELL) continue;
      T0  = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*CONST_mu; //in cgs units
      Ytn = YCool(T0);
      n_e = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(CONST_mue*CONST_amu);
      n_i = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(const_mui*CONST_amu);//cgs number densities      
      n_tot = n_e + n_i;
      scrh = Ytn - (g_gamma-1.)*(n_e*n_i/n_tot)*dtsub*unit_t/CONST_kB; //cgs unit
      scrh = YinvCool(scrh);
      if (scrh<=g_minCoolingTemp) scrh = g_minCoolingTemp;
      d->Vc[PRS][k][j][i] = scrh*d->Vc[RHO][k][j][i]/(KELVIN*CONST_mu); //in code units
    } /* -- end loop on points -- */
  } /* cnt */
}
