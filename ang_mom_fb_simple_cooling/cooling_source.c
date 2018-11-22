/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Simplify original PLUTO cooling routine to be consistent with the 
  semi-implicit treatment of Sharma et al. 2010.

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
     - "THERMAL INSTABILITY WITH ANISOTROPIC THERMAL CONDUCTION AND ADIABATIC 
     COSMIC RAYS: IMPLICATIONS FOR COLD FILAMENTS IN GALAXY CLUSTERS" \n 
     Sharma, Parrish and Quataert, ApJ (2010) 720, 652
     
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
  double scrh, dtsub, q_m;
  double T0, n_e, n_i, unit_q, const_mui = 1./(1./CONST_mu - 1./CONST_mue);
  double LamCool(double);
  double xloc, xglob;
  
/* --------------------------------------------------------
    Set number and indices of the time-dependent variables
   -------------------------------------------------------- */

  unit_q = UNIT_DENSITY*pow(UNIT_VELOCITY,3.0); unit_q /= UNIT_LENGTH;

  DOM_LOOP(k,j,i){  /* -- span the computational domain -- */

    T0  = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*CONST_mu; //in cgs units
    scrh = LamCool(T0);
    n_e = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(CONST_mue*CONST_mp);
    n_i = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(const_mui*CONST_mp);//cgs number densities 
    scrh = n_e*n_i*scrh/unit_q; //code units
  /* ------------------------------------------
      Suggest next time step
     ------------------------------------------ */
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
    DOM_LOOP(k,j,i){  /* -- span the computational domain -- */ //first apply cooling
      #if INTERNAL_BOUNDARY == YES
       if (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY) continue;
      #endif
      if (d->flag[k][j][i] & FLAG_SPLIT_CELL) continue;

      T0 = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*CONST_mu; //in cgs units
      scrh = LamCool(T0);
      n_e = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(CONST_mue*CONST_amu);
      n_i = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(const_mui*CONST_amu);//cgs number densities
      q_m = n_e*n_i*scrh/unit_q; //code units
      d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i]/(1. + q_m*dtsub*(g_gamma-1.)/d->Vc[PRS][k][j][i]);
    }
  }
}
