/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the HD equations.
  
  The ConsToPrim() converts an array of conservative quantities to 
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total 
  energy unless zone has been tagged with FLAG_ENTROPY.
  In this case we recover pressure from conserved entropy:
  
      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)
  
  \author A. Mignone (mignone@ph.unito.it)
          B. Vaidya
  \date   April 14, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void PrimToCons (double **uprim, double **ucons, int ibeg, int iend)
/*!
 * Convert primitive variables to conservative variables. 
 *
 * \param [in]  uprim array of primitive variables
 * \param [out] ucons array of conservative variables
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final index of computation
 *
 *********************************************************************** */
{
  int  i, nv, status;
  double *v, *u;
  double rho, rhoe, T, gmm1;

  #if EOS == IDEAL
   gmm1 = g_gamma - 1.0;
  #endif
  for (i = ibeg; i <= iend; i++) {
  
    v = uprim[i];
    u = ucons[i];

    u[RHO] = rho = v[RHO];
    EXPAND (u[MX1] = rho*v[VX1];  ,
            u[MX2] = rho*v[VX2];  ,
            u[MX3] = rho*v[VX3];)

    #if EOS == IDEAL
     u[ENG] = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
     u[ENG] = 0.5*rho*u[ENG] + v[PRS]/gmm1;
    #elif EOS == PVTE_LAW
     status = GetPV_Temperature(v, &T);
     if (status != 0){
       T      = T_CUT_RHOE;
       v[PRS] = Pressure(v, T);
     }
     rhoe = InternalEnergy(v, T);

     u[ENG] = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
     u[ENG] = 0.5*rho*u[ENG] + rhoe;

     if (u[ENG] != u[ENG]){
       printf("PrimToCons\n");
       printf("KE:%12.6e uRHO : %12.6e, m2 : %12.6e \n",rhoe,v[RHO],u[ENG]);
       QUIT_PLUTO(1);
     }
    #endif

    for (nv = NFLX; nv < (NFLX + NSCL); nv++) u[nv] = rho*v[nv];
  }

}
/* ********************************************************************* */
int ConsToPrim (double **ucons, double **uprim, int ibeg, int iend, 
                unsigned char *flag, Grid *grid, int j, int k)
/*!
 * Convert from conservative to primitive variables.
 *
 * \param [in]  ucons  array of conservative variables
 * \param [out] uprim  array of primitive variables
 * \param [in]  beg    starting index of computation
 * \param [in]  end    final index of computation
 * \param [out] flag   array of flags tagging zones where conversion
 *                     went wrong.
 * 
 * \return Return (0) if conversion was succesful in every zone 
 *         [ibeg,iend]. 
 *         Otherwise, return a non-zero integer number giving the bit 
 *         flag(s) turned on during the conversion process.
 *         In this case, flag contains the failure codes of those
 *         zones where where conversion did not go through.
 *
 *********************************************************************** */
{
  int  i, nv, status=0, use_energy;
  double tau, rho, gmm1, rhoe, T;
  double kin, m2, rhog1, d3;
  double *u, *v;
  double *x1, *x2, *x3, *dx1, *dx2, *dx3;

  x1=grid[IDIR].x;
  x2=grid[JDIR].x;
  x3=grid[KDIR].x;
  dx1=grid[IDIR].dx;
  dx2=grid[JDIR].dx;
  dx3=grid[KDIR].dx;
  if(DIMENSIONS == 3){
    d3 = dx3[k];
  }else{
    d3= 2.0*CONST_PI;
  }
 
  #if EOS == IDEAL
   gmm1 = g_gamma - 1.0;
  #endif
 
  for (i = ibeg; i <= iend; i++) {

    flag[i] = 0;
    u = ucons[i];
    v = uprim[i];
   
    m2  = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);
    
  /* -------------------------------------------
           Check density positivity 
     ------------------------------------------- */
  
    if (u[RHO] < g_smallDensity) {
//      print("! ConsToPrim: rho < g_smallDensity (%8.2e), ", u[RHO]);
      Where (i, NULL);
      madd += (g_smallDensity - u[RHO])*x1[i]*x1[i]*dx1[i]*sin(x2[j])*dx2[j]*d3; 
      u[RHO]   = g_smallDensity;
      flag[i] |= RHO_FAIL;
    }

    v[RHO] = rho = u[RHO];
    tau = 1.0/u[RHO];
    EXPAND(v[VX1] = u[MX1]*tau;  ,
           v[VX2] = u[MX2]*tau;  ,
           v[VX3] = u[MX3]*tau;)
      
  /* --------------------------------------------
      Recover pressure from energy or entropy:
      1. IDEAL Equation of state
     -------------------------------------------- */

    #if EOS == IDEAL
   
     kin = 0.5*m2/u[RHO];
     use_energy = 1;     

     if (use_energy){
       if (u[ENG] < (g_smallPressure/gmm1 + kin)) {
         WARNING(
           print("! ConsToPrim: E < 0 (%8.2e), ", u[ENG]);
           Where (i, NULL);
         )
         u[ENG]   = g_smallPressure/gmm1 + kin;
         flag[i] |= ENG_FAIL;
       }else{
         v[PRS] = gmm1*(u[ENG] - kin);
         if (v[PRS] < g_smallPressure){
           WARNING(
             print("! ConsToPrim: p(E) < 0 (%8.2e), ", v[PRS]);
             Where (i, NULL);
           )
           v[PRS]   = g_smallPressure;
           flag[i] |= PRS_FAIL;
           u[ENG]   = v[PRS]/gmm1 + kin; /* -- redefine energy -- */
         }
       }
/*       if (u[ENG] > (g_largePressure/gmm1 + kin)) {
         WARNING(
           print("! ConsToPrim: E > 0 (%8.2e), ", u[ENG]);
           Where (i, NULL);
         )
         u[ENG]   = g_largePressure/gmm1 + kin;
         flag[i] |= ENG_FAIL;
       }else{
         v[PRS] = gmm1*(u[ENG] - kin);
         if (v[PRS] > g_largePressure){
           WARNING(
             print("! ConsToPrim: p(E) > 0 (%8.2e), ", v[PRS]);
             Where (i, NULL);
           )
           v[PRS]   = g_largePressure;
           flag[i] |= PRS_FAIL;
           u[ENG]   = v[PRS]/gmm1 + kin; 
         }
       } */ 

       //if (v[PRS]==0.0) print("here\n");
       //if (v[PRS]>g_largePressure) print("v[PRS], g_largePressure=%20.10e %20.10e\n", v[PRS], g_largePressure);

        v[PRS] = fmin(v[PRS],g_largePressure);

     }

    #endif  /* EOS == IDEAL */

  /* -- compute scalars  -- */

    for (nv = NFLX; nv < NVAR; nv++) v[nv] = u[nv]*tau;

    status |= flag[i];
  }
  return(status);
}
