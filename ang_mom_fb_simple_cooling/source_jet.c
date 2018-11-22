/* ///////////////////////////////////////////////////////////////////// */
/*! 
 *   \file  
 *     \brief Contains basic function for AGN jet implementation.
 *
 *       The source_jet.c file contains source function for jet implementation  
 *       useful for problem configuration.
 *       It is automatically searched for by the makefile.
 *
 *       \author Deovrat Prasad (deovrat@physics.iisc.ernet.in)
 *               \date   May 12, 2015
 *               */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ********************************************************************* */
void Src_jet (const Data *d, Grid *grid)
{
  int i,j,k;
  double *x1, *x2, *x3, *dx1, *dx2, *dx3;
  double psi_jet, mdot_acc, darea, norm_jet;
  double sendarray[1], recvarray[1];
  double RJET, SIGRJ, THJ, d3, TinK, tau_dep=200.0*CONST_Myr;

    x1=grid[IDIR].x; x2=grid[JDIR].x; x3=grid[KDIR].x;
    dx1=grid[IDIR].dx; dx2=grid[JDIR].dx; dx3=grid[KDIR].dx;
    
    RJET = g_inputParam[RJETKPC]*CONST_pc*1.e3/UNIT_LENGTH;  
    SIGRJ = g_inputParam[SIGR]*CONST_pc*1.e3/UNIT_LENGTH;
    THJ = g_inputParam[THJET]*CONST_PI/180.0;  
    norm_jet = 3.0/(2.0*CONST_PI*pow(RJET,3)*(1.0-cos(THJ)));
    
/* ***************Applying the mass source************** */    
    mdot_acc =0.0;
/* ************** Mass accretion through the inner radius*************** */
    #ifdef PARALLEL
    if(grid[0].rank_coord == 0){
    #endif
    KDOM_LOOP(k){
    JDOM_LOOP(j){
    mdot_acc += x1[IBEG]*x1[IBEG]*sin(x2[j])*dx2[j]*dx3[k] \
              *d->Vc[RHO][k][j][IBEG]*d->Vc[VX1][k][j][IBEG];
    }}
    #ifdef PARALLEL
    }
    sendarray[0] = mdot_acc;
    MPI_Allreduce (sendarray, recvarray, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    mdot_acc = recvarray[0];
    #endif

    if (mdot_acc < 0.0){ //only inflowing mdot
      mdot_jet = -g_inputParam[EFF]*mdot_acc;
    }else{
      mdot_jet = 0.0;
    }
    #ifdef PARALLEL
    if (prank==0){
    #endif
      print ("%20.10e %20.10e\n", mdot_acc, mdot_jet);
    #ifdef PARALLEL
    }
    #endif 
/* *********************************************************************** */
    DOM_LOOP(k,j,i){
       psi_jet = 0.25*(2.0 + tanh((THJ - x2[j])/g_inputParam[SIGTH]) \
       + tanh((THJ + x2[j]- CONST_PI)/g_inputParam[SIGTH])) \
       *(1.0 + tanh((RJET - x1[i])/SIGRJ));
       d->Vc[RHO][k][j][i] += g_dt*norm_jet*mdot_jet*psi_jet;
       d->Vc[VX1][k][j][i] += g_dt*norm_jet*mdot_jet*g_inputParam[VJET]*psi_jet\
       /(d->Vc[RHO][k][j][i]*UNIT_VELOCITY);

       d->Vc[TRC][k][j][i] += g_dt*norm_jet*mdot_jet*psi_jet/d->Vc[RHO][k][j][i];
       TinK   = d->Vc[PRS][k][j][i]*CONST_mu*CONST_mp/(d->Vc[RHO][k][j][i]*CONST_kB);
       if (TinK < 5.e4 && d->Vc[RHO][k][j][i] > 1.e-24){
         msub += x1[i]*x1[i]*sin(x2[j])*dx2[j]*dx3[k]*dx1[i]*d->Vc[RHO][k][j][i]*( 1.0- exp(-g_dt/tau_dep) );
         d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i]*exp(-g_dt/tau_dep);
       }
    }
}
