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
  double RJET, SIGRJ, THJ, d3;

    x1=grid[IDIR].x; x2=grid[JDIR].x; x3=grid[KDIR].x;
    dx1=grid[IDIR].dx; dx2=grid[JDIR].dx; dx3=grid[KDIR].dx;
    
    RJET = g_inputParam[RJETKPC]*CONST_pc*1.e3/UNIT_LENGTH;  
    SIGRJ = g_inputParam[SIGR]*CONST_pc*1.e3/UNIT_LENGTH;
    THJ = g_inputParam[THJET]*CONST_PI/180.0;  
    norm_jet = 3.0/(2.0*CONST_PI*pow(RJET,3)*(1.0-cos(THJ)));
    
/* ***************Applying the mass source************** */    
    mdot_acc =0.0;
/* ************** Mass accretion based on PL extrapolation *************** */
/* get the density, temperature, grid data first */
    #ifdef PARALLEL
      MPI_Comm cart_comm, comm2d, comm1d;
      int nd, ldim[DIMENSIONS], periods[DIMENSIONS], reord, remain_dims[DIMENSIONS];
      reord = 1;
      for(nd=0;nd<DIMENSIONS;nd++) ldim[nd]=grid[nd].nproc;
      D_EXPAND( periods[0]=0; ,
                periods[1]=0; ,
                periods[2]=1;  )
      MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, ldim, periods, reord, &cart_comm);
      D_EXPAND( remain_dims[0]=0; ,
                remain_dims[1]=1; ,
                remain_dims[2]=1; )
      MPI_Cart_sub( cart_comm, remain_dims, &comm2d);
      D_EXPAND( remain_dims[0]=1; ,
                remain_dims[1]=0; ,
                remain_dims[2]=0; )
      MPI_Cart_sub( cart_comm, remain_dims, &comm1d);
    #endif
      double d_ew[NX1_TOT*grid[0].nproc], T_ew[NX1_TOT*grid[0].nproc], sh_vol[NX1_TOT*grid[0].nproc], dvol, k1[NVAR], T_keV, v0[NVAR];
      int nv;
 
      IDOM_LOOP(i){
	d_ew[i] = 0.0; T_ew[i] = 0.0; sh_vol[i] = 0.0;
        KDOM_LOOP(k){
        JDOM_LOOP(j){
           VAR_LOOP(nv) v0[nv] = d->Vc[nv][k][j][i];
           T_keV = (d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*KELVIN*CONST_mu/1.16e7; 
           if (T_keV > 0.5 && T_keV<8.0) {
             v0[RHOE] = d->Vc[PRS][k][j][i]/(g_gamma-1.);
             Radiat(v0, k1);
             dvol = -k1[RHOE]*D_SELECT (4.*CONST_PI*x1[i]*x1[i]*dx1[i]; ,
                            2.*CONST_PI*x1[i]*x1[i]*dx1[i]*sin(x2[j])*dx2[j]; ,
                            x1[i]*x1[i]*dx1[i]*sin(x2[j])*dx2[j]*dx3[k]; )
             sh_vol[i] += dvol;
             d_ew[i] += d->Vc[RHO][k][j][i]*dvol; //d in code units
             T_ew[i] += T_keV*dvol; //in keV
           }
        }
        }
      }
    double buf_grid[NX1_TOT*grid[0].nproc]; 
    #ifdef PARALLEL 
      double buf_in[3*NX1_TOT], buf_out[4*NX1_TOT], buf_out1[4*NX1_TOT*grid[0].nproc];
      for (i=0; i<(IEND-IBEG+1); i++){
        buf_in[i] = d_ew[i+IBEG]; buf_in[i+IEND-IBEG+1] = T_ew[i+IBEG]; 
        buf_in[i+2*(IEND-IBEG+1)] = sh_vol[i+IBEG];
      }
      MPI_Reduce( &buf_in[0], &buf_out[0], 3*(IEND-IBEG+1), MPI_DOUBLE, MPI_SUM, 0, comm2d );
      for (i=0; i<(IEND-IBEG+1); i++){
        buf_out[3*(IEND-IBEG+1)+i] = x1[i+IBEG];
      } 
      MPI_Gather( &buf_out[0], 4*(IEND-IBEG+1), MPI_DOUBLE, &buf_out1[0], 4*(IEND-IBEG+1), MPI_DOUBLE, 0, comm1d);
      int l;
      if (prank==0) {
        for (l=0; l<grid[0].nproc; l++) {
        for (i=0; i<(IEND-IBEG+1); i++){
          d_ew[i+IBEG+l*(IEND-IBEG+1)] = buf_out1[i+l*4*(IEND-IBEG+1)]/buf_out1[i+2*(IEND-IBEG+1)+l*4*(IEND-IBEG+1)];
          T_ew[i+IBEG+l*(IEND-IBEG+1)] = buf_out1[i+IEND-IBEG+1+l*4*(IEND-IBEG+1)]/buf_out1[i+2*(IEND-IBEG+1)+l*4*(IEND-IBEG+1)];
          buf_grid[i+IBEG+l*(IEND-IBEG+1)] = buf_out1[i+3*(IEND-IBEG+1)+l*4*(IEND-IBEG+1)];
        }}
      } 
    #else
      IDOM_LOOP(i){
        d_ew[i] /= sh_vol[i];
        T_ew[i] /= sh_vol[i];
        buf_grid[i] = x1[i];	
      } 
    #endif
/* now find the best-fit power law parameters*/

    if (prank==0) {
      double Delta, S, Sxx, Sx, Sxy, Sy, a_param, b_param;
      int  i_10kpc=IBEG;
      while ( buf_grid[i_10kpc]*UNIT_LENGTH < 1.e4*CONST_pc ) ++i_10kpc;
      S=Sx=Sy=Sxx=Sxy=0.0;
      for (i=IBEG; i<=i_10kpc; i++){
	S += 1.0; Sx += log(buf_grid[i]); Sy += log(d_ew[i]);
        Sxx += log(buf_grid[i])*log(buf_grid[i]); Sxy += log(buf_grid[i])*log(d_ew[i]); 
      } 
      Delta = S*Sxx - Sx*Sx; a_param = (Sxx*Sy-Sx*Sxy)/Delta; b_param = (S*Sxy-Sx*Sy)/Delta; //rho = exp(a)*x1**b 
      double r_Bondi, d_Bondi, T_Bondi=T_ew[IBEG]*1.16e7, c_s=sqrt(g_gamma*CONST_kB*T_Bondi/(CONST_mu*CONST_mp)); //isothermal
      double mdot_Bondi;
      r_Bondi = 2.*CONST_G*g_inputParam[M_smbh]/(c_s*c_s)/UNIT_LENGTH; //code units
      d_Bondi = exp(a_param)*pow(r_Bondi,b_param);
      mdot_Bondi = CONST_PI*c_s*d_Bondi*r_Bondi*r_Bondi/UNIT_VELOCITY; //code units
      mdot_jet = g_inputParam[EFF]*mdot_Bondi; // includes efficiency
    }
    #ifdef PARALLEL //broadcast mdot_jet
    MPI_Bcast (&mdot_jet, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Comm_free(&cart_comm);
    MPI_Comm_free(&comm2d);
    MPI_Comm_free(&comm1d);
    #endif

/*    KDOM_LOOP(k){
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
    
    if (mdot_acc < 0.0){
      mdot_jet = -g_inputParam[EFF]*mdot_acc; //this is actually epsilon*Mdot
    }else{
      mdot_jet = 0.0;
    } */
/* *********************************************************************** */
    DOM_LOOP(k,j,i){
       psi_jet = 0.25*(2.0 + tanh((THJ - x2[j])/g_inputParam[SIGTH]) \
       + tanh((THJ + x2[j]- CONST_PI)/g_inputParam[SIGTH])) \
       *(1.0 + tanh((RJET - x1[i])/SIGRJ));
       d->Vc[RHO][k][j][i] += g_dt*norm_jet*mdot_jet*psi_jet;
       d->Vc[VX1][k][j][i] += g_dt*norm_jet*mdot_jet*g_inputParam[VJET]*psi_jet\
       /(d->Vc[RHO][k][j][i]*UNIT_VELOCITY);
    }
}
