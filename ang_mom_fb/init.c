/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sepy 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

double KentX (double);
float ran2 (long *);

/* ********************************************************************* */
void Init ( double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  v[RHO] = 0.0;
  v[PRS] = 0.0;
  EXPAND( v[VX1] = 0.0;  ,
          v[VX2] = 0.0;  ,
          v[VX3] = 0.0; )
}
/* ************************************************************************ */
void ICM_init(Data *d, Grid *grid)
/* Initializes the density, velocity and pressure by calculating the pressure
 * using Newton-Raphson. 
 * Written by Deovrat Prasad on 6th April, 2015.
 ************************************************************************** */
{
     int i, j, k, nv, ll, rank, cord, jj, kk;
     int mreq = 600;
     int n1p, n1m, ldim[DIMENSIONS], periods[DIMENSIONS], reord, nd;
     double *x1, *x2, *x3;
     double KENT[NX1_TOT], DM_POT[NX1_TOT];
     double const1, pg, fn, fnp, err, den, pout_r, pout_s;

//     print("%d %d %d %d\n", IBEG, IEND, NX1_TOT, prank);
     madd = 0.0; eadd =0.0; msub = 0.0;
#ifdef PARALLEL     
     MPI_Comm cart_comm;
     MPI_Request req[mreq];
     MPI_Status status;
     reord = 1;
     for(nd=0;nd<DIMENSIONS;nd++) ldim[nd]=grid[nd].nproc;
     D_EXPAND( periods[0]=0; ,
               periods[1]=0; ,
               periods[2]=1;  )
     MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, ldim, periods, reord, &cart_comm);
     MPI_Cart_shift( cart_comm , 0, 1, &n1m, &n1p);
#endif
     const1 = CONST_mu*CONST_mp*pow((CONST_mue*CONST_mp),(g_gamma-1.0))/(CONST_kB*1.16e7);
     x1 = grid[IDIR].x; x2 = grid[JDIR].x; x3 = grid[KDIR].x;
     ITOT_LOOP(i){
      DM_POT[i] = UNIT_VELOCITY*UNIT_VELOCITY*BodyForcePotential(x1[i], 0.0, 0.0); //convert to cgs
      KENT[i] = KentX( x1[i]*UNIT_LENGTH ); 
     }

     k = j = 0;
#ifdef PARALLEL
     if (grid[0].rank_coord == grid[0].nproc - 1){
#endif
     d->Vc[PRS][k][j][NX1_TOT-1] = KENT[NX1_TOT-1]*pow((g_inputParam[ne_out]*CONST_mue*CONST_mp),g_gamma)/const1;
#ifdef PARALLEL
     }

     for (ll = grid[0].nproc - 1; ll>=0; ll--){
       if (ll==grid[0].rank_coord){
          if (ll!= grid[0].nproc-1){
             MPI_Recv(&pout_r, 1, MPI_DOUBLE, n1p,12345, cart_comm, &status);
             d->Vc[PRS][k][j][IEND+3] = pout_r;
	     print("%20.10e %d %d \n", pout_r, grid[0].rank_coord, grid[1].rank_coord);
          }   
#endif
       for(i = NX1_TOT-2; i >=0; i--){
         pg = d->Vc[PRS][k][j][i+1];
         err = 1.e20;
         while(err >= 1.e-10){ //N-R to solve for pressure profile
         
         fn = d->Vc[PRS][k][j][i+1]-pg + 0.5*(DM_POT[i+1] - DM_POT[i])*pow(const1,1.0/g_gamma) \
         *(pow(d->Vc[PRS][k][j][i+1]/KENT[i+1],1.0/g_gamma)+pow(pg/KENT[i],1.0/g_gamma));

         fnp = -1.0 + 0.5*(DM_POT[i+1]-DM_POT[i])*(1.0/g_gamma) \
         *(pow(pg,(1.0/g_gamma-1.0))*(pow(const1/KENT[i],1.0/g_gamma)));

        pg = pg - fn/fnp;
        err = abs(fn/fnp)/pg;
       }
       d->Vc[PRS][k][j][i] = pg;
       }

       ITOT_LOOP(i) d->Vc[RHO][k][j][i] = pow(const1*d->Vc[PRS][k][j][i]/KENT[i],1.0/g_gamma);

#ifdef PARALLEL
    pout_s = d->Vc[PRS][k][j][IBEG+2];
    if (ll != 0) {
       MPI_Send(&pout_s, 1, MPI_DOUBLE, n1m, 12345, cart_comm);
    }
   }
  }  
#endif
   ITOT_LOOP(i) {
    d->Vc[RHO][0][0][i] /= UNIT_DENSITY; 
    d->Vc[PRS][0][0][i] /= (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
   }
   g_d_outer = d->Vc[RHO][0][0][IEND];
   g_p_outer = d->Vc[PRS][0][0][IEND];
   TOT_LOOP(k,j,i){
    d->Vc[RHO][k][j][i] = d->Vc[RHO][0][0][i];
    d->Vc[PRS][k][j][i] = d->Vc[PRS][0][0][i];
    //print("%20.10e %20.10e %d %d %d %d %d\n", d->Vc[RHO][k][j][i], d->Vc[PRS][k][j][i], k, j, i, grid[0].rank_coord, grid[1].rank_coord);
    //print("%20.10e %20.10e %d %d %d\n", d->Vc[RHO][k][j][i], d->Vc[PRS][k][j][i], k, j, i); 
    EXPAND( d->Vc[VX1][k][j][i]=0.0;  ,
               d->Vc[VX2][k][j][i]=0.0;  ,
               d->Vc[VX3][k][j][i]=0.0; )
   }

/* Introduction of Large Amplitude Random Pertubation in 
 * density at small scales.
 */
 int n,n1, l, l1, m, m1;
 long iseed = -2001;
 double kx, ky, kz;
 double phi, aklm, xs, ys, zs, pert;
 double delrho[NX3_TOT][NX2_TOT][NX1_TOT];

 memset(delrho, 0, sizeof delrho);
/* for(n = 4; n<=20; n++){
 for(n1 = -n; n1<=n; n1+=2*n){
    kx = 2.0*CONST_PI*n1/(2.*g_domEnd[IDIR]);
    for(l = 4; l<=20; l++){
    for(l1 = -l; l1<=l; l1+=2*l){
       ky = 2.0*CONST_PI*l1/(2.*g_domEnd[IDIR]);
       for(m = 4; m<=20; m++){
       for(m1 = -m; m1<=m; m1+=2*m){
          kz = 2.0*CONST_PI*m1/(2.*g_domEnd[IDIR]);
          pert=ran2(&iseed);
          phi = 2.0*CONST_PI*pert;
          pert=ran2(&iseed);          
          aklm = 0.15*g_inputParam[amp]*(0.5-pert)/sqrt(1.*pow(n,2)+pow(l,2)+pow(m,2));
          
          TOT_LOOP(k,j,i){
          xs = x1[i]*sin(x2[j])*cos(x3[k]);
          ys = x1[i]*sin(x2[j])*sin(x3[k]);
          zs = x1[i]*cos(x2[j]);
  
          delrho[k][j][i] += aklm*d->Vc[RHO][k][j][i]*cos(phi + kx*xs + ky*ys + kz*zs);
          }
       }} 
    }} 
 }}*/ 
 TOT_LOOP(k,j,i){
   d->Vc[RHO][k][j][i] += delrho[k][j][i];
   d->Vc[TRC][k][j][i] = 0.0;
   #if PHYSICS == MHD || PHYSICS == RMHD
   d->Vc[BX1][k][j][i] = 0.0;
   d->Vc[BX2][k][j][i] = 0.0;
   d->Vc[BX3][k][j][i] = 0.0;
   #endif
 } 
#ifdef PARALLEL
   MPI_Comm_free(&cart_comm);
#endif
}
/* ********************************************************************* */
double KentX (double r) //r is in cgs; returns X-ray entropy
{
  return g_inputParam[K0] + g_inputParam[K100]*pow(r/(CONST_pc*1.e5),g_inputParam[kalpha]);
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
  int k, j, i;
  double g_mass, g_TE, g_mdot,mdot_x, mdot_c, darea; 
  double g_KE1, g_KE2, g_KE3, g_mom1, g_mom2, g_mom3;
  double *dx1, *dx2, *dx3, *x1, *x2, *x3, d3;
  double sendarray[14], recvarray[14], dvol, TkeV;
  FILE *hist_file;

  x1 = grid[IDIR].x; x2 = grid[JDIR].x; x3 = grid[KDIR].x;
  dx1 = grid[IDIR].dx; dx2 = grid[JDIR].dx; dx3 = grid[KDIR].dx;
  #ifdef PARALLEL
  if (prank==0)
  #endif
  if (g_stepNumber==0) {
    hist_file = fopen ("pluto_hst.out", "w");
    fprintf(hist_file,"#time g_dt mass TE KE1 KE2 KE3 MOM1 MOM2 MOM3 mdot mdot_x mdot_c madd eadd mdot_jet msub\n ");
  }
  else hist_file = fopen ("pluto_hst.out", "a");

  g_mass=0.0; g_TE=0.0; g_KE1=0.0; g_KE2=0.0; g_KE3=0.0;
  g_mom1=0.0; g_mom2=0.0; g_mom3=0.0;
  g_mdot=0.0; mdot_x = 0.0; mdot_c = 0.0;

  DOM_LOOP(k,j,i){
        dvol = x1[i]*x1[i]*sin(x2[j])*dx1[i]*dx2[j]*dx3[k];
        g_mass += d->Vc[RHO][k][j][i]*dvol;
        g_TE += d->Vc[PRS][k][j][i]*dvol/(g_gamma-1.0);
        EXPAND( g_KE1 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i]*dvol;  ,
                g_KE2 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i]*dvol;  ,
                g_KE3 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*d->Vc[VX3][k][j][i]*dvol; )
        EXPAND( g_mom1 += d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*dvol;  ,
                g_mom2 += d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*dvol;  ,
                g_mom3 += d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*dvol; )
  }
  #ifdef PARALLEL
  if (grid[0].rank_coord == 0 ){
  #endif
  KDOM_LOOP(k){
  JDOM_LOOP(j){
  EXPAND( darea = 4.0*CONST_PI*x1[IBEG]*x1[IBEG];,
          darea = 2.0*CONST_PI*x1[IBEG]*x1[IBEG]*sin(x2[j])*dx2[j];,
          darea = x1[IBEG]*x1[IBEG]*sin(x2[j])*dx2[j]*dx3[k];)

  TkeV = (d->Vc[PRS][k][j][IBEG]*CONST_mu*CONST_mp)/(d->Vc[RHO][k][j][IBEG]*CONST_kB*1.16e9);

  g_mdot += d->Vc[RHO][k][j][IBEG]*d->Vc[VX1][k][j][IBEG]*darea;
  if(TkeV >= 0.5 && TkeV <=8.0){
    mdot_x += d->Vc[RHO][k][j][IBEG]*d->Vc[VX1][k][j][IBEG]*darea;
  }else if( TkeV < 0.5){
    mdot_c += d->Vc[RHO][k][j][IBEG]*d->Vc[VX1][k][j][IBEG]*darea;
  } 
  }}
  #ifdef PARALLEL
  } 
  #endif

  #ifdef PARALLEL
   sendarray[0]=g_mass; sendarray[1]=g_TE; sendarray[2]=g_KE1; sendarray[3]=g_KE2;
   sendarray[4]=g_KE3; sendarray[5]=g_mom1; sendarray[6]=g_mom2; sendarray[7]=g_mom3; 
   sendarray[8]=g_mdot; sendarray[9] = mdot_x; sendarray[10] = mdot_c;
   sendarray[11] = madd; sendarray[12] = eadd; sendarray[13] = mdot_jet; 
   sendarray[14] = msub;  
   MPI_Reduce (sendarray, recvarray, 15, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
   if (prank == 0){
     g_mass=recvarray[0]; g_TE=recvarray[1]; g_KE1=recvarray[2]; g_KE2=recvarray[3];
     g_KE3=recvarray[4]; g_mom1=recvarray[5]; g_mom2=recvarray[6]; g_mom3=recvarray[7]; 
     g_mdot=recvarray[8]; mdot_x = recvarray[9]; mdot_c = recvarray[10];
     madd = recvarray[11]; eadd = recvarray[12]; mdot_jet = recvarray[13];
     msub = recvarray[14];
  #endif
    fprintf(hist_file,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n ", g_time, g_dt, g_mass, g_TE, g_KE1, g_KE2, g_KE3, g_mom1, g_mom2, g_mom3, g_mdot, mdot_x, mdot_c, madd, eadd, mdot_jet, msub);

    fclose(hist_file);
  #ifdef PARALLEL
    }
  #endif
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3, *dx1, *dx2, *dx3;
  double DM_POT[NX1_TOT];
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;
  dx1 = grid[IDIR].dx;
  dx2 = grid[JDIR].dx;
  dx3 = grid[KDIR].dx;
 
   ITOT_LOOP(i){
   DM_POT[i] = BodyForcePotential(x1[i], 0.0, 0.0); 
   }
   if (side == 0) {    /* -- check solution inside domain -- */
    #ifdef PARALLEL
    if (grid[0].rank_coord == grid[0].nproc - 1){
    #endif
    KTOT_LOOP(k){
    JTOT_LOOP(j){
    d->Vc[RHO][k][j][IEND] = g_d_outer;
    d->Vc[PRS][k][j][IEND] = g_p_outer;
    }}
    #ifdef PARALLEL
    }
    #endif
    #ifdef PARALLEL
    if (grid[0].rank_coord == 0){
    #endif
    KTOT_LOOP(k){
    JTOT_LOOP(j){
    if( d->Vc[VX1][k][j][IBEG+1] >= 0.0){
            EXPAND( d->Vc[VX1][k][j][IBEG] = 0.0 ; ,
                    d->Vc[VX2][k][j][IBEG] = 0.0 ; ,
                    d->Vc[VX3][k][j][IBEG] = 0.0 ; )
        } 
    }}
    #ifdef PARALLEL
    }
    #endif
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
    #ifdef PARALLEL
    if (grid[0].rank_coord == 0){
    #endif
       BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i+1];
        if( d->Vc[VX1][k][j][IBEG] < 0.0){
           EXPAND( d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][i+1] ;  , 
                   d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][i+1] ;  ,
                   d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][i+1] ; )
        }else{
           EXPAND(d->Vc[VX1][k][j][i] = 0.0 ;,
                  d->Vc[VX2][k][j][i] = 0.0 ;,
                  d->Vc[VX3][k][j][i] = 0.0 ;)
        }
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i+1] + (d->Vc[RHO][k][j][i+1])*(DM_POT[i+1]-DM_POT[i]);
       }
    #ifdef PARALLEL
    }
    #endif
  }}
  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      #ifdef PARALLEL
      if (grid[0].rank_coord == grid[0].nproc - 1){
      #endif    
      BOX_LOOP(box,k,j,i){ 
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i-1];
        EXPAND( d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][2*IEND-i+1];,
               d->Vc[VX2][k][j][i] =  d->Vc[VX2][k][j][2*IEND-i+1];, 
               d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][2*IEND-i+1];)
        d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][j][i-1] - (d->Vc[RHO][k][j][i-1])*(DM_POT[i]-DM_POT[i-1]);
      }
     #ifdef PARALLEL
     }
     #endif
  }}

}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  double H0 = 2.19701879455606e-18, dcrit=3.*H0*H0/(8.*CONST_PI*CONST_G);             /* 390 kpc */
  double R200=pow(3.*g_inputParam[M200]/(4.*CONST_PI*200.*dcrit),(1./3.));
  double Phi, x;
  x = x1*UNIT_LENGTH/R200;
  Phi = -CONST_G*g_inputParam[M200]*log(1.+g_inputParam[c200]*x) \
         /((log(1.+g_inputParam[c200])-g_inputParam[c200]/(1.+g_inputParam[c200]))*R200*x);
  Phi += g_inputParam[Vc_sis]*g_inputParam[Vc_sis]*log(x1*UNIT_LENGTH/g_inputParam[r0_sis]);
  Phi -= CONST_G*g_inputParam[M_smbh]/(x1*UNIT_LENGTH);
  Phi /= UNIT_VELOCITY*UNIT_VELOCITY;
  return Phi;
}
#endif
/* ************************************************************************ */
float ran2(long *idum)

/*!
 *  Random number generator taken from Numerical Recipes in C.
 *  Returns a uniform deviate between 0.0 and 1.0.
 *  Call with idum a negative integer to initialize; 
 *  thereafter do not alter idum between successive deviates in a sequence.
 *  RNMX should approximate the largest floating value that is less than 1.
 ************************************************************************* */
{
 #define IM1 2147483563
 #define IM2 2147483399
 #define AM (1.0/IM1)
 #define IMM1 (IM1-1)
 #define IA1 40014
 #define IA2 40692
 #define IQ1 53668
 #define IQ2 52774
 #define IR1 12211
 #define IR2 3791
 #define NTAB 32
 #define NDIV (1+IMM1/NTAB)
 #define EPS 1.2e-7
 #define RNMX (1.0-EPS)

    int j ;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0){
       if (-(*idum) < 1) *idum =1;
       else *idum = -(*idum);
       idum2 = (*idum);
       for (j=NTAB + 7; j>=0; j--) {
           k=(*idum)/IQ1;
           *idum=IA1*(*idum-k*IQ1) - k*IR1;
           if (*idum < 0) *idum += IM1;
           if (j < NTAB) iv[j] = *idum;
       }
       iy = iv[0];
    }
    k=(*idum)/IQ1;
    *idum =IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 <0) idum2 += IM2;
    j = iy/NDIV;
    iy=iv[j]-idum2;
    iv[j]= *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}










