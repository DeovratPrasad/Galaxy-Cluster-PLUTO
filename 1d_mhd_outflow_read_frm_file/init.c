/* ///////////////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
double ran2(long *);
void Get_ICMDP(double radius, double *, double *);
void Get_Del(double radius, double theta, double phi, double *);

void Get_Del( double r, double th, double ph, double *del){ //returns del at a given (r,t,p)
  static long iseed=-2001;
  double pert, kx, ky, kz, xs, ys, zs, den=0.0;
  int nmin=4, nmax=20;
  int l, m, n, l1, m1, n1;
  static double ***A_lmn, ***phi_lmn; //3-D array of amplitudes, phases

  if (A_lmn==NULL && phi_lmn==NULL){ //calculate the coeffts. in the first step 
    A_lmn = ARRAY_3D(2*(nmax-nmin+1), 2*(nmax-nmin+1), 2*(nmax-nmin+1),  double);
    phi_lmn = ARRAY_3D(2*(nmax-nmin+1), 2*(nmax-nmin+1), 2*(nmax-nmin+1),  double);

    for(l=0; l<2*(nmax-nmin)+2; l++){
      if (l<=nmax-nmin){ //l1 gives correct wavenumber
        l1 = l + nmin;
      } else {
        l1 = l-2*nmax+nmin-1;
      }
    for(m=0; m<2*(nmax-nmin)+2; m++){
      if (m<=nmax-nmin){ //m1 gives correct wavenumber
        m1 = m + nmin;
      } else {
        m1 = m-2*nmax+nmin-1;
      }
    for(n=0; n<2*(nmax-nmin)+2; n++){
      if (n<=nmax-nmin){ //n1 gives correct wavenumber
        n1 = n + nmin;
      } else {
        n1 = n-2*nmax+nmin-1;
      }

      pert = ran2(&iseed);
      phi_lmn[n][m][l] = 2.0*CONST_PI*pert;
      pert = ran2(&iseed);
      A_lmn[n][m][l] = 0.15*g_inputParam[amp]*(0.5-pert)/sqrt(1.*(n1*n1+l1*l1+m1*m1));

    }}}
  }

    for(l=0; l<2*(nmax-nmin)+2; l++){
      if (l<=nmax-nmin){ //l1 gives correct wavenumber
        l1 = l + nmin;
      } else {
        l1 = l-2*nmax+nmin-1;
      }
      kx = 2.0*CONST_PI*l1/(2.*g_inputParam[RMAX]);
    for(m=0; m<2*(nmax-nmin)+2; m++){
      if (m<=nmax-nmin){ //m1 gives correct wavenumber
        m1 = m + nmin;
      } else {
        m1 = m-2*nmax+nmin-1;
      }
      ky = 2.0*CONST_PI*m1/(2.*g_inputParam[RMAX]);
    for(n=0; n<2*(nmax-nmin)+2; n++){
      if (n<=nmax-nmin){ //n1 gives correct wavenumber
        n1 = n + nmin;
      } else {
        n1 = n-2*nmax+nmin-1;
      }
      kz = 2.*CONST_PI*n1/(2.*g_inputParam[RMAX]);


      xs = r*sin(th)*cos(ph);
      ys = r*sin(th)*sin(ph);
      zs = r*cos(th);

      den += A_lmn[n][m][l]*cos(phi_lmn[n][m][l] + kx*xs + ky*ys + kz*zs);
    }}}

    *del = den;
}

void Get_ICMDP( double r, double *rho, double *prs){ //returns the smooth profiles; pert added separately
  char fname[256] = "profiles.dat";
  int i, j, c, cnt, r_indx, NVals = 10000;
  double wt_r, lnx, lnR_beg, lnR_end, dln_R;
  static double *Rin, *Pin, *Rhoin;
  FILE *fp;

  /* Read the data from the profiles.dat file generated from Python */
  if (Rin == NULL){
    fp = fopen(fname, "r");
    Rin = ARRAY_1D(NVals, double);
    Pin = ARRAY_1D(NVals, double);
    Rhoin = ARRAY_1D(NVals, double);
    cnt = 0;
    while(cnt < NVals && fscanf(fp,"%le %le %le",&Rin[cnt], &Rhoin[cnt], &Pin[cnt]) == 3){
      cnt += 1;
    }
    fclose(fp);
  }

  /* Do Linear Interpolation */
  lnR_beg = log(Rin[0]);
  lnR_end = log(Rin[NVals-1]);
  dln_R = (lnR_end - lnR_beg)/(NVals-1); // Uniform in Logspace.
  lnx = log(r); //table is in code units  
  if (lnx < lnR_beg || lnx>lnR_end) {
    print("initial radius out of bounds, lnx, lnR_beg, lnR_end = %le %le %le\n", lnx, lnR_beg, lnR_end);
    QUIT_PLUTO(1);
  }
  else {
    r_indx = (int) ceil((lnx - lnR_beg)/dln_R);
    wt_r = (log(Rin[r_indx]) - lnx)/dln_R;
  }
  *rho = Rhoin[r_indx]*(1.0 - wt_r) + Rhoin[r_indx-1]*wt_r;
  *prs = Pin[r_indx]*(1.0 - wt_r) + Pin[r_indx-1]*wt_r;
}
/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 *********************************************************************** */
{
  double dicm, picm, delta;

  Get_Del(x1, x2, x3, &delta); //returns perturbations at the present grid center 

  Get_ICMDP(x1, &dicm, &picm); //returns density, p at a given location in code units
  v[RHO] = dicm*(1.+delta);
  #if HAVE_ENERGY
   v[PRS] = picm;
  #endif

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[TRC] = 0.0;

  #if PHYSICS == MHD
   v[BX1] = 0.0; v[BX2] = 0.0; v[BX3] = 0.0;
   v[AX1] = 0.0; v[AX2] = 0.0; v[AX3] = 0.0;
  #endif
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
  double g_mcold, g_mhot;
  double *dx1, *dx2, *dx3, *x1, *x2, *x3, d3;
  double sendarray[17], recvarray[17], dvol, TkeV;
  FILE *hist_file;

  x1 = grid[IDIR].x; x2 = grid[JDIR].x; x3 = grid[KDIR].x;
  dx1 = grid[IDIR].dx; dx2 = grid[JDIR].dx; dx3 = grid[KDIR].dx;
  #ifdef PARALLEL
  if (prank==0)
  #endif
  if (g_stepNumber==0) {
    hist_file = fopen ("pluto_hst.out", "w");
    fprintf(hist_file,"#[1]time [2]g_dt [3]mass [4]TE [5]KE1 [6]KE2 [7]KE3 [8]MOM1 [9]MOM2 [10]MOM3 [11]mdot [12]mdot_x [13]mdot_c [14]madd [15]eadd [16]mdot_jet [17]msub [18]mcold [19]mhot\n ");
  }
  else hist_file = fopen ("pluto_hst.out", "a");

  g_mass=0.0; g_TE=0.0; g_KE1=0.0; g_KE2=0.0; g_KE3=0.0;
  g_mom1=0.0; g_mom2=0.0; g_mom3=0.0;
  g_mdot=0.0; mdot_x = 0.0; mdot_c = 0.0; g_mcold = 0.0; g_mhot=0.0;

  DOM_LOOP(k,j,i){
        EXPAND ( dvol = 4.*CONST_PI*x1[i]*x1[i]*dx1[i];,
                 dvol = 2.*CONST_PI*x1[i]*x1[i]*sin(x2[j])*dx1[i]*dx2[j];,
                 dvol = x1[i]*x1[i]*sin(x2[j])*dx1[i]*dx2[j]*dx3[k];)
        g_mass += d->Vc[RHO][k][j][i]*dvol;
        g_TE += d->Vc[PRS][k][j][i]*dvol/(g_gamma-1.0);
        EXPAND( g_KE1 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i]*dvol;  ,
                g_KE2 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i]*dvol;  ,
                g_KE3 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*d->Vc[VX3][k][j][i]*dvol; )
        EXPAND( g_mom1 += d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*dvol;  ,
                g_mom2 += d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*dvol;  ,
                g_mom3 += d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*dvol; )
        TkeV = (d->Vc[PRS][k][j][IBEG]*CONST_mu*CONST_mp*UNIT_VELOCITY*UNIT_VELOCITY)/(d->Vc[RHO][k][j][IBEG]*CONST_kB*1.16e9);
        if (TkeV >= 0.5 && TkeV <=8.0){
          g_mhot += d->Vc[RHO][k][j][i]*dvol;

        } else if (TkeV < 0.5){
           g_mcold += d->Vc[RHO][k][j][i]*dvol;
        }
  }
  #ifdef PARALLEL
  if (grid[0].rank_coord == 0 ){
  #endif
  KDOM_LOOP(k){
  JDOM_LOOP(j){
  EXPAND( darea = 4.0*CONST_PI*x1[IBEG]*x1[IBEG];,
          darea = 2.0*CONST_PI*x1[IBEG]*x1[IBEG]*sin(x2[j])*dx2[j];,
          darea = x1[IBEG]*x1[IBEG]*sin(x2[j])*dx2[j]*dx3[k];)

  TkeV = (d->Vc[PRS][k][j][IBEG]*CONST_mu*CONST_mp*UNIT_VELOCITY*UNIT_VELOCITY)/(d->Vc[RHO][k][j][IBEG]*CONST_kB*1.16e9);

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
   sendarray[14] = msub; sendarray[15] = g_mcold; sendarray[16] = g_mhot;
   MPI_Reduce (sendarray, recvarray, 17, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
   if (prank == 0){
     g_mass=recvarray[0]; g_TE=recvarray[1]; g_KE1=recvarray[2]; g_KE2=recvarray[3];
     g_KE3=recvarray[4]; g_mom1=recvarray[5]; g_mom2=recvarray[6]; g_mom3=recvarray[7];
     g_mdot=recvarray[8]; mdot_x = recvarray[9]; mdot_c = recvarray[10];
     madd = recvarray[11]; eadd = recvarray[12]; mdot_jet = recvarray[13];
     msub = recvarray[14]; g_mcold = recvarray[15]; g_mhot = recvarray[16];
  #endif
    fprintf(hist_file,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n ", g_time, g_dt, g_mass, g_TE, g_KE1, g_KE2, g_KE3, g_mom1, g_mom2, g_mom3, g_mdot, mdot_x, mdot_c, madd, eadd, mdot_jet, msub, g_mcold, g_mhot);

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
  double dend, pend;
  double  *x1, *x2, *x3, *dx1, *dx2, *dx3;
  double DM_POT[NX1_TOT];
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;
  dx1 = grid[IDIR].dx;
  dx2 = grid[JDIR].dx;
  dx3 = grid[KDIR].dx;

   if (g_stepNumber==0){
     Get_ICMDP( x1[IEND], &dend, &pend); //fixing the outer d, p to the smooth values 
     g_d_outer = d->Vc[RHO][0][0][IEND];
     g_p_outer = d->Vc[PRS][0][0][IEND];
     ITOT_LOOP(i){
     DM_POT[i] = BodyForcePotential(x1[i], 0.0, 0.0);
     }
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

double ran2(long *idum)

/*!
 *  *  Random number generator taken from Numerical Recipes in C.
 *   *  Returns a uniform deviate between 0.0 and 1.0.
 *    *  Call with idum a negative integer to initialize; 
 *     *  thereafter do not alter idum between successive deviates in a sequence.
 *      *  RNMX should approximate the largest floating value that is less than 1.
 *       ************************************************************************* */
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
