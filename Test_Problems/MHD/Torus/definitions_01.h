#define  PHYSICS                 MHD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                CYLINDRICAL
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 1
#define  USER_DEF_PARAMETERS     6
#define  USER_DEF_CONSTANTS      1

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  MHD_FORMULATION         CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        NO
#define  RESISTIVE_MHD           NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  RMIN                    0
#define  RMAX                    1
#define  RHO_CUT                 2
#define  BETA                    3
#define  ETA                     4
#define  SCALE_HEIGHT            5

/* -- user-defined symbolic constants -- */

#define  USE_DIPOLE              NO

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          YES
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         NO
#define  SHOCK_FLATTENING          NO
#define  ARTIFICIAL_VISCOSITY      NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   MC_LIM
#define  CT_EMF_AVERAGE            UCT_HLL
#define  CT_EN_CORRECTION          NO
#define  ASSIGN_VECTOR_POTENTIAL   YES
#define  UPDATE_VECTOR_POTENTIAL   NO
