#define  PHYSICS                 RMHD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           HANCOCK
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 1
#define  USER_DEF_PARAMETERS     4
#define  USER_DEF_CONSTANTS      0

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  MHD_FORMULATION         CONSTRAINED_TRANSPORT

/* -- user-defined parameters (labels) -- */

#define  SIGMA_TOR               0
#define  SIGMA_POL               1
#define  VEL0                    2
#define  MACH                    3

/* -- user-defined symbolic constants -- */


/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          YES
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         NO
#define  SHOCK_FLATTENING          NO
#define  ARTIFICIAL_VISCOSITY      NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   MC_LIM
#define  CT_EMF_AVERAGE            ARITHMETIC
#define  CT_EN_CORRECTION          NO
#define  ASSIGN_VECTOR_POTENTIAL   YES
#define  UPDATE_VECTOR_POTENTIAL   NO
#define  PRIMITIVE_HANCOCK         NO
