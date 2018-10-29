#define  PHYSICS                 MHD
#define  DIMENSIONS              1
#define  COMPONENTS              1
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 TABULATED
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 1
#define  USER_DEF_PARAMETERS     17
#define  USER_DEF_CONSTANTS      3

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  MHD_FORMULATION         EIGHT_WAVES
#define  BACKGROUND_FIELD        NO
#define  RESISTIVE_MHD           NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  K0                      0
#define  K100                    1
#define  ne_out                  2
#define  kalpha                  3
#define  amp                     4
#define  EFF                     5
#define  SIGR                    6
#define  SIGTH                   7
#define  RJETKPC                 8
#define  THJET                   9
#define  VJET                    10
#define  M_smbh                  11
#define  c200                    12
#define  M200                    13
#define  Vc_sis                  14
#define  r0_sis                  15
#define  RMAX                    16

/* -- user-defined symbolic constants -- */

#define  UNIT_DENSITY            (CONST_mp)
#define  UNIT_LENGTH             (CONST_pc*1.e5)
#define  UNIT_VELOCITY           (1.e8)

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             NO
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          MULTID
#define  ARTIFICIAL_VISCOSITY      NO
#define  CHAR_LIMITING             YES
#define  LIMITER                   DEFAULT
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  UPDATE_VECTOR_POTENTIAL   NO