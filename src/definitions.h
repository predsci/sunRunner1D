#define  PHYSICS                        MHD
#define  DIMENSIONS                     1
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LimO3
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            12

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   EIGHT_WAVES
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  GAMMA                          0
#define  MU_MMW                         1
#define  T_0                            2
#define  RHO_0                          3
#define  RHO_PERT                       4
#define  V_0                            5
#define  V_PERT                         6
#define  BP_0                           7
#define  BP_PERT                        8
#define  CME_START_TIME                 9
#define  CME_DURATION                   10
#define  PROFILE                        11

/* [Beg] user-defined constants (do not change this line) */

#define  SHOCK_FLATTENING               ONED
#define  CHAR_LIMITING                  YES
#define  LIMITER                        MINMOD_LIM

/* [End] user-defined constants (do not change this line) */
