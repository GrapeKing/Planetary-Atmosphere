#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  COMPONENTS                     1
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     (VECTOR+POTENTIAL)
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            4
#define  INTERNAL_BOUNDARY        	    YES

/* -- physics dependent declarations -- */

#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  RHOINF                         0
#define  BONDI                          1
#define  G                              2
#define  CS                             3

/* [Beg] user-defined constants (do not change this line) */


/* [End] user-defined constants (do not change this line) */
