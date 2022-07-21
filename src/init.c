/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
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
  double temp0 = g_inputParam[T_0];
  double mu    = g_inputParam[MU_MMW];

  g_gamma = g_inputParam[GAMMA];
  v[RHO]  = g_inputParam[RHO_0];
  v[VX1]  = g_inputParam[V_0];
  v[VX2]  = 0.0;
  v[VX3]  = 0.0;
  v[PRS]  = v[RHO]*temp0/(KELVIN*mu);
  v[TRC]  = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD

   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = g_inputParam[BP_0];

   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = 0.0;

  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
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
    int i;
    static int idx_1au;
    
    double r0,r1,dr,dr_test;
    double mu = g_inputParam[MU_MMW];
    double *r = grid->x[IDIR];

    char  fname[512];
    char fname0[512];
    FILE *fp, *fp0;
    
    char fname1[512];
    FILE *fp1;
    static int first_call = 1;
    double cme_start_time = g_inputParam[CME_START_TIME];
    double cme_duration = g_inputParam[CME_DURATION];
    double cme_end_time;
    static double r_x0[2];
    double r_x1[2], v_x0[2], h[2];
    double *** vr = d ->Vc[VX1];
    int ix0, ix1;
    int idx0[2], idx1[2];
    
    sprintf (fname , "%s/obs_1au.dat",RuntimeGet()->output_dir);
    sprintf (fname0,  "%s/obs_r0.dat",RuntimeGet()->output_dir);    
    sprintf (fname1,  "%s/tracers.dat",RuntimeGet()->output_dir);
    
    if (g_stepNumber == 0){ /* Open file for writing and compute idxat start only */

      /* ---- Get grid point index closest to 1 AU ---- */

      idx_1au = 0;
      dr=1e99;
      for (i=0; i<= NX1_TOT; i++){
        dr_test = fabs(r[i]-1.0);
        if (dr_test<dr){
          dr=dr_test;
          idx_1au=i;
        }
      }
      fp = fopen(fname,"w");
      fprintf (fp,"Observer point at r[%d]=%12.6e AU\n",idx_1au,r[idx_1au]);
      fprintf (fp,"%7s %12s %12s %12s %12s %12s\n","time","vr","rho","temp","Bp","prs");

      fp0 = fopen(fname0,"w");
      fprintf (fp0,"Observer point at r[%d]=%12.6e AU\n",0,r[0]);
      fprintf (fp0,"%7s %12s %12s %12s %12s %12s\n","time","vr","rho","temp","Bp","prs");      
      
    }else{
      fp  = fopen(fname ,"a");
      fp0 = fopen(fname0,"a");
    }

    fprintf (fp,"%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",g_time,
        d->Vc[VX1][0][0][idx_1au],
        d->Vc[RHO][0][0][idx_1au],
        (KELVIN*mu)*(d->Vc[PRS][0][0][idx_1au])/(d->Vc[RHO][0][0][idx_1au]),
	     d->Vc[BX3][0][0][idx_1au],
	     d->Vc[PRS][0][0][idx_1au]);   
    fclose(fp);

    fprintf (fp0,"%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",g_time,
        d->Vc[VX1][0][0][0],
        d->Vc[RHO][0][0][0],
        (KELVIN*mu)*(d->Vc[PRS][0][0][0])/(d->Vc[RHO][0][0][0]),
	     d->Vc[BX3][0][0][0],
	     d->Vc[PRS][0][0][0]);    
    fclose(fp0);
    

/*  Initial position for tracers is at inner boundary. Start advecting at the beginning/end of CME */
    if (first_call) {
        fp1 = fopen(fname1,"w");

        r_x0[0] = 0.1395998; /* 30 Rs */
        r_x0[1] = 0.1395998;
        
        fprintf (fp1,"%7s %12s %12s\n","time","tracer1","tracer2");
        fclose(fp1);
        first_call = 0;
        
    }
    if (g_time >= cme_start_time) {
        /* Locate the tracers on the grid, needed for the velocity interpolation */
        for (i=0; i <= NX1_TOT; i++) {
            if (r_x0[0] > r[i] && r_x0[0] < r[(i+1)]) {
                idx0[0] = i-1;
                idx1[0] = i;
                idx0[0] = MAX(0,idx0[0]);
                idx0[0] = MIN(NX1_TOT-1, idx0[0]);
                idx1[0] = MIN(NX1_TOT-1, idx1[0]);
            }
            if (r_x0[1] > r[i] && r_x0[1] < r[(i+1)]) {
                idx0[1] = i-1;
                idx1[1] = i;
                idx0[1] = MAX(0,idx0[1]);
                idx0[1] = MIN(NX1_TOT-1, idx0[1]);
                idx1[1] = MIN(NX1_TOT-1, idx1[1]);
            }
        }
        cme_end_time = cme_start_time + cme_duration;
        
        h[0] = (r_x0[0] - r[idx0[0]]) / (r[idx1[0]] - r[idx0[0]]);
	h[1] = (r_x0[1] - r[idx0[1]]) / (r[idx1[1]] - r[idx0[1]]);

	/* Interpolate velocity */
	
        v_x0[0] = (1.0 - h[0]) * vr[0][0][idx0[0]] + h[0] * vr[0][0][idx1[0]];
        v_x0[1] = (1.0 - h[1]) * vr[0][0][idx0[1]] + h[1] * vr[0][0][idx1[1]];

	
        if (g_time < cme_end_time) {
            v_x0[1] = 0.0;
        }

	/* advect */
	
        r_x1[0] = r_x0[0] + v_x0[0] * g_dt;
        r_x1[1] = r_x0[1] + v_x0[1] * g_dt;
        
        /* update location */
	
        r_x0[0] = r_x1[0];
        r_x0[1] = r_x1[1];

	/* ensure that tracers do not leave the computational domain */
	
        r_x0[0] = MIN(r_x0[0], r[(NX1_TOT-1)]);
	r_x0[1] = MIN(r_x0[1], r[(NX1_TOT-1)]);

        fp1  = fopen(fname1 ,"a");
        fprintf (fp1,"%12.6e %12.6e %12.6e\n",g_time, r_x0[0], r_x0[1]);
        fclose(fp1);
                 
    }
    
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

  double time_cme_stop, fact, factB;
  double time_cme_start = g_inputParam[CME_START_TIME];
  double dt_cme         = g_inputParam[CME_DURATION];
  double v0             = g_inputParam[V_0];
  double v_pert         = g_inputParam[V_PERT];
  double rho0           = g_inputParam[RHO_0];
  double rho_pert       = g_inputParam[RHO_PERT];
  double bp0            = g_inputParam[BP_0];
  double bp_pert        = g_inputParam[BP_PERT];
  double temp0          = g_inputParam[T_0];
  double mu             = g_inputParam[MU_MMW];
  int    profile        = g_inputParam[PROFILE];

  time_cme_stop  = time_cme_start + dt_cme;

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
        if ( (g_time < time_cme_start) || (g_time > time_cme_stop) ){
          X1_BEG_LOOP(k,j,i) {
            d->Vc[RHO][k][j][i] = rho0; // Density at r=r0
            d->Vc[VX1][k][j][i] = v0;   // Speed at r=r0 km/s
            d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*temp0/(KELVIN*mu);   // Thermal pressure at r=r0
            d->Vc[BX3][k][j][i] = bp0; // B-phi at r=r0
          }
        }else{

	  fact = pow(sin((g_time-time_cme_start)*CONST_PI/dt_cme),2);
	  if (profile == 0) {
	    factB = fact; 
	  } else {
	    factB = sin((g_time-time_cme_start)*2*CONST_PI/dt_cme);
	  } 
          /*fact = pow(sin((g_time-time_cme_start)*CONST_PI/dt_cme),2);
	    factB = sin((g_time-time_cme_start)*2*CONST_PI/dt_cme); */
          X1_BEG_LOOP(k,j,i) {
            d->Vc[RHO][k][j][i] = rho0 + fact*rho_pert;  // Density at r=r0
            d->Vc[VX1][k][j][i] = v0   + fact*v_pert;    // Velocity at r=r0
            d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*temp0/(KELVIN*mu); // Thermal pressure at r=r0
            d->Vc[BX3][k][j][i] = bp0  + factB*bp_pert; //B-phi at r=r0
	  }
	}
  }
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
  return 0.0;
}
#endif
