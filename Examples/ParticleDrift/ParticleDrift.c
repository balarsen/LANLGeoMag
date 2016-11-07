

/********************************************************************
 *  Particle Drift
 *  Brian Larsen
 *  Nov 6, 2016
 *
 * This code takes in a position and pitch angle and tracks that 
 * particle around the globe using drift averaged quantities. 
 *
 * Input
 * Position: GSM input position (Re)
 * datetime: date and time of the calcualtion
 * Pitch angle: pitch angle at the postion (deg)
 * B model: the magneti field model to use (T89 to start)
 * Number of drifts to calculate: (start with 1)
 * Output filename: JSON headed ASCII filename to write
 *
 *
 * Algorithm
 * - given position, calculate L, MLT (done)
 * - calculate mirror point positions (done)
 * - calculate magnetic equator position (done)
 * - trace a line out from the south mirror point to north atmosphere (done)
 * - trim the trace at the north mirror point (done)
 * - compute the drift velocity at each point along the trace (done)
 * - create and write out a JSON headed file with all the data (done)
 * - move the particle to a new position and repeat
 *     - how to decide where to move it to? Maybe at mag equator...
 ********************************************************************/


#include <stdio.h>
#include <string.h>
#include <math.h>

#include <Lgm_Vec.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_Metadata.h>

// in Coulombs
const double E_CHARGE=1.6021766208e-19;
// MeV/c**2
const double E_RESTMASS=0.5109989461;
//  MeV/c**2
const double P_RESTMASS=938.2720813;


void printTraceResult(int traceStatus, const char *msg ) {
  switch (traceStatus) {
  case LGM_OPEN_IMF :
    printf("%s is LGM_OPEN_IMF\n", msg);
    exit(-1);
  case LGM_CLOSED :
    printf("%s is LGM_CLOSED\n", msg);
    break;
  case LGM_OPEN_N_LOBE :
    printf("%s is LGM_OPEN_N_LOBE\n", msg);
    exit(-2);
  case LGM_OPEN_S_LOBE :
    printf("%s is LGM_OPEN_S_LOBE\n", msg);
    exit(-3);
  case LGM_INSIDE_EARTH :
    printf("%s is LGM_INSIDE_EARTH\n", msg);
    exit(-4);
  case LGM_TARGET_HEIGHT_UNREACHABLE :
    printf("%s is LGM_TARGET_HEIGHT_UNREACHABLE\n", msg);
    exit(-5);
  case LGM_BAD_TRACE :
    printf("%s is LGM_BAD_TRACE\n", msg);
    exit(-6);
  default :
    printf("%s trace was undefined\n", msg);
    exit(-7);
  }
  return;
}



int main(void) {
  long int            Date;
  double              UTC, initA, initL, I, Bm, M, Sm;
  Lgm_Vector          initPos, initMirrorN, initMirrorS, initBmin, lineMirrorN, tmpV;
  // Lgm_Vector          lineMirrorS;
  Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
  Lgm_CTrans          *c = Lgm_init_ctrans( 0 );
  Lgm_ElapsedTimeInfo  tInfo;
  int                 traceStatus;
  double              *sN;
  //  double              *sS;
  double              *BmagN;
  /* double              *BmagS; */
  double              *PxN, *PyN, *PzN;
  /* double              *PxS, *PyS, *PzS; */
  FILE                *fp_out;
  int                 nPtsN;
  /* int                 nPtsS; */
  const double        trace_cut_tol=0.01; // TODO should this be passed in
  double              energy; // MeV
  Lgm_Vector          eVel, pVel;
  Lgm_Vector          *eVelarr, *pVelarr;

  
  energy = 0.050; // TODO should be input, 50keV
  
  Lgm_InitMagInfoDefaults(mInfo);
  //  MagInfo->Lgm_TraceToMirrorPoint_Tol
  
  Lgm_ElapsedTimeInit( &tInfo, 255, 95, 0 );

  // TODO input
  Date = 20100101;
  // TODO input
  UTC = 12.4;
  // TODO input
  initA = 80;
  // TODO input
  Lgm_SetVecElements(&initPos, 5.5, 1, 0.3); // GSM
  // TODO input
  
  
  Lgm_Set_Lgm_B_T89c(mInfo);

  initL = Lgm_McIlwain_L( Date, UTC, &initPos, initA, 0, &I, &Bm, &M, mInfo );

  printf("Pitch Angle: %g    McIlwain L  = %.10g   ( I, Bm, M = %g %g %g )\n", initA, initL, I, Bm, M);

  // int  Lgm_TraceToMirrorPoint( Lgm_Vector *u, Lgm_Vector *v, double *Sm, double Bm, double sgn, double tol, Lgm_MagModelInfo *Info );

  printf("Initial point is: ");
  Lgm_PrintVector(&initPos);


  //  printf("MagInfo->Lgm_TraceToMirrorPoint_Tol = %lf   %lf \n", mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo->P);

  // TODO, Bm is set to 1000, reasonable?
  traceStatus = Lgm_TraceToMirrorPoint(&initPos, &initMirrorN, &Sm, 1000, 1, mInfo->Lgm_TraceToMirrorPoint_Tol,
				       mInfo);
  printTraceResult( traceStatus, "North mirror");
  printf("North mirror point is: ");
  Lgm_PrintVector(&initMirrorN);

  traceStatus = Lgm_TraceToMirrorPoint(&initPos, &initMirrorS, &Sm, 1000, -1, mInfo->Lgm_TraceToMirrorPoint_Tol,
				       mInfo);
  printTraceResult( traceStatus, "South mirror");
  printf("South mirror point is: ");
  Lgm_PrintVector(&initMirrorS);

  // TODO are SM or the 1000 good values?
  // int  Lgm_TraceToMinBSurf( Lgm_Vector *, Lgm_Vector *, double, double, Lgm_MagModelInfo * );
  traceStatus = Lgm_TraceToMinBSurf(&initPos, &initBmin,  Sm, 1000, mInfo);
  printTraceResult( traceStatus, "Equator");
  printf("Equator point is: ");
  Lgm_PrintVector(&initBmin);

  // trace out the line north
  //  int Lgm_TraceLine( Lgm_Vector *u, Lgm_Vector *v, double H0, double sgn, double tol, int AddBminPoint, Lgm_MagModelInfo *Info ) {
  // TODO H0 is km above earth to trace to, 100, this could be a parameter
  traceStatus = Lgm_TraceLine(&initMirrorS, &lineMirrorN, 100, 1, 1e-5, 1, mInfo);
  nPtsN = mInfo->nPnts;
  /* double              s[LGM_MAX_INTERP_PNTS];    // distance along FL */
  /* double              Px[LGM_MAX_INTERP_PNTS];   // Px along FL  (in GSM) */
  /* double              Py[LGM_MAX_INTERP_PNTS];   // Py along FL  (in GSM) */
  /* double              Pz[LGM_MAX_INTERP_PNTS];   // Pz along FL  (in GSM) */
  /* Lgm_Vector          Bvec[LGM_MAX_INTERP_PNTS]; // 3D B-field vector   (in GSM) */
  /* double              Bmag[LGM_MAX_INTERP_PNTS]; // magnitude of B */
  /* double              BminusBcdip[LGM_MAX_INTERP_PNTS]; // magnitude of B minus magnitude of Cent. Dipole */
  /* double              MaxDiv;     // Dont subdivide the FL length with steps bigger than this. */
  /* int                 nDivs;      // Number of divisions of FL length to try to make (actual number of points defined may be different as MAxDiv mux be respected.) */
  /* int                 nPnts;      // actual number of points defined */

  // copy the array from mInfo to an array and store it
  /* void *memcpy(void *str1, const void *str2, size_t n) */
  /*   str1 -- This is pointer to the destination array where the content is to be copied, type-casted to a pointer of type void*. */
  /*   str2 -- This is pointer to the source of data to be copied, type-casted to a pointer of type void*. */
  /*   n -- This is the number of bytes to be copied. */
  sN = calloc(nPtsN, sizeof (double));
  PxN = calloc(nPtsN, sizeof (double));
  PyN = calloc(nPtsN, sizeof (double));
  PzN = calloc(nPtsN, sizeof (double));
  BmagN = calloc(nPtsN, sizeof (double));
  eVelarr = calloc(nPtsN, sizeof (Lgm_Vector));
  pVelarr = calloc(nPtsN, sizeof (Lgm_Vector));
  memcpy(sN, mInfo->s, nPtsN*sizeof(double));
  memcpy(PxN, mInfo->Px, nPtsN*sizeof(double));
  memcpy(PyN, mInfo->Py, nPtsN*sizeof(double));
  memcpy(PzN, mInfo->Pz, nPtsN*sizeof(double));
  memcpy(BmagN, mInfo->Bmag, nPtsN*sizeof(double));
  

  printTraceResult( traceStatus, "Equator");
  printf("Traced from: ");
  Lgm_PrintVector(&initMirrorS);  
  printf("Traced to: ");
  Lgm_PrintVector(&lineMirrorN);
  printf("There are %d points along the line\n", nPtsN);

  printf("Trimming the trace off at the North mirror point\n");
  for (int i=0; i<nPtsN; i++) {
    // how do we know we have passed initMirrorN?
    // TODO assumimng that we can just look at GSM-Z and stop then
    if (  fabs(PzN[i]-initMirrorN.z)<trace_cut_tol &&
	  fabs(PxN[i]-initMirrorN.x)<trace_cut_tol &&
	  fabs(PyN[i]-initMirrorN.y)<trace_cut_tol ) {
      printf("Found the north cutoff at: ");
      Lgm_SetVecElements(&tmpV, PxN[i], PyN[i], PzN[i]);
      Lgm_PrintVector(&tmpV);
      printf("which should look like: ");
      Lgm_PrintVector(&initMirrorN);
      nPtsN = i+1;
      break;
    }
  }
  printf("This leaves %d points along the line\n", nPtsN);


  /*************************/
  // compute the velocity at each point
  /*      \param[in]      u0          Position (in GSM) to use. [Re]
   *      \param[out]     Vel         The computed drift velocity in GSM coords. [km/s]
   *      \param[in]      q           Charge of the particle. [C]
   *      \param[in]      T           Kinetic energy of the particle. [MeV]
   *      \param[in]      E0          Rest energy of the particle. [MeV]
   *      \param[in]      Bm          Magnitude of B-field at mirror point for particle. [nT]
   *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
   *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
   *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
   *
   *
   *      \author         Mike Henderson
   *      \date           2011
   *
   */
  //int Lgm_GradAndCurvDriftVel( Lgm_Vector *u0, Lgm_Vector *Vel, Lgm_MagModelInfo *m ) {

  {
    double km_s;
    Lgm_Vector veltmp;

    printf("Computing drift velocity at input point\n");
    // compute electron drift velocity
    mInfo->Lgm_VelStep_E0 = E_RESTMASS;
    mInfo->Lgm_VelStep_q = -1*E_CHARGE;
    mInfo->Lgm_VelStep_T = energy;
    mInfo->Lgm_VelStep_Bm = BmagN[0]; // Magnitude of B-field at mirror point for particle. [nT]
    mInfo->Lgm_VelStep_h = 0.01; //TODO should this be input
    mInfo->Lgm_VelStep_DerivScheme = LGM_DERIV_SIX_POINT;
    Lgm_GradAndCurvDriftVel(&initPos, &eVel, mInfo );
    printf("electron drift velocity has been computed and is: ");
    Lgm_PrintVector(&eVel);
    // also print in km/s
    km_s = Lgm_Magnitude(&eVel)*Re;
    printf("   which is %lf km/s\n", km_s);
    
    // compute proton drift velocity
    mInfo->Lgm_VelStep_E0 = P_RESTMASS;
    mInfo->Lgm_VelStep_q = E_CHARGE;
    mInfo->Lgm_VelStep_T = energy;
    mInfo->Lgm_VelStep_Bm = BmagN[0]; // Magnitude of B-field at mirror point for particle. [nT]
    mInfo->Lgm_VelStep_h = 0.01; //TODO should this be input
    mInfo->Lgm_VelStep_DerivScheme = LGM_DERIV_SIX_POINT;
    Lgm_GradAndCurvDriftVel(&initPos, &pVel, mInfo );
    printf("proton drift velocity has been computed and is: ");
    Lgm_PrintVector(&pVel);
    km_s = Lgm_Magnitude(&pVel)*Re;
    printf("   which is %lf km/s\n", km_s);

    for (int i=0; i<nPtsN; i++) {
      mInfo->Lgm_VelStep_E0 = E_RESTMASS;
      mInfo->Lgm_VelStep_q = -1*E_CHARGE;
      mInfo->Lgm_VelStep_T = energy;
      mInfo->Lgm_VelStep_Bm = BmagN[0]; // Magnitude of B-field at mirror point for particle. [nT]
      mInfo->Lgm_VelStep_h = 0.01; //TODO should this be input
      mInfo->Lgm_VelStep_DerivScheme = LGM_DERIV_SIX_POINT;
      Lgm_SetVecElements(&veltmp, PxN[i], PyN[i], PzN[i]);
      Lgm_GradAndCurvDriftVel(&veltmp, &eVelarr[i], mInfo );
    }
    for (int i=0; i<nPtsN; i++) {
      mInfo->Lgm_VelStep_E0 = P_RESTMASS;
      mInfo->Lgm_VelStep_q = E_CHARGE;
      mInfo->Lgm_VelStep_T = energy;
      mInfo->Lgm_VelStep_Bm = BmagN[0]; // Magnitude of B-field at mirror point for particle. [nT]
      mInfo->Lgm_VelStep_h = 0.01; //TODO should this be input
      mInfo->Lgm_VelStep_DerivScheme = LGM_DERIV_SIX_POINT;
      Lgm_SetVecElements(&veltmp, PxN[i], PyN[i], PzN[i]);
      Lgm_GradAndCurvDriftVel(&veltmp, &pVelarr[i], mInfo );
    }

    
  }


  

  

  {
    /* double *eVelx, *eVely, *eVelz; */
    /* double *pVelx, *pVely, *pVelz; */

    /* eVelx = calloc(nPtsN, sizeof (double)); */
    /* eVely = calloc(nPtsN, sizeof (double)); */
    /* eVelz = calloc(nPtsN, sizeof (double)); */
    
    /* pVelx = calloc(nPtsN, sizeof (double)); */
    /* pVely = calloc(nPtsN, sizeof (double)); */
    /* pVelz = calloc(nPtsN, sizeof (double)); */
    
    /* need a meatadata variable for each entry */
    Lgm_metadata_variable sN_meta;
    Lgm_metadata_variable PxN_meta;
    Lgm_metadata_variable PyN_meta;
    Lgm_metadata_variable PzN_meta;
    Lgm_metadata_variable BmagN_meta;

    Lgm_metadata_variable eVelx_meta;
    Lgm_metadata_variable eVely_meta;
    Lgm_metadata_variable eVelz_meta;

    Lgm_metadata_variable pVelx_meta;
    Lgm_metadata_variable pVely_meta;
    Lgm_metadata_variable pVelz_meta;

    
    /* Lgm_metadata_variable meta_meta_meta; */
    // and a string to put the answer
    char* json_header;

    // create a 1 dimensional variable, with data, a name, and the variable
    Lgm_metadata_initvar(&sN_meta, // the variable to store it
			 1,            // Dimension
			 1,            // has data (bool)
			 "sN"      // Name 
			 );
    // add an attribute to the variable (all arguments must be strings)
    Lgm_metadata_addStringAttr(&sN_meta, // the variable to add it to
			       "UNITS",      // name
			       "unknown",        // value
			       0             // boolean it is not an array 
			       );
    Lgm_metadata_addStringAttr(&sN_meta, "LABEL", "s", 0);
    Lgm_metadata_addStringAttr(&sN_meta, "COMMENT", "Distance along the file line from South to North mirror point", 0);

    Lgm_metadata_initvar(&PxN_meta, 1, 1, "X");
    Lgm_metadata_addStringAttr(&PxN_meta, "Units", "Re", 0);
    Lgm_metadata_addStringAttr(&PxN_meta, "LABEL", "X", 0);
    Lgm_metadata_addStringAttr(&PxN_meta, "COMMENT", "x position of the field line going north", 0);

    Lgm_metadata_initvar(&PyN_meta, 1, 1, "Y");
    Lgm_metadata_addStringAttr(&PyN_meta, "Units", "Re", 0);
    Lgm_metadata_addStringAttr(&PyN_meta, "LABEL", "Y", 0);
    Lgm_metadata_addStringAttr(&PyN_meta, "COMMENT", "y position of the field line going north", 0);

    Lgm_metadata_initvar(&PzN_meta, 1, 1, "Z");
    Lgm_metadata_addStringAttr(&PzN_meta, "Units", "Re", 0);
    Lgm_metadata_addStringAttr(&PzN_meta, "LABEL", "Z", 0);
    Lgm_metadata_addStringAttr(&PzN_meta, "COMMENT", "z position of the field line going north", 0);

    Lgm_metadata_initvar(&BmagN_meta, 1, 1, "B");
    Lgm_metadata_addStringAttr(&BmagN_meta, "Units", "nT", 0);
    Lgm_metadata_addStringAttr(&BmagN_meta, "LABEL", "Bmag", 0);
    Lgm_metadata_addStringAttr(&BmagN_meta, "COMMENT", "Bmag at each point along the line", 0);

    Lgm_metadata_initvar(&eVelx_meta, 1, 1, "eVx");
    Lgm_metadata_addStringAttr(&eVelx_meta, "Units", "km/s", 0);
    Lgm_metadata_addStringAttr(&eVelx_meta, "LABEL", "eVx", 0);
    Lgm_metadata_addStringAttr(&eVelx_meta, "COMMENT", "Electron X velocity", 0);

    Lgm_metadata_initvar(&eVely_meta, 1, 1, "eVy");
    Lgm_metadata_addStringAttr(&eVely_meta, "Units", "km/s", 0);
    Lgm_metadata_addStringAttr(&eVely_meta, "LABEL", "eVy", 0);
    Lgm_metadata_addStringAttr(&eVely_meta, "COMMENT", "Electron Y velocity", 0);

    Lgm_metadata_initvar(&eVelz_meta, 1, 1, "eVz");
    Lgm_metadata_addStringAttr(&eVelz_meta, "Units", "km/s", 0);
    Lgm_metadata_addStringAttr(&eVelz_meta, "LABEL", "eVz", 0);
    Lgm_metadata_addStringAttr(&eVelz_meta, "COMMENT", "Electron Z velocity", 0);

    Lgm_metadata_initvar(&pVelx_meta, 1, 1, "pVx");
    Lgm_metadata_addStringAttr(&pVelx_meta, "Units", "km/s", 0);
    Lgm_metadata_addStringAttr(&pVelx_meta, "LABEL", "pVx", 0);
    Lgm_metadata_addStringAttr(&pVelx_meta, "COMMENT", "Proton X velocity", 0);

    Lgm_metadata_initvar(&pVely_meta, 1, 1, "pVy");
    Lgm_metadata_addStringAttr(&pVely_meta, "Units", "km/s", 0);
    Lgm_metadata_addStringAttr(&pVely_meta, "LABEL", "pVy", 0);
    Lgm_metadata_addStringAttr(&pVely_meta, "COMMENT", "Proton Y velocity", 0);

    Lgm_metadata_initvar(&pVelz_meta, 1, 1, "pVz");
    Lgm_metadata_addStringAttr(&pVelz_meta, "Units", "km/s", 0);
    Lgm_metadata_addStringAttr(&pVelz_meta, "LABEL", "pVz", 0);
    Lgm_metadata_addStringAttr(&pVelz_meta, "COMMENT", "Proton Z velocity", 0);
    
    // create the JSON header
    //   this is a variadic function, number of vars and all the vars
    json_header = Lgm_metadata_JSONheader(11, &sN_meta, &PxN_meta, &PyN_meta, &PzN_meta, &BmagN_meta,
					  &eVelx_meta, &eVely_meta, &eVelz_meta,
					  &pVelx_meta, &pVely_meta, &pVelz_meta);
    printf("%s", json_header);


    fp_out = fopen("file.txt", "w");
    if (fp_out == NULL)
      {
	printf("Error opening file!\n");
	exit(1);
      }
    fprintf(fp_out, "%s", json_header);
    for (int i=0; i<nPtsN; i++) {
      fprintf(fp_out, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", sN[i], PxN[i], PyN[i], PzN[i],
	      BmagN[i],
	      eVelarr[i].x*Re, eVelarr[i].y*Re, eVelarr[i].z*Re,
	      pVelarr[i].x*Re, pVelarr[i].y*Re, pVelarr[i].z*Re);
    }
    fclose(fp_out);
 

    

  }









  
  
  Lgm_PrintElapsedTime( &tInfo );
  Lgm_FreeMagInfo(mInfo);
  Lgm_free_ctrans(c);
  
  // TODO shold we free() or just let the exit occur?
  
  return(0);
}


