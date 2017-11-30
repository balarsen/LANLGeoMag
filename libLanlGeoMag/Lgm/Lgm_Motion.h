#ifndef LGM_MOTION_H
#define LGM_MOTION_H
#include <math.h>
#include "Lgm_Constants.h"
// #define LGM_c               (2.99792458e8)          // Speed of light  m/s
// #define LGM_Ee0            (0.510998910)                // Electron rest energy in MeV
// #define LGM_Ep0             (938.27201323)              // Proton rest energy in MeV

double Lgm_getVelocity(double energy, double mass);
  /***********************
   * energy in eV
   * mass in MeV/c**2
   *
   * get the veocity of the particle in m/s
   * KE = mc**2*(1/gamma -1), gamma = sqrt(1-v**2/c**2)
   * Solve this for v gives
   * v = c*sqrt(1-(mc**2/E+1)**2)
   ***********************/

double Lgm_getVelocity_nonrel(double energy, double mass); 
  /***********************
   * energy in eV
   * mass in MeV/c**2
   *
   * non relatavistic calculation
   * E = 1/2 m v**2
   * v = sqrt(2E/m)
   ***********************/


#endif
