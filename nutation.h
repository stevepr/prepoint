/*
   Naval Observatory Vector Astrometry Software (NOVAS)
   C Edition, Version 3.1
   
   nutation.h: Header file for nutation models

   U. S. Naval Observatory
   Astronomical Applications Dept.
   Washington, DC 
   http://www.usno.navy.mil/USNO/astronomical-applications
*/

#ifndef _NUTATION_
   #define _NUTATION_

// sbp ... for handling calls from c++ main routine
//
#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/*
   Function prototypes
*/

    EXTERNC void iau2000a (double jd_high, double jd_low,

                  double *dpsi, double *deps);

    EXTERNC void iau2000b (double jd_high, double jd_low,

                  double *dpsi, double *deps);

    EXTERNC void nu2000k (double jd_high, double jd_low,

                 double *dpsi, double *deps);

#endif
