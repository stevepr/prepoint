// prepoint application
//	Generates a list of PrePointing ICRS RA,DE positions for a desired Target ICRS RA,DE position and Target Time
//		When non-moving scope is pointed at a given prepoint location at the prepoint time, the target position will
//		drift into the center of the FOV at the specified target time.
//	
//	Usage: prepoint --ra=ddd.ddd|hh:mm:ss.s --de=ddd.ddd|+dd:mm:ss.s --time==jjjj.jjj|yyyy-mm-dd_hh:mm:ss [--hours=h.hhhh | --days=d.dddd] [--interval=sec] 
//						[--start=jjjj.jjj|yyyy-mm-dd_hh:mm:ss]
//			ra,de = ICRS target position (deg or hms,dms)
//			time = Target acquisition time (UTC)
//			hours or days = duration of time prior to the target time for which prepoint times are generated
//			interval = interval between each prepoint time in the output list of prepoint times & positions (sec)
//			start = latest prepoint time in the output list
// 
//	OUTPUT
//		Header lines and a CSV table is sent to the standard output.  TETE positions are essentially JNOW for mounts.
// 
//		example:
//			Target Date (UTC): 2460629.7708333 JD = 2024-11-15 06:30:00
//			Target GCRS position: 05:26:17.50, +28:36:28.30
//			Target TETE position: 05:27:51.67, +28:37:48.33
//			Prepoint DateTime, GCRS RA, GCRS DE, TETE RA, TETE DE
//			2024-11-15 06:29:50, 05:26:07.48, +28:36:27.92, 05:27:41.65, +28:37:48.33
// 
//	NOTES
//		This program uses the NOVAS library for astrometric computations as developed by
//		the Astronomical Applications Department of the U.S. Naval Observatory
//

#include <iostream>
#include <string.h>
#include <math.h>
#include "novas.h"
#include "novascon.h"
#include "nutation.h"


#define USAGE   "Usage: prepoint --ra=ddd.ddd|hh:mm:ss.s --de=ddd.ddd|+dd:mm:ss.s\n"\
				"			--time==jjjj.jjj|yyyy-mm-dd_hh:mm:ss [--hours=h.hhhh | --days=d.dddd]\n"\
				"			[--interval=sec] [--start=jjjj.jjj|yyyy-mm-dd_hh:mm:ss]\n"\
				"	ra,de = ICRS target position (deg or hms,dms)\n"\
				"	time = target acquistion time (UTC)\n"\
				"	hours or days = duration of time prior to the target time\n"\
				"	                for which prepoint times are generated\n"\
				"	interval = interval between each prepoint time in the output list\n"\
				"	           of prepoint times & positions (sec)\n"\
				"	start = latest prepoint time in the output list.\n"

#define EXIT_FAILURE	1
#define cLF "\n"
#define cCRLF "\r\n"
#define	C_IMAXLINE		1000					/* Max chars in an input line */

#define MJD_OFFSET 2400000.5

// function prototypes
//
int gcrs2tete(double jd_tdb, double ra_gcrs, double de_gcrs,
	double* ra_tete, double* de_tete);

int tete2gcrs(double jd_tdb, double ra_tete, double de_tete,
	double* ra_gcrs, double* de_gcrs);

void tete2altaz(double jd_ut1, double delta_t, double ra_tete, double de_tete, double obs_lat, double obs_long, double obs_h,
	double* alt, double* az);

void d2dms(double de, char *sigh, int* deg, int* min, double* sec);
void r2hms(double ra, int* hrs, int* min, double* sec);
int parseRA(char* arg, double* ra);
int parseDE(char* arg, double* de);
int parseDateTime(char* arg, double* jd);


//**********************************************
//
//   MAIN
//
//**********************************************


int main(int argc, char* argv[])
{
	char* arg;
	double dblTmp;

	// Inputs
	//
	double ra_target = -100.0;			// target RA in deg (ICRS)
	double de_target = -100.0;			// target DE in deg (ICRS)
	double jd_target = -1.0;			// target time for prepoint location (UTC)
	double duration_hours = 24.0;		// duration of prepoint times list (hours)
	double interval_sec = 10.0;			// interval between prepoint times (seconds)
	double deltaT = 69.5;
	double jd_start = -1.0;				// < 0 => start with first interval before target time

	// outputs
	//
	double ra_jnow;						// target RA deg - approx JNOW
	double de_jnow;						// target DE deg - approx JNOW
	short int yr, mon, day;
	double hr;
	int ih, im, isec;
	int dd, dm;
	double ds;
	int rh, rm;
	double rs;
	char sign;


	// Prepoint data
	//
	double H_target;			// Greenwich Hour Angle (hours)
	double jd_tdb;
	double mjd_prepoint;		// prepoint time UTC, MJD to avoid rounding errors due to double precision
	double mjd_end;
	double ra_gcrs;				// GCRS coordinate at prepoint time
	double de_gcrs;
	double ra_tete;				// TETE coordinate at prepoint time (approx jnow)
	double de_tete;
	double gast;				// GAST (hours)

	// misc
	//

    //************
    //  General init
    //


	/**********************************************
	 *	Get Command Line parameters
	 *		
	 */

	argc--;
	argv++;
	while (argc-- > 0) {
		arg = *argv++;

		// parse argument
		//
		if ((arg[0] != '-') || (arg[1] != '-'))
		{
			// no argument prefix
			//
			printf("argument syntax error\n");
			printf(USAGE);
			return(EXIT_FAILURE);
		}

		if (strncmp(arg + 2, "ra=", 3) == 0)				// --ra=
		{
			arg += 5;		// first char of numeric value
			if (parseRA(arg,&ra_target) != 0)
			{
				printf("error parsing ra= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}
		else if (strncmp(arg + 2, "de=", 3) == 0)			// --de=
		{
			arg += 5;		// first char of numeric value
			if (parseDE(arg,&de_target) != 0)
			{
				printf("error parsing de= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}
		else if (strncmp(arg + 2, "time=", 5) == 0)			// --time=
		{
			arg += 7;		// first char of numeric value
			if (parseDateTime(arg,&jd_target) != 0)
			{
				printf("error parsing time= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}
		else if (strncmp(arg + 2, "hours=", 6) == 0)		// --hours=
		{
			arg += 8;		// first char of numeric value
			if (sscanf(arg, "%lf", &dblTmp) == 0)
			{
				printf("error parsing hours= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
			duration_hours = dblTmp;
		}
		else if (strncmp(arg + 2, "days=", 5) == 0)			// --days=
		{
			arg += 7;		// first char of numeric value
			if (sscanf(arg, "%lf", &dblTmp) == 0)
			{
				printf("error parsing days= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
			duration_hours = dblTmp * 24.0;		// days to hours
		}
		else if (strncmp(arg + 2, "interval=", 9) == 0)		// --interval=
		{
			arg += 11;		// first char of numeric value
			if (sscanf(arg, "%lf", &interval_sec) == 0)
			{
				printf("error parsing interval= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}
		else if (strncmp(arg + 2, "start=", 6) == 0)		// --start=
		{
			arg += 8;		// first char of numeric value
			if (parseDateTime(arg, &jd_start) != 0)
			{
				printf("error parsing start= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}
		else if (strncmp(arg + 2, "deltaT=", 7) == 0)		// --deltaT=
		{
			arg += 9;		// first char of numeric value
			if (sscanf(arg, "%lf", &deltaT) == 0)
			{
				printf("error parsing deltaT= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}


	} // end of command line parsing

	// validate input arguments
	//

	if (ra_target < -90.0)
	{
		printf("argument error: no target RA specified.\n");
		printf(USAGE);
		return(EXIT_FAILURE);
	}
	if ((ra_target < 0.0) || (ra_target > 360.0))
	{
		printf("argument error: target RA is out of range.\n");
		printf(USAGE);
		return(EXIT_FAILURE);
	}

	if (de_target < -90.0)
	{
		printf("argument error: no target DE specified.\n");
		printf(USAGE);
		return(EXIT_FAILURE);
	}
	if ((de_target < -90.0) || (de_target > 90.0))
	{
		printf("argument error: target DE is out of range.\n");
		printf(USAGE);
		return(EXIT_FAILURE);
	}

	if (jd_target <= 0.0)
	{
		printf("argument error: no target date/time specified.\n");
		printf(USAGE);
		return(EXIT_FAILURE);
	}

	if (interval_sec < 0)
	{
		printf("argument error: interval < 0.\n");
		printf(USAGE);
		return(EXIT_FAILURE);
	}

	if (jd_start < 0.0)
	{
		// default -> start prepoint list with first interval before target time
		//
		jd_start = jd_target - interval_sec / 86400.0;
	}

	// DONE gather input parameters
	//

	/**********************************************
	 *	Compute prepoint positions
	 *		- compute TETE position of target at target time
	 *		- output TETE position at target time
	 *		- for each prepoint time
	 *			- output TETE and GCRS positions at prepoint time
	 *
	 */
	jd_tdb = jd_target + (deltaT/86400.0);

	if (gcrs2tete(jd_tdb, ra_target, de_target, &ra_jnow, &de_jnow) != 0)
	{
		printf("Error from gcrs2tete.");
		return (EXIT_FAILURE);
	}

	// output target data
	//	Target Date/Time
	//  Target GCRS coords
	//  Target TETE coords
	//
	cal_date(jd_target, &yr, &mon, &day, &hr);
	ih = (int)hr;
	dblTmp = hr - ih;
	dblTmp = dblTmp * 60.0;
	im = (int)dblTmp;
	dblTmp = dblTmp - im;
	isec = (int)(dblTmp * 60);
	printf("Target Date (UTC): %15.7lf JD = %04d-%02d-%02d %02d:%02d:%02d\n",jd_target, yr,mon,day,ih,im,isec);

	r2hms(ra_target, &rh, &rm, &rs);
	d2dms(de_target, &sign, &dd, &dm, &ds);
	printf("Target GCRS position: %02d:%02d:%05.2f, %c%02d:%02d:%05.2f\n", rh,rm,rs,sign,dd,dm,ds);

	r2hms(ra_jnow, &rh, &rm, &rs);
	d2dms(de_jnow, &sign, &dd, &dm, &ds);
	printf("Target TETE position: %02d:%02d:%05.2f, %c%02d:%02d:%05.2f\n", rh, rm, rs, sign, dd, dm, ds);

	printf("deltaT = %f\n", deltaT);


	//  now compute prepoint positions/times
	//	
	//		de_tete = de_jnow (will remain constant)
	//		ra_tete = GAST(jd_prepoint) - H_target
	//

	printf("Prepoint DateTime, GCRS RA, GCRS DE, TETE RA, TETE DE\n");		// output header line

	de_tete = de_jnow;										// DE in TETE frame does not change
	sidereal_time(jd_target, 0.0, deltaT, 1, 1, 1, &gast);	// gast at Target time
	H_target = gast - ra_jnow / 15;							// target H -> this will remain constant in the TETE frame

	mjd_prepoint = jd_start - MJD_OFFSET;
	mjd_end = jd_start - MJD_OFFSET - (duration_hours / 24.0);
	while (mjd_prepoint > mjd_end)
	{

		// compute new ra_tete from H
		//
		sidereal_time(mjd_prepoint + MJD_OFFSET, 0.0, deltaT, 1, 1, 1, &gast);		// gast at prepoint time (hours)
		ra_tete = (gast - H_target) * 15;								// compute ra_tete (hours -> deg)
		if (ra_tete < 0)
		{
			ra_tete =360.0 + fmod(ra_tete, 360.0);			// module 360.0 degrees
		}
	

		// compute GCRS position from TETE position
		//
		jd_tdb = mjd_prepoint + MJD_OFFSET + (deltaT / 86400.0);
		if (tete2gcrs(jd_tdb, ra_tete, de_tete, &ra_gcrs, &de_gcrs) != 0)
		{
			printf("Error from tete2gcrs.");
			return (EXIT_FAILURE);
		}

		// output prepoint data
		//
		cal_date(mjd_prepoint + MJD_OFFSET, &yr, &mon, &day, &hr);
		ih = (int)hr;
		dblTmp = hr - ih;
		dblTmp = dblTmp * 60.0;
		im = (int)dblTmp;
		dblTmp = dblTmp - im;
		isec = (int)(dblTmp * 60);
		printf("%04d-%02d-%02d %02d:%02d:%02d,", yr, mon, day, ih, im, isec);

		r2hms(ra_gcrs, &rh, &rm, &rs);
		d2dms(de_gcrs, &sign, &dd, &dm, &ds);
		printf(" %02d:%02d:%05.2f, %c%02d:%02d:%05.2f,", rh, rm, rs, sign, dd, dm, ds);
		r2hms(ra_tete, &rh, &rm, &rs);
		d2dms(de_tete, &sign, &dd, &dm, &ds);
		printf(" %02d:%02d:%05.2f, %c%02d:%02d:%05.2f\n", rh, rm, rs, sign, dd, dm, ds);


		// next time
		//
		mjd_prepoint -= interval_sec / 86400.0;

	} // end of loop through prepoint times

	//**********************
	// DONE - no errors
	//**********************
	return(EXIT_SUCCESS);

} // end of main


/*
------------------------------------------------------------------------

	gcrs2tete()

   PURPOSE:
	  This function converts a star's GCRS position to a True Equator/equinox of date position.
	  * Assumes TDB = TT

   REFERENCES:
	  

   INPUT
   ARGUMENTS:
	jd_tdb (double)
		TDB time for TETE (true equator/equinox)
	ra_gcrs (double)
		RA of star in GCRS
	de_gcrs (double)
		DE of star in GCRS

   OUTPUT
   ARGUMENTS:
	ra_tete (double)
		RA of star in TETE
	de_tete
		DE of star in TETE

   RETURNED
   VALUE:
	(int)
		= 0 => no problems	 
		<  0 ... error from function 'vector2radec'

   GLOBALS
   USED:
	  T0                 novascon.c

   FUNCTIONS
   CALLED:
	  frame_tie          novas.c
	  precession         novas.c
	  nutation           novas.c
	  vector2radec       novas.c
	  sin                math.h
	  cos                math.h

   VER./DATE/
   PROGRAMMER:
	  2024-11-13, sbp

   NOTES:


------------------------------------------------------------------------
*/
int gcrs2tete(double jd_tdb, double ra_gcrs, double de_gcrs,
					double* ra_tete, double* de_tete)
{

	double r, d;
	double pos1[3];
	double pos2[3];
	double pos3[3];
	double pos4[3];
	int error;

	/*
	   Form position vector in equatorial system from input coordinates.
	*/

	r = ra_gcrs * DEG2RAD;
	d = de_gcrs * DEG2RAD;

	pos1[0] = cos(d) * cos(r);
	pos1[1] = cos(d) * sin(r);
	pos1[2] = sin(d);

	/*
	   Transform the position vector from GCRS to true equator and equinox
	   of date.
	*/

	frame_tie(pos1, 1, pos2);
	precession(T0, pos2, jd_tdb, pos3);
	nutation(jd_tdb, 0, 0, pos3, pos4);

	/*
	   Convert the position vector to equatorial spherical coordinates.
	*/

	if ((error = vector2radec(pos4, ra_tete, de_tete)) != 0)
	{
		*ra_tete = 0.0;
		*de_tete = 0.0;
		return (-1 * error);
	}
	*ra_tete = *ra_tete * 15;	// convert ra hrs to deg

	return (error);


} // end of gcrs2tete



/*
------------------------------------------------------------------------

	tete2gcrs()

   PURPOSE:
	  This function converts a star's True equator/equinox (of date) position to a GCRS position

	  * Assumes TDB = TT
	  * Accuracy should be better than 10 arc seconds
	  
   REFERENCES:


   INPUT
   ARGUMENTS:
	jd_tdb (double)
		TDB time for equator/equinox
	ra_tete (double)
		RA of star in True Equator/Equniox coordinates
	de_tete
		DE of star in TETE coordinates

   OUTPUT
   ARGUMENTS:
	ra_gcrs (double)
		RA of star in GCRS
	de_gcrs (double)
		DE of star in GCRS

   RETURNED
   VALUE:
	(int)
		= 0 => no problems
		<  0 ... error from function 'vector2radec'

   GLOBALS
   USED:
	  T0                 novascon.c

   FUNCTIONS
   CALLED:
	  frame_tie          novas.c
	  precession         novas.c
	  nutation           novas.c
	  vector2radec       novas.c
	  sin                math.h
	  cos                math.h

   VER./DATE/
   PROGRAMMER:
	  2024-11-13, sbp

   NOTES:


------------------------------------------------------------------------
*/
int tete2gcrs(double jd_tdb, double ra_tete, double de_tete,
				double* ra_gcrs, double* de_gcrs
				)
{

	double r, d;
	double pos1[3];
	double pos2[3];
	double pos3[3];
	double pos4[3];
	int error;

	/*
	   Form position vector in equatorial system from input coordinates.
	*/

	r = ra_tete * DEG2RAD;
	d = de_tete * DEG2RAD;

	pos1[0] = cos(d) * cos(r);
	pos1[1] = cos(d) * sin(r);
	pos1[2] = sin(d);

	/*
	   Transform the position vector from true equator and equinox of date to GCRS
	*/

	nutation(jd_tdb, 1, 0, pos1, pos2);
	precession(jd_tdb, pos2, T0, pos3);
	frame_tie(pos3, 0, pos4);

	/*
	   Convert the position vector to equatorial spherical coordinates.
	*/

	if ((error = vector2radec(pos4, ra_gcrs, de_gcrs)) != 0)
	{
		*ra_gcrs = 0.0;
		*de_gcrs = 0.0;
		return (-1 * error);
	}
	*ra_gcrs = *ra_gcrs * 15;		// convert ra hrs to deg

	return (error);


} // end of tete2gcrs



/*
------------------------------------------------------------------------

	tete2altaz()

   PURPOSE:
	  This function converts a star's True equator/equinox (of date) position to an approximate Alt/Az 
	  * Assumes TDB = TT
	  * Aberration, parallax, and gravitational deflection are ignored.
	  
   REFERENCES:
	* Explanatory Supplement (3rd edition), section 7.1.3

   INPUT
   ARGUMENTS:
	jd_tdb (double)
		TDB time for TETE
	ra_jnow (double)
		RA of star in TETE
	de_jnow
		DE of star in TETE
	obs_lat
		geodetic latitude of observer (north positive, degrees)
	obs_long
		geodetic longitude of observer (east positive, degrees)
	obs_h
		geodetic height of observer (meters)

   OUTPUT
   ARGUMENTS:
	alt (double)
		altitiude of positon
	az (double)
		azimuth of position

   RETURNED
   VALUE:
	(int)
		= 0 => no problems


   GLOBALS
   USED:
	  T0                 novascon.c

   FUNCTIONS
   CALLED:
	  frame_tie          novas.c
	  precession         novas.c
	  nutation           novas.c
	  vector2radec       novas.c
	  sin                math.h
	  cos                math.h

   VER./DATE/
   PROGRAMMER:
	  2024-11-13, sbp

   NOTES:
   Uses TETE as an approximation of topocentric position (ignoring parallax, aberration and other corrections).  But this is
   good enough for pre-pointing and good enough to match JNOW used in telescope mounts.


------------------------------------------------------------------------
*/
void tete2altaz(double jd_ut1, double delta_t, double ra_tete, double de_tete, double obs_lat, double obs_long, double obs_h,
						double* alt, double* az)
{
	double rar, decr;

	on_surface geo_loc = { obs_lat, obs_long, obs_h, 10.0, 1010.0 };

	// compute alt,az
	//
	equ2hor(jd_ut1, delta_t, 1, 0.0, 0.0, &geo_loc, ra_tete, de_tete, 0, alt, az, &rar, &decr);


} // end of tete2altaz

/*
------------------------------------------------------------------------

	d2dms

   PURPOSE:
	  Converts declination degrees to declination deg, min, sec

   REFERENCES:

   INPUT
   ARGUMENTS:
	de (double)

   OUTPUT
   ARGUMENTS:
	dms (text)

   RETURNED
   VALUE:
	(int)
		= 0 => no problems


   GLOBALS
   USED:


   FUNCTIONS
   CALLED:

   VER./DATE/
   PROGRAMMER:
	  2024-11-13, sbp

   NOTES:

------------------------------------------------------------------------
*/
void d2dms(double de, char *sign, int *deg, int *min, double *sec)
{
	double dblTmp;


	*sign = '+';
	if (de < 0)
	{
		*sign = '-';
		de = -de;			// make it positive
	}

	*deg = (int)de;
	dblTmp = de - *deg;
	dblTmp = dblTmp * 60.0;
	*min = (int)dblTmp;
	dblTmp = dblTmp - *min;
	*sec = dblTmp * 60;


} // end of d2dms

/*
------------------------------------------------------------------------

	r2hms

   PURPOSE:
	  Converts RA degrees to hours, min, sec

   REFERENCES:

   INPUT
   ARGUMENTS:
	ra (double)

   OUTPUT
   ARGUMENTS:
	dms (text)

   RETURNED
   VALUE:
	(int)
		= 0 => no problems


   GLOBALS
   USED:


   FUNCTIONS
   CALLED:

   VER./DATE/
   PROGRAMMER:
	  2024-11-13, sbp

   NOTES:

------------------------------------------------------------------------
*/
void r2hms(double ra, int* hrs, int* min, double* sec)
{
	double dblTmp;

	dblTmp = ra / 15;			// convert to hours (0 to 24)
	*hrs = (int)dblTmp;
	dblTmp = dblTmp - *hrs;
	dblTmp = dblTmp * 60.0;
	*min = (int)dblTmp;
	dblTmp = dblTmp - *min;
	*sec = dblTmp * 60;


} // end of r2hms

/*
------------------------------------------------------------------------

	parseRA

   PURPOSE:
	  parse RA input value
	  accepts either deg.dddd or hh:mm:ss.s

   REFERENCES:

   INPUT
   ARGUMENTS:
	arg (char *)

   OUTPUT
   ARGUMENTS:
	ra (deg)

   RETURNED
   VALUE:
	(int)
		= 0 => no problems


   GLOBALS
   USED:


   FUNCTIONS
   CALLED:

   VER./DATE/
   PROGRAMMER:
	  2024-11-13, sbp

   NOTES:

------------------------------------------------------------------------
*/
int parseRA(char* arg, double* ra)
{
	int hh, mm;
	double ss;

	// look for ":" signify sexagesimal HH:MM:SS format
	//
	if (arg[2] == ':')
	{
		// decode hms
		//
		if (sscanf(arg, "%2d", &hh) == 0)
		{
			return (-1);	// error
		}
		arg += 3;
		if (sscanf(arg, "%2d", &mm) == 0)
		{
			return (-1);	// error
		}
		arg += 3;
		if (sscanf(arg, "%4lf", &ss) == 0)
		{
			return(-1);		// error
		}
		
		*ra = (double)((double)hh * 3600 + (double)mm * 60 + ss)*15.0/3600.0;		// ra in degrees
	}
	else
	{
		// decode as degrees ddd.dddd
		//
		if (sscanf(arg, "%lf", ra) == 0)
		{
			return(-1);
		}

	}

	// all done, no issues
	//
	return(0);

} // end of parseRA


/*
------------------------------------------------------------------------

	parseDE

   PURPOSE:
	  parse DE input value.
	  accepts either -/+ deg.ddd  or -/+dd:mm:ss.s

   REFERENCES:

   INPUT
   ARGUMENTS:
	arg (char *)

   OUTPUT
   ARGUMENTS:
	de (deg)

   RETURNED
   VALUE:
	(int)
		= 0 => no problems


   GLOBALS
   USED:


   FUNCTIONS
   CALLED:

   VER./DATE/
   PROGRAMMER:
	  2024-11-13, sbp

   NOTES:

------------------------------------------------------------------------
*/
int parseDE(char* arg, double* de)
{
	double sign;
	int dd, mm;
	double ss;

	// look for ":" signify sexagesimal sDD:MM:SS format
	//
	if (arg[3] == ':')
	{
		// decode hms
		//
		if (arg[0] == '-')
		{
			sign = -1.0;
		}
		else if (arg[0] == '+')
		{
			sign = 1.0;
		}
		else
		{
			return (-1);	// error
		}

		arg++;
		if (sscanf(arg, "%2d", &dd) == 0)
		{
			return (-1);	// error
		}
		arg += 3;
		if (sscanf(arg, "%2d", &mm) == 0)
		{
			return (-1);	// error
		}
		arg += 3;
		if (sscanf(arg, "%4lf", &ss) == 0)
		{
			return(-1);		// error
		}
		*de = sign * (double)((double)dd * 3600 + (double)mm * 60 + ss)  / 3600.0;		// de in degrees

	}
	else
	{
		// decode as degrees ddd.dddd
		//
		if (sscanf(arg, "%lf", de) == 0)
		{
			return(-1);
		}

	}

	// all done, no issues
	//
	return(0);

} // end of parseDE

/*
------------------------------------------------------------------------

	parseDateTime

   PURPOSE:
	  parse date+time and returns julian date
	  format can be either
		jjjjjjj.fffff style julian date
		yyyy-mm-dd_hh:mm:ss

   REFERENCES:

   INPUT
   ARGUMENTS:
	arg (char *)

   OUTPUT
   ARGUMENTS:
	jd = julian date

   RETURNED
   VALUE:
	(int)
		= 0 => no problems


   GLOBALS
   USED:


   FUNCTIONS
   CALLED:

   VER./DATE/
   PROGRAMMER:
	  2024-11-13, sbp

   NOTES:

------------------------------------------------------------------------
*/
int parseDateTime(char* arg, double* jd)
{
	int yy, mon, dd;
	int hh, mm, ss;
	double hr;

	// look for first "-" to signify sexagesimal yyyy-mm-dd_hh:mm:ss format
	//
	if (arg[4] == '-')
	{
		// decode calendar date format
		//
		if (sscanf(arg, "%4d", &yy) == 0)
		{
			return (-1);	// error
		}

		arg += 5;
		if (sscanf(arg, "%2d", &mon) == 0)
		{
			return (-1);	// error
		}
		arg += 3;
		if (sscanf(arg, "%2d", &dd) == 0)
		{
			return (-1);	// error
		}

		arg += 3;
		if (sscanf(arg, "%2d", &hh) == 0)
		{
			return (-1);	// error
		}
		arg += 3;
		if (sscanf(arg, "%2d", &mm) == 0)
		{
			return (-1);	// error
		}
		arg += 3;
		if (sscanf(arg, "%2d", &ss) == 0)
		{
			return(-1);		// error
		}

		hr = (double)((double)hh * 3600 + (double)mm * 60 + ss) / 3600.0;		// hours of day

		*jd = julian_date(yy, mon, dd, hr);			// convert to julian date
	}
	else
	{
		// decode as decimal julian date (e.g. jjjjjjj.jjj)
		//
		if (sscanf(arg, "%lf", jd) == 0)
		{
			return(-1);
		}

	}

	// all done, no issues
	//
	return(0);

} // end of parseDE

