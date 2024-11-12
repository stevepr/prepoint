// prepoint.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string.h>
#include <math.h>


#define USAGE   "Usage: prepoint --ra=rrrr.rrr --de=dddd.ddd --app --jd==jjjj.jjj [--hours=h.hhhh | --days=d.dddd] --interval=sec \n"
#define EXIT_FAILURE	1
#define cLF "\n"
#define cCRLF "\r\n"
#define	C_IMAXLINE		1000					/* Max chars in an input line */

// externals
//
extern "C" int checkout_stars(void);

int main(int argc, char* argv[])
{
	char* arg;
	double dblTmp;

	// Inputs
	//
	double ra_target = -100.0;			// target RA in deg (ICRS)
	double de_target = -100.0;			// target DE in deg (ICRS)
	double jd_target = 0.0;				// target time for prepoint location
	double duration_hours = 24.0;		// duration of prepoint times list (hours)
	double interval_sec = 10.0;			// interval between prepoint times (seconds)

	// outputs
	//
	double ra_app;						// Apparent place of target RA deg
	double de_app;						// Apparent place of target DE deg


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

		if (strncmp(arg + 2, "ra=", 3) == 0)
		{
			arg += 5;		// first char of numeric value
			if (sscanf(arg, "%lf", &ra_target) == 0)
			{
				printf("error parsing ra= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}
		else if (strncmp(arg + 2, "de=", 3) == 0)
		{
			arg += 5;		// first char of numeric value
			if (sscanf(arg, "%lf", &de_target) == 0)
			{
				printf("error parsing de= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}
		else if (strncmp(arg + 2, "jd=", 3) == 0)
		{
			arg += 5;		// first char of numeric value
			if (sscanf(arg, "%lf", &jd_target) == 0)
			{
				printf("error parsing jd= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}
		else if (strncmp(arg + 2, "hours=", 6) == 0)
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
		else if (strncmp(arg + 2, "days=", 5) == 0)
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
		else if (strncmp(arg + 2, "interval=", 9) == 0)
		{
			arg += 11;		// first char of numeric value
			if (sscanf(arg, "%lf", &interval_sec) == 0)
			{
				printf("error parsing interval= value.\n");
				printf(USAGE);
				return(EXIT_FAILURE);
			}
		}


	} // end of command line parsing

	// validate input arguments
	//

	if (interval_sec < 0)
	{
		printf("argument error: interval < 0.\n");
		printf(USAGE);
		return(EXIT_FAILURE);
	}


    // run basic star validation
    //
    checkout_stars();


}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
