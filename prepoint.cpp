// prepoint.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string.h>
#include <math.h>
#include "novas.h"

#define USAGE   "Usage: prepoint --ra=rrrr --de=dddd --app --target_time=tttt [--hours=h.hhhh | --days=d.dddd] --interval=sec \n"
#define EXIT_FAILURE	1
#define cLF "\n"
#define cCRLF "\r\n"
#define	C_IMAXLINE		1000					/* Max chars in an input line */

// externals
//
extern "C" int checkout_stars(void);

int main(int argc, char* argv[])
{
    std::cout << "Hello World!\n";

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
