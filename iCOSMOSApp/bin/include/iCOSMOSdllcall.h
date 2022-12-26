#pragma once

/*****************************************************************************************/
/*  iCOSMOS																		Å@       */
/*	"improved Computation of Sensitivities in Model ODE Systems"						 */
/*                                                                                       */
/*	  Atsuko Miyawaki-Kuwakado                                   (2018, July 17)         */
/*	                                                             (2022, May  24)         */
/*****************************************************************************************/

//#include "stdafx.h"//#include <stdio.h>

//#include <stdlib.h>
//#include <malloc.h>
//#include <math.h>
#include <complex>
using namespace std;
//#include <windows.h>
//#include <time.h> 

//--*---*---*---*---*---*---*---*---*--//
//#include <iostream>
//#include <iomanip>
//#include<fstream>
//--*---*---*---*---*---*---*---*---*--//

extern "C" complex<double>*_cdecl C_d_vector(int n);
extern "C" complex<double>**_cdecl C_d_matrix(int m, int n);
typedef void(*po0)(complex<double>* x, complex<double>* fx, complex<double>** v, complex<double>* sv, int nd1, int ni, int cv, complex<double>** stoin, double t);
extern "C" void _cdecl iCOSMOSmain(po0 po);
