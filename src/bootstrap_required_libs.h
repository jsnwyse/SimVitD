/*Functions for the fitting of Latent Class Analysis models using MCMC
	methods. Two implementations of the model are included in the Bayesian
	formulation: collapsed and not collapsed.
	
	Author:	Jason Wyse,
			School of Computer Science and Statistics,
			Lloyd Institute,
			Trinity College,
			Dublin 2,
			Ireland.
			mailto: wyseja@tcd.ie
			
	Last modification of this code: Mon 09 May 2016 15:38:03 IST   */

#ifndef __SIMVITD_REQUIRED_LIBS__
#define __SIMVITD_REQUIRED_LIBS__

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#include <R_ext/Visibility.h>

//definitions

#define TRUE 1
#define FALSE 0
#define log_2_pi 1.83787706640934533908
#define M_EULER   0.57721566490153286060651209008 /* Euler constant */

#endif
