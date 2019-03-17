
#ifndef _CMB_DEFS_H_
#define _CMB_DEFS_H_

/*useful macros*/
#define DMIN(a,b) GSL_MIN_DBL(a,b)
#define DMAX(a,b) GSL_MAX_DBL(a,b)
#define IMAX(a,b) GSL_MAX_INT(a,b)
#define IMIN(a,b) GSL_MIN_INT(a,b)
#define ISQR(x) ((x)*(x))

/******************************************/


/******************************************/
#define TRUE 1
#define FALSE 0
#define CMB_success 1
#define log_2_pi 1.837877066409345
#define log_det_Q 1. /*log generalized determinant of
											precision matrix for IGMRF on lattice*/
#define SYNC_LOWER -3.3 /*lower limit for the free parameter on synchotron*/
#define SYNC_UPPER -2.7 /*upper limit for the free param on synchotron*/
#define DUST_LOWER 1. /*lower limit for the free parameter on galactic dust*/ 
#define DUST_UPPER 2. /*upper limit for the free parameter on galactic dust*/ 

/*experimental*/
//#define free_lower -3.
//#define free_upper -2.3

#define SHAPE_PRECISION 1.
#define RATE_PRECISION 0.0001

/*for testing- placing priors on the real line*/
#define mu_params 0.
#define precision_params 1.


//const double PLANCK = 6.626068*pow(10.,-34.); /*Planck's constant*/
#define PLANCK 6.626068*pow(10.,-34.)
//const double BOLTZMANN = 1.3806503*pow(10.,-23.); /*Boltzmann's constant*/
#define BOLTZMANN 1.3806503*pow(10.,-23.)

//(log10(6.626) - 34) - (log10(1.38) - 23) - log10(T_0)

#define log10PLANCK -33.1787441114855639057 //- 34.
#define log10BOLTZMANN -22.859916308396282858 //- 23.

//#define PB (6.626068/1.3806503)*pow(10.,-11.)
//const double T0 = 2.725; /*average CMB temp*/
#define T0 2.725
#define log10T0 0.435366506612661297027
//const double T1 = 18.1; /*Galactic dust temperature*/
#define T1 18.1
#define log10T1 1.25767857486918455123
#define C0  (PLANCK/(BOLTZMANN*T0))
#define C1  (PLANCK/(BOLTZMANN*T1))
#define log10c 8.4768207029279274423

#define FITS_NULL_VAL 3800000.00

#define FINITE_DIFFERENCE_STEPSIZE .01
#define SQRT_EPS 1.490116*pow(10,-8.)

#define BLAS_LEVEL2 2
#define BLAS_LEVEL3 3
#define BLAS_level 3


#define WMAP_DATA 0
#define WMAP_DATA_9yr 1
#define PLANCK_DATA 2
#define CUSTOM 3
#define TEST_DATA 4

static const double free_free_factors[5] = {4.661,2.130, 1.361, .605, .273 };

#endif
