/*===========================================================================*/
/*                        Generalized Zeta Function                          */
/*                                                                           */
/*    This program calculates the generalized Zeta function as defined       */
/*    in Leskovec, Prelovsek, arXiv:hep-lat/1202.2145                        */
/*                                                                           */
/*    It works for scattering of two hadrons (either equal od unequal mass)  */
/*    at zero total momentum and non-zero total momentum.                    */
/*                                                                           */
/*    This is a C implementation trying to use as general as possible code   */
/*    so that it is portable to different systems.                           */
/*                                                                           */
/*    Code is based on Christian B. Lang's fortran implementation            */
/*    swapping some routines for library routines.                           */
/*                                                                           */
/*===========================================================================*/
/*      double _Complex generalized_zeta                                     */
/*             calculates the generalized zeta function                      */
/*             input:                                                        */
/*             int l,m,dv,sz                                                 */
/*             double m1,m2,q2                                               */
/*             output:                                                       */
/*             double _Complex Z_lm^dv                                       */
/*===========================================================================*/

#include <complex.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169
#endif

extern int gen_zeta_error;

struct arg1_params{ 
	int size;
	int l;
	int m;
	int dvec[3];
	double vvec[3];
	double afactor;
	double qsq;
};

struct arg2_params
{
	double qsq;	
};

double _Complex generalized_zeta(int const l, int const m, int const *const dvec, double const mp1, double const mp2, double const q2, int const size);
double _Complex generalized_zeta_PL(int const l, int const m, int const d1, int const d2, int const d3, double const m1, double const m2, double const q2, int const N_L, int const size);
double _Complex first_term(int const size, int const l, int const m, int const *const dvec, double const *const vvec, double const afactor, double const qsq);
double _Complex second_term(int const l, int const m, double const *const vvec, const double qsq);
double _Complex third_term(int const size, int const l, int const m, int const *const dvec, double const *const vvec, double const afactor, double const qsq);
double argument1_void(double t, void *p);
double argument2_void(double t, void *p);

