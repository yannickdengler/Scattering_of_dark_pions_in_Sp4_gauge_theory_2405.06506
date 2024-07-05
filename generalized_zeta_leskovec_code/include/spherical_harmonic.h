/*===========================================================================*/
/*    void spherical harmonic                                                */
/*          calculates the spherical polynomial r^l Y_lm                     */
/*          using the GSL library for Y_lm                                   */
/*          input double: rvec[3]                                            */
/*          input int: l, m                                                  */
/*          output double _Complex: Y_lm                                     */
/*===========================================================================*/

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169
#endif

void spherical_poly(int const l_par, int const m_par, double const *const rvec, double _Complex *Y_lm);
void coord(double const *const cart, double *rad);
void spherical_harmonic(int const l_par, int const m_par, double const *const rvec, double _Complex *Y_lm);
double rl(int const l, double const *const rvec);
