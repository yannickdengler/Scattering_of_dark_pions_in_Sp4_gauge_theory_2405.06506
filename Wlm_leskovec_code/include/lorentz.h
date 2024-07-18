/*===========================================================================*/
/*    double lorentz_gamma_inv                                               */
/*          inverse lorentz transform                                        */
/*          input double: a[3], v[3]                                         */
/*          g = gammainv(v) acting on a                                      */
/*          g = gammainv*a_parallel + a_transversal                          */
/*          gamma=1/sqrt(1-v.v)                                              */
/*          output double: g[3]                                              */
/*===========================================================================*/
/*    double lorentz_gamma                                                   */
/*          lorentz transform                                                */
/*          input double: a[3], v[3]                                         */
/*          g = gamma(v) acting on a                                         */
/*          g = gamma*a_parallel + a_transversal                             */
/*          gamma=1/sqrt(1-v.v)                                              */
/*          output double: g[3]                                              */
/*===========================================================================*/

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169
#endif

void lorentz_gamma_inv(double const *const a, double const *const v, double *g);
void lorentz_gamma(double const *const a, double const *const v, double *g);
double lorentz_gamma_scalar(double const E, int const *const d, int const Nx);