/*                        Generalized Zeta Function                          */

#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>
#include <lorentz.h>
#include <gen_zeta.h>
#include <spherical_harmonic.h>

#ifndef M_PI 
#define M_PI 3.141592653589793238462643383279502884197169
#endif

 int main(int argc, char **argv){
	 /*getting input from line*/
	 if (argc != 11){
 		 printf("Argument count incorrect. Fix eet. [%d]\n",argc);
		 printf("proper: m1 m2 q2 l m space size d[0] d[1] d[2]\n");
		 exit(42);
   }

     /*parse input*/
	 double m1=atof(argv[1]);
	 double m2=atof(argv[2]);
	 double q2=atof(argv[3]);
	 int l=atoi(argv[4]);
	 int m=atoi(argv[5]);
	 int lspace=atoi(argv[6]);
	 int sz=atoi(argv[7]);
	 int dv[3];
	 dv[0]=atoi(argv[8]);
	 dv[1]=atoi(argv[9]);
	 dv[2]=atoi(argv[10]);

	 double mp1= m1 * lspace / (2.0 * M_PI);
	 double mp2= m2 * lspace / (2.0 * M_PI);

	 /*turning off GSL error handle*/
	 gsl_set_error_handler_off();

	 double _Complex zeta;
	 zeta=generalized_zeta(l,m,dv,mp1,mp2,q2,sz);

	 if (gen_zeta_error) {
		 printf(" Error [%d] in generalized zeta function\n",gen_zeta_error);
		 return gen_zeta_error;
	 }

	 /*get w_lm*/
 	 double twopiL = 2.0*M_PI/((double) lspace);	 

	 double dbl_l = (double) l;
	 double d2 = pow(dv[0],2)+pow(dv[1],2)+pow(dv[2],2);

	 double Est = (2.0*M_PI/lspace)*( sqrt(pow(mp1,2) + q2) + sqrt(pow(mp2,2) + q2) );
	 double E = sqrt(pow(Est,2) + pow(twopiL,2)*d2);
	 double gam_val = lorentz_gamma_scalar(E, dv,lspace);
	 double complex q = csqrt(q2);
	 double complex denom = pow(M_PI,3./2.)*sqrt(2*dbl_l+1.0)*gam_val*cpow(q,dbl_l+1.0);

     printf("pi32: %f\n", pow(M_PI,3./2.));
     printf("sqrt2lp1: %f\n", sqrt(2*dbl_l+1.0));
     printf("E_Est: %f\n",E/Est );
     printf("gamm: %f\n",gam_val );
     printf("qlp1: %f %f\n", creal(cpow(q,dbl_l+1.0)), cimag(cpow(q,dbl_l+1.0)) );

	 double complex wlm = zeta/denom;				/* eq 96 Rummukainen Gottlieb*/

	 printf("m1= %f,  m2= %f,  l= %d,  m= %d,\nlspace= %d,  sz=% d,  dv= [%d, %d, %d]\n",m1,m2,l,m,lspace,sz,dv[0],dv[1],dv[2]);
     printf("E, %f, q2: %f, Zlm: %18.14f,%18.14f\n",E, q2, creal(zeta), cimag(zeta));
     printf("E, %f, q2: %f, wlm: %18.14f,%18.14f\n",E, q2, creal(wlm), cimag(wlm));
   
   return 0;
 }