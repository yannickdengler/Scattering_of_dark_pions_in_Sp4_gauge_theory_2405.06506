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
	 /*zeta=generalized_zeta_PL(l,m,dv[0],dv[1],dv[2],m1,m2,q2,lspace,sz);*/

	 if (gen_zeta_error) {
		 printf(" Error [%d] in generalized zeta function\n",gen_zeta_error);
		 return gen_zeta_error;
	 }
/* 
   printf("m1= %f,  m2= %f,  l= %d,  m= %d,\nlspace= %d,  sz=% d,  dv= [%d, %d, %d]\n",m1,m2,l,m,lspace,sz,dv[0],dv[1],dv[2]);
   printf("q2: %f, zeta: %18.14f,%18.14f\n",q2, creal(zeta), cimag(zeta)); */

   printf("%18.14f,%18.14f\n", creal(zeta), cimag(zeta));

   return 0;
 }


/* below is code variant with ECM as input -- this is not a good idea as it does not 
allow for different dispersion relation to be used. Although one should not use them anyway.*/

/*	if (argc != 11){
		printf("Argument count incorrect. Fix eet. [%d]\n",argc);
		printf("proper: m1 m2 E_CM l m space size d[0] d[1] d[2]\n");
		exit(42);
	}
	double m1=atof(argv[1]);
	double m2=atof(argv[2]);
	double E_CM=atof(argv[3]);
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

	double pstar2 = ( pow(E_CM,4) + pow(( pow(m1,2) - pow(m2,2) ),2) - 2.0*pow(E_CM,2)*( pow(m1,2) + pow(m2,2) ) ) / (4.0 * pow(E_CM,2) );
	double q2 = pstar2 * pow(lspace/(2.0*M_PI),2);

	double _Complex zeta;
	zeta=generalized_zeta(l,m,dv,mp1,mp2,q2,sz);

	if (gen_zeta_error) {
		printf(" Error [%d] in generalized zeta function\n",gen_zeta_error);
		return gen_zeta_error;
	}

  printf("m1= %f,  m2= %f,  l= %d,  m= %d,\nlspace= %d,  sz=% d,  dv= [%d, %d, %d]\n",m1,m2,l,m,lspace,sz,dv[0],dv[1],dv[2]);
	printf("E_CM: %f, zeta: %18.14f,%18.14f\n",E_CM, creal(zeta), cimag(zeta));

	return 0;*/