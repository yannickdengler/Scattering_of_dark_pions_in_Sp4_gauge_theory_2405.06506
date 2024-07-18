#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <spherical_harmonic.h>
#include <gsl/gsl_sf_legendre.h>

void coord(double const *const cart, double *rad)
{
	/*first compute spherical coordinates from cartesian*/
	for(int i=0;i<3;i++){
		rad[i]=0.;
	}

	rad[0]=sqrt( pow(cart[0],2) + pow(cart[1],2) + pow(cart[2],2) );
	if (rad[0] < DBL_EPSILON){
		return;
	}
	rad[1]=acos(cart[2]/rad[0]);
	if(rad[0] < DBL_EPSILON && rad[1] < DBL_EPSILON){
		return;
	}
	rad[2]=atan2(cart[1],cart[0]);
	return;
}

void spherical_poly(int const l_par, int const m_par, double const *const rvec, double _Complex *Y_lm)
{
	/*while this says spherical harmonic it means harmonic polynomial*/
	int l = l_par;
	int m = abs(m_par);

	double r[3];
	coord(rvec,r);

	/*get the legendre polynomial*/
	double const plm=gsl_sf_legendre_sphPlm(l,m,cos(r[1]));
	/*construct spherical polynomial*/
	if (m_par >= 0){
		*Y_lm = pow(r[0],l) * plm * cexp(I * m * r[2]);
	}
	else{
		*Y_lm = pow(-1,m) * pow(r[0],l) * plm * cexp( - I * m * r[2]);
	}
	return;
}

void spherical_harmonic(int const l_par, int const m_par, double const *const rvec, double _Complex *Y_lm)
{
	/*while this says spherical harmonic it means harmonic polynomial*/
	int l = l_par;
	int m = abs(m_par);

	double r[3];
	coord(rvec,r);

	/*get the legendre polynomial*/
	double const plm=gsl_sf_legendre_sphPlm(l,m,cos(r[1]));
	/*construct spherical polynomial*/
	if (m_par >= 0){
		*Y_lm = plm * cexp(I * m * r[2]);
	}
	else{
		*Y_lm = pow(-1,m) * plm * cexp( - I * m * r[2]);
	}
	return;
}

double rl(int const l, double const *const rvec)
{
	/*evalualing r^l -- part of spherical polynomial*/
	double rad[3];
	coord(rvec, rad);
	double rl=pow(rad[0],l);
	return rl;
}
