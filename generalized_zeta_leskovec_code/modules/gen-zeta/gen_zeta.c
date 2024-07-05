#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <float.h>
#include <lorentz.h>
#include <gen_zeta.h>
#include <spherical_harmonic.h>

int gen_zeta_error = 0;

double Z_argument1_void(double t, void *p)
{
	struct arg1_params *params = (struct arg1_params *)p;

	int size = (params->size);
	int l = (params->l);
	int m = (params->m);
	double afactor = (params->afactor);
	double qsq = (params->qsq);
	int *dvec = (params->dvec);
	double *vvec = (params->vvec);

	double _Complex a=0.0;
	double _Complex add;
	double _Complex rYlm;
	double uvec[3];
	double rtvec[3];
	double gvec[3];
	double rvec[3];
	double r2;
	double arg1;

	for(int uz=-size; uz<=size; uz++){
		for(int uy=-size;uy<=size; uy++){
			for(int ux=-size;ux<=size; ux++){
				uvec[0]=ux;
				uvec[1]=uy;
				uvec[2]=uz;
				lorentz_gamma(uvec,vvec,gvec);

				for(int j=0;j<3;j++){
					rvec[j] = M_PI * gvec[j];
					rtvec[j] = - M_PI * gvec[j]/t;
				}
				r2=pow(rvec[0],2) + pow(rvec[1],2) + pow(rvec[2],2);
				if (r2>0.0){
					add=cexp(I * M_PI * afactor * (ux*dvec[0] + uy*dvec[1] + uz*dvec[2]));
					spherical_poly(l, m, rtvec, &rYlm);
					a += add * rYlm * exp(-r2/t);
				}

			}
		}
	}
	/* cabs(a) gives wrong sign some places. needs to be creal*/
	arg1 = exp(t * qsq) * creal(a) / sqrt(pow(t,3));
	return arg1;
}

double _Complex Z_first_term(int const size, int const l, int const m, int const *const dvec, double const *const vvec, double const afactor, double const qsq)
{
	const double fac1=1.0/sqrt(1.0 - pow(vvec[0],2) - pow(vvec[1],2) - pow(vvec[2],2));
	const double gfac=fac1 * sqrt(pow(M_PI,3)) * cpow(-I,l);

	/*gsl expects this kind of passing of the functions it integrates*/
	struct arg1_params params;
	params.size=size;
	params.l=l;
	params.m=m;
	for(int i=0; i<3; i++){
		params.dvec[i]=dvec[i];
		params.vvec[i]=vvec[i];
	}
	params.afactor=afactor;
	params.qsq=qsq;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(4000);
  	double result, error;
	gsl_function F;
  	F.function = &Z_argument1_void;
  	F.params = &params;
	gsl_integration_qags (&F, 0.0, 1.0, 1e-11, 1e-13, 4000, w, &result, &error);
    gsl_integration_workspace_free (w);
	double first;
	first= gfac * result;

	return first;
}

double Z_argument2_void(double t, void *p)
{
	struct arg2_params *params = (struct arg2_params *)p;
	double qsq = (params->qsq);
	double arg2;
	arg2 = (exp(t * qsq) - 1.0)/sqrt(pow(t,3));
	return arg2;
}

double _Complex Z_second_term(int const l, int const m, double const *const vvec, const double qsq)
{
	if (l != 0 || m!=0){
		return 0+I*0;
	}

	const double gfac=1.0/sqrt(1.0 - pow(vvec[0],2) - pow(vvec[1],2) - pow(vvec[2],2));

	struct arg2_params params;
	params.qsq=qsq;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(4000);

  	double result, error;
	gsl_function F;
  	F.function = &Z_argument2_void;
  	F.params = &params;
	gsl_integration_qags (&F, 0.0, 1.0, 1e-11, 1e-13, 4000, w, &result, &error);
    gsl_integration_workspace_free (w);
	double second;
	second = gfac * sqrt(pow(M_PI,3)) / sqrt(4.0 * M_PI) * result - gfac * M_PI;
	return second;
}

double _Complex Z_third_term(int const size, int const l, int const m, int const *const dvec, double const *const vvec, double const afactor, double const qsq)
{
	double _Complex third=0.0 + I*0.0;
	double _Complex hp=0.0 + I * 0.0;
	double nvec[3];
	double rvec[3];
	double r2;

	/*note: it used to be ni+0.5*afactor*dvec[i] which was a bug*/
	/*this is now fixed*/
	for(int nz = -size; nz <= size; nz++){
		/*nvec[2]=nz + 0.50 * afactor * dvec[2];*/
		nvec[2]=nz - 0.50 * afactor * dvec[2];
		for(int ny = -size; ny <= size; ny++){
			nvec[1]=ny - 0.50 * afactor * dvec[1];
			for(int nx = -size; nx <= size; nx++){
				nvec[0]=nx - 0.50 * afactor * dvec[0];
				lorentz_gamma_inv(nvec,vvec,rvec);
				r2 = pow(rvec[0],2) + pow(rvec[1],2) + pow(rvec[2],2);
				if (fabs(r2 - qsq) > DBL_EPSILON ){ /*same as if (r2 != qsq) {*/
					spherical_poly(l, m, rvec,&hp);
					third += hp * exp( - (r2 - qsq) ) / (r2 - qsq);
				}
				else {
/* 					printf("issue in third term: infinities\n");
 */					gen_zeta_error = 1;
				}

			}
		}
	}

	return third;
}

double _Complex generalized_zeta(int const l, int const m, int const *const dvec, double const mp1, double const mp2, double const q2, int const size)
{
	/*E_CM_2: center of mass energy squared*/
	double E_CM_2=pow( sqrt(pow(mp1,2)+q2) + sqrt(pow(mp2,2)+q2) ,2);

	/*normalized d vector*/
	double vvec[3];
	for(int i=0; i<3; i++){
		vvec[i] = dvec[i] / sqrt(pow(dvec[0],2) + pow(dvec[1],2) + pow(dvec[2],2) + E_CM_2);
	}

	/*the A factor*/
	double afactor= 1.0 + (pow(mp1,2) - pow(mp2,2)) / E_CM_2;

	double qsq=q2;
	double _Complex t1,t2,t3;
	t1=Z_first_term(size, l, m, dvec, vvec, afactor, qsq);
	t2=Z_second_term(l, m, vvec, qsq);
	t3=Z_third_term(size, l, m, dvec, vvec, afactor, qsq);
	double _Complex gen_zeta=t1+t2+t3;

	return gen_zeta;
}

/* double _Complex generalized_zeta(int const l, int const m, int const *const dvec, double const mp1, double const mp2, double const q2, int const size)*/
double _Complex generalized_zeta_PL(int const l, int const m, int const d1, int const d2, int const d3, double const m1, double const m2, double const q2, int const N_L, int const size)
{

	double mp1= m1 * N_L / (2.0 * M_PI);
	double mp2= m2 * N_L / (2.0 * M_PI);
	/*E_CM_2: center of mass energy squared*/
	double E_CM_2=pow( sqrt(pow(mp1,2)+q2) + sqrt(pow(mp2,2)+q2) ,2);

	/*normalized d vector*/
	double vvec[3];
	vvec[0] = d1 / sqrt(pow(d1,2) + pow(d2,2) + pow(d3,2) + E_CM_2);
	vvec[1] = d2 / sqrt(pow(d1,2) + pow(d2,2) + pow(d3,2) + E_CM_2);
	vvec[2] = d3 / sqrt(pow(d1,2) + pow(d2,2) + pow(d3,2) + E_CM_2);
	/* double vvec[3];
	 for(int i=0; i<3; i++){
	 	vvec[i] = dvec[i] / sqrt(pow(dvec[0],2) + pow(dvec[1],2) + pow(dvec[2],2) + E_CM_2);
	 }*/

	int dvec[3];
	dvec[0] = d1;
	dvec[1] = d2;
	dvec[2] = d3;
	
	/*the A factor*/
	double afactor= 1.0 + (pow(mp1,2) - pow(mp2,2)) / E_CM_2;

	double qsq=q2;
	double _Complex t1,t2,t3;
	t1=Z_first_term(size, l, m, dvec, vvec, afactor, qsq);
	t2=Z_second_term(l, m, vvec, qsq);
	t3=Z_third_term(size, l, m, dvec, vvec, afactor, qsq);
	double _Complex gen_zeta=t1+t2+t3;

	return gen_zeta;
}
