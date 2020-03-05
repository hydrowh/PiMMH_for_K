#pragma once
#define _USE_MATH_DEFINES
#ifdef _MSC_VER
#include <math.h>
#endif

#ifdef __unix
#include <cmath>
#endif
#include <gmp.h>
#include <mpfr.h>
#include <float.h>
#include "prob_dis.h"
#include "Init.h"

void Q_theta_d(vector<double> &theta) {
	double r0, r1, r2;
	r0 = rand_tnormal(-1.0, 1.0, 0.9, 0.1);
	r1 = rand_gamma(1.0, 0.01);
	r2 = rand_gamma(1.0, 1.0);

	theta[0] = r0;
	theta[1] = (1.0 / r1);
	theta[2] = (1.0 / r2);
}

double P_theta_p(vector<double> &theta) {
	double p;
	p = prob_normal(theta[0], 0.9, 0.1)*prob_gamma(1.0 / theta[1], 1.0, 0.01)*prob_gamma(1.0 / theta[2], 1.0, 1.0);
	if (isnan(p)) {
		p = DBL_MIN;
	}
	return p;
}

double Q_theta_p(vector<double> &theta) {
        double p;
	p = prob_normal(theta[0],0.9,0.1)*prob_gamma(1.0/theta[1], 1.0, 0.01)*prob_gamma(1.0/theta[2], 1.0, 1.0);
        if (isnan(p)) {
		p = DBL_MIN;
	}
	return p;
}

void Xtk_init(vector<double> &X, vector<double> &Y, vector<double> &theta) {
	X[0] = rand_normal(0.0, 1.0);
}

double Wtk_init(vector<double> &X, vector<double> &Y, vector<double> &theta) {
	double p1, p2, p3;
	double tmp1, tmp2;
	tmp1 = Y[0] * exp(-X[0]);
	tmp2 = theta[2];
	p1 = exp(-X[0]) * prob_normal(tmp1, 0.0, tmp2);

	tmp1 = X[0];
	p2 = prob_normal(tmp1, 0.0, 1.0);
	
	tmp1 = X[0];
	p3 = prob_normal(tmp1, 0.0, 1.0);

	return p1 * p2 / p3;
}


void rand_Xt1__yt1_xakt_theta(vector<double> &X_t1, vector<double> &Y, vector<double> &X_akt, vector<double> &theta)
{
	double tmp1, tmp2;
	tmp1 = theta[0] * X_akt[0];
	tmp2 = theta[1];
	X_t1[0] = tmp1 + rand_normal(0.0, tmp2);
}

double P_yt1__xt1_theta(vector<double>&Y, vector<double> &X_t1, vector<double> &theta) {
	double p;
	double tmp1, tmp2;
	tmp1 = Y[0] * exp(-X_t1[0]);
	tmp2 = theta[2];
	p = exp(-X_t1[0]) * prob_normal(tmp1, 0.0, tmp2);
	if (isnan(p)) {
		p = DBL_MIN;
	}
	return p;
}

double P_xt1__xakt_theta(vector<double> &X_t1, vector<double> &X_akt, vector<double> &theta) {
	double p;
	double tmp1, tmp2;
	tmp1 = X_t1[0] - theta[0] * X_akt[0];
	tmp2 = theta[1];
	p = prob_normal(tmp1, 0.0, tmp2);
	if (isnan(p)) {
		p = DBL_MIN;
	}
	return p;
}

double Q_xt1__yt1_xakt_theta(vector<double> &X_t1, vector<double>&Y, vector<double> &X_akt, vector<double> &theta){
	double p;
	double tmp1, tmp2;
	tmp1 = X_t1[0] - theta[0] * X_akt[0];
	tmp2 = theta[1];
	p = prob_normal(tmp1, 0.0, tmp2);
	if(isnan(p)){
		p = DBL_MIN;
	}
	return p;
}

int D_gmp(mpf_t *Wij, mpf_t RP, int Ix) {
	int ans;
	int i;
	mpf_t sum, sum_t, rand;
	
	mpf_init(sum);
	mpf_init(sum_t);
	mpf_init(rand);
	
	mpf_set_d(sum, 0.0);
	mpf_set_d(sum_t, 0.0);
	for (i = 0; i < Ix ; i++) {
		mpf_set(sum_t, sum);
		mpf_add(sum, sum_t, Wij[i]);
	}
	mpf_mul(rand, RP, sum);
	mpf_set_d(sum, 0.0);
	mpf_set_d(sum_t, 0.0);
	ans = 0;
	for (i = 0; i < Ix ; i++) {
		ans = i;
		mpf_set(sum_t, sum);
		mpf_add(sum, sum_t, Wij[i]);
		if(0 <= mpf_cmp(rand, sum_t) && 0 > mpf_cmp(rand, sum)){
			break;
		}
	}

	mpf_clear(sum);
	mpf_clear(sum_t);
	mpf_clear(rand);
	
	return ans;
}

int Draw_akt(vector<double> &W) {
	double sum;
	double tmp, unir;
	int ans;
	sum = 0.0;
	tmp = 0.0;
	int i;
	for (i = 0; i < W.size() ; i++) {
		sum += W[i];
	}
	unir = Uniform() * sum;
	ans = 0;
	sum = 0.0;
	for (i = 0; i < W.size(); i++) {
		ans = i;
		tmp = sum;
		sum += W[i];
		if (tmp <= unir && unir < sum) {
			break;
		}
	}
	return ans;
}

int W_resample(vector<vector<double> > &Wij, int j_i) {
	double sum;
	double tmp, unir;
	int ans;
	sum = 0.0;
	int i;
	for (i = 0; i < I; i++) {
		sum += Wij[i][j_i];
	}
	unir = Uniform() * sum;
	ans = 0;
	sum = 0.0;
	for (i = 0; i < I; i++) {
		ans = i;
		tmp = sum;
		sum += Wij[i][j_i];
		if (tmp <= unir && unir < sum) {
			break;
		}
	}
	return ans;
}


