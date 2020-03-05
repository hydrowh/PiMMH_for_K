#pragma once
#define MT_N 624
#define MT_M 397
#define MATRIX_A (unsigned long)2567483615  /* constant vector a */
#define UPPER_MASK (unsigned long)2147483648 /* most significant w-r bits */
#define LOWER_MASK (unsigned long)2147483647 /* least significant r bits */

static unsigned long mt[MT_N]; /* the array for the state vector  */
static int mti = MT_N + 1; /* mti==MT_N+1 means mt[MT_N] is not initialized */

						   /* initializes mt[MT_N] with a seed */
void init_genrand(unsigned long s)
{
	mt[0] = s & (unsigned long)4294967295;
	for (mti = 1; mti<MT_N; mti++) {
		mt[mti] =
			(1812433253 * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		mt[mti] &= (unsigned long)4294967295;
	}
}

void init_by_array(unsigned long init_key[], int key_length)
{
	int i, j, k;
	init_genrand(19650218);
	i = 1; j = 0;
	k = (MT_N>key_length ? MT_N : key_length);
	for (; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525))
			+ init_key[j] + j; /* non linear */
		mt[i] &= (unsigned long)4294967295; /* for WORDSIZE > 32 machines */
		i++; j++;
		if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
		if (j >= key_length) j = 0;
	}
	for (k = MT_N - 1; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941)) - i; /* non linear */
		mt[i] &= (unsigned long)4294967295; /* for WORDSIZE > 32 machines */
		i++;
		if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
	}

	mt[0] = (unsigned long)2147483648; /* MSB is 1; assuring non-zero initial array */
}

unsigned long genrand_int32(void)
{
	unsigned long y;
	static unsigned long mag01[2] = { (unsigned long)0, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= MT_N) { /* generate N words at one time */
		int kk;

		if (mti == MT_N + 1)   /* if init_genrand() has not been called, */
			init_genrand(5489); /* a default initial seed is used */

		for (kk = 0; kk<MT_N - MT_M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 1];
		}
		for (; kk<MT_N - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 1];
		}
		y = (mt[MT_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 1];

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & (unsigned long)2636928640;
	y ^= (y << 15) & (unsigned long)4022730752;
	y ^= (y >> 18);

	return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
	return (long)(genrand_int32() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real01(void)
{
	return genrand_int32()*(1.0 / 4294967295.0);
}

/* generates a random number on [0,1)-real-interval */
double genrand_real01_(void)
{
	return genrand_int32()*(1.0 / 4294967296.0);
}

/* generates a random number on (0,1)-real-interval */
double genrand_real_01_(void)
{
	return (((double)genrand_int32()) + 0.5)*(1.0 / 4294967296.0);
}

/* generates a random number on [0,1) with 53-bit resolution*/
double Uniform(void)
{
	unsigned long a = genrand_int32() >> 5, b = genrand_int32() >> 6;
	return (a*67108864.0 + b)*(1.0 / 9007199254740992.0);
}

double rand_exp(double lambda) { //パラメータλの指数分布に従う乱数
	return -log(Uniform()) / lambda;
}

double rand_normal(double mu, double sigma) {//平均mu, 分散sigmaの正規分布に従う乱数
	double z = sqrt(-2.0*log(Uniform())) * sin(2.0*M_PI*Uniform());
	return mu + sqrt(sigma)*z;
}

double rand_gamma(double kappa, double lamda) {//形状母数 kappa、尺度母数 theta ガンマ分布に従う乱数

	int int_kappa;
	double frac_kappa;
	double theta; 
	theta = 1.0 / lamda;
	int_kappa = (int)kappa;
	frac_kappa = kappa - (double)int_kappa;

	double u, uu;
	double b, p, x_frac, x_int;
	int i;

	/*integer part*/
	x_int = 0;
	for (i = 0; i<int_kappa; i++) {
		x_int += -log(Uniform()); // add expnential random number with mean 1
	}

	/*fractional part*/
	if (fabs(frac_kappa) < 0.01) x_frac = 0;

	else {
		b = (exp(1.0) + frac_kappa) / exp(1.0);
		while (1) {

			u = Uniform();
			p = b*u;

			uu = Uniform();

			if (p <= 1.0) {
				x_frac = pow(p, 1.0 / frac_kappa);
				if (uu <= exp(-x_frac)) break;
			}

			else {
				x_frac = -log((b - p) / frac_kappa);
				if (uu <= pow(x_frac, frac_kappa - 1.0)) break;
			}

		}
	}

	return (x_int + x_frac)*theta;
}



double rand_beta(double A, double B) {
	double tmp = rand_gamma(A, 1.0);
	return tmp / (tmp + rand_gamma(B, 1.0));
}


double prob_normal(double x, double m, double s) {
	return exp(-1.0 * (x - m) * (x - m) / (2.0 * s)) / sqrt(2.0 * M_PI * s);
}

double prob_gamma(double x, double kappa, double lamda) {
	return  pow(x, kappa - 1) * exp(-lamda * x) * pow(lamda, kappa) / tgamma(kappa);
}

double prob_beta(double x, double A, double B) {
	double tmp = tgamma(A)*tgamma(B) / tgamma(A + B);
	return pow(x, A - 1)*pow(1 - x, B - 1) / tmp;
}

double rand_tnormal(double left, double right, double m, double s) {
	double GaussNodeX[21];
	double GaussNodeW[20];
	double integral_l;
	double integral_r;
	double mean = m;
	double sigma = sqrt(s);
	double left_v = left;
	double right_v = right;

	const double c_0 = 2.515517;
	const double c_1 = 0.802853;
	const double c_2 = 0.010328;
	const double d_0 = 1.432788;
	const double d_1 = 0.189269;
	const double d_2 = 0.001308;
	double t;
	double ans;

	double alpha = (left - mean) / sigma;
	double beta = (right - mean) / sigma;

	double a, b, c, center;
	double x, tmp;
	double u = Uniform();

	GaussNodeX[0] = -9.931285e-1;
	GaussNodeX[1] = -9.639719e-1;
	GaussNodeX[2] = -9.122344e-1;
	GaussNodeX[3] = -8.391169e-1;
	GaussNodeX[4] = -7.463319e-1;
	GaussNodeX[5] = -6.360536e-1;
	GaussNodeX[6] = -5.108670e-1;
	GaussNodeX[7] = -3.737060e-1;
	GaussNodeX[8] = -2.277858e-1;
	GaussNodeX[9] = -7.652652e-2;
	GaussNodeX[10] = 7.652652e-2;
	GaussNodeX[11] = 2.277858e-1;
	GaussNodeX[12] = 3.737060e-1;
	GaussNodeX[13] = 5.108670e-1;
	GaussNodeX[14] = 6.360536e-1;
	GaussNodeX[15] = 7.463319e-1;
	GaussNodeX[16] = 8.391169e-1;
	GaussNodeX[17] = 9.122344e-1;
	GaussNodeX[18] = 9.639719e-1;
	GaussNodeX[19] = 9.931285e-1;

	GaussNodeW[0] = 1.761400e-2;
	GaussNodeW[1] = 4.060142e-2;
	GaussNodeW[2] = 6.267204e-2;
	GaussNodeW[3] = 8.327674e-2;
	GaussNodeW[4] = 1.019301e-1;
	GaussNodeW[5] = 1.181945e-1;
	GaussNodeW[6] = 1.316886e-1;
	GaussNodeW[7] = 1.420961e-1;
	GaussNodeW[8] = 1.491729e-1;
	GaussNodeW[9] = 1.527533e-1;
	GaussNodeW[10] = 1.527533e-1;
	GaussNodeW[11] = 1.491729e-1;
	GaussNodeW[12] = 1.420961e-1;
	GaussNodeW[13] = 1.316886e-1;
	GaussNodeW[14] = 1.181945e-1;
	GaussNodeW[15] = 1.019301e-1;
	GaussNodeW[16] = 8.327674e-2;
	GaussNodeW[17] = 6.267204e-2;
	GaussNodeW[18] = 4.060142e-2;
	GaussNodeW[19] = 1.761400e-2;

	integral_l = 0.0;
	if (alpha < 0.0) {
		a = alpha;
		b = 0.0;
		c = -1.0;
	}
	else {
		a = 0.0;
		b = alpha;
		c = 1.0;
	}
	center = (b - a) * 0.5;
	for (int i = 0; i < 20; i++) {
		x = center * GaussNodeX[i] + (a + b) * 0.5;
		tmp = exp(-0.5 * (x * x)) / sqrt(2.0 * M_PI);
		integral_l += GaussNodeW[i] * tmp;
	}
	if (alpha == 0.0) {
		integral_l = 0.5;
	}
	integral_l = 0.5 + c * center * integral_l;


	integral_r = 0.0;
	if (beta < 0.0) {
		a = beta;
		b = 0.0;
		c = -1.0;
	}
	else {
		a = 0.0;
		b = beta;
		c = 1.0;
	}
	center = (b - a) * 0.5;
	for (int i = 0; i < 20; i++) {
		x = center * GaussNodeX[i] + (a + b) * 0.5;
		tmp = exp(-0.5 * (x * x)) / sqrt(2.0 * M_PI);
		integral_r += GaussNodeW[i] * tmp;
	}
	if (beta == 0.0) {
		integral_r = 0.5;
	}
	integral_r = 0.5 + c * center * integral_r;


	tmp = integral_l + u * (integral_r - integral_l);

	if (tmp < 0.0)tmp = 0.0;
	if (tmp > 1.0)tmp = 1.0;
	if (tmp < 0.5) {
		t = sqrt(-2.0 * log(tmp));
		t -= ((c_2 * t + c_1) * t + c_0) / (((d_2 * t + d_1) * t + d_0) * t + 1.0);
		t *= -1.0;
	}
	else {
		tmp = 1.0 - tmp;
		t = sqrt(-2.0 * log(tmp));
		t -= ((c_2 * t + c_1) * t + c_0) / (((d_2 * t + d_1) * t + d_0) * t + 1.0);
	}
	ans = t * sigma + mean;
	if (ans < left_v)ans = left_v;
	if (ans > right_v)ans = right_v;

	return ans;
}

