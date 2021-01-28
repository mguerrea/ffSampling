#include "ffsampling.h"

#define LEN 16

int main()
{
	double sigmin = 1.277833697;
	int z;
	double mu = -91.90471153063714000;
	double sigma = 1.703799041475491800;

	// init_RCDT();
	// z = SamplerZ(mu, sigma, sigmin);
	// printf("z = %d\n", z);

	// t_pol *f = NULL;

	// for (int i = 0; i < 5; i++)
	// 	addCoeff(&f, newCoeff(i));
	// printPol(f);

	t_pol f;
	f.len = LEN;
	f.coeffs = malloc(sizeof(double)*LEN);
	for (int i = 0; i < LEN; i++)
		f.coeffs[i] = i;

	t_pol_fft g = fft(f);
	for (int i = 0; i < g.len; i++)
		printf("%.20f%+.20fi\n", crealf(g.coeffs[i]), cimagf(g.coeffs[i]));
	return (0);
}