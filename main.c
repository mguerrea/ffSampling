#include "ffsampling.h"

#define LEN 8

int main()
{
	double sigmin = 1.277833697;
	int z;
	double mu = -91.90471153063714000;
	double sigma = 1.703799041475491800;

	// init_RCDT();
	// z = SamplerZ(mu, sigma, sigmin);
	// printf("z = %d\n", z);

	// t_pol f;
	// f.len = LEN;
	// f.coeffs = malloc(sizeof(double)*LEN);
	// for (int i = 0; i < LEN; i++)
	// 	f.coeffs[i] = i;

	// t_pol_fft g = fft(f);
	// for (int i = 0; i < g.len; i++)
	// 	printf("%.20f%+.20fi\n", crealf(g.coeffs[i]), cimagf(g.coeffs[i]));
	// free(f.coeffs);
	// free(g.coeffs);
	
	float tab[][2] = {{72, -37}, {73,-4}, {-80, 11}, {65, 49}};
	t_pol f = {2, tab[0]}, g = {2, tab[1]}, F = {2, tab[2]}, G = {2, tab[3]};
	t_sk key = gen_sk(f, g, F, G);
	print_mat(key.gram);
	return (0);
}