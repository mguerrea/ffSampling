#include "ffsampling.h"

#define LEN 4

extern t_params params[];

int main()
{
	double sigmin = 1.277833697;
	int z;
	double mu = -91.90471153063714000;
	double sigma = 1.703799041475491800;

	// init_RCDT();
	// z = SamplerZ(mu, sigma, sigmin);
	// printf("z = %d\n", z);

	float tab[][2] = {{72, -37}, {73,-4}, {-80, 11}, {65, 49}};
	t_pol f = {2, tab[0]}, g = {2, tab[1]}, F = {2, tab[2]}, G = {2, tab[3]};
	sigma = params[0].sigma;
	t_sk key = gen_sk(f, g, F, G, sigma);
	print_tree(key.T);
	printf("\n");
	return (0);
}