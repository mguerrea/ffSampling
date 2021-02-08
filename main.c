#include "ffsampling.h"
#include <time.h>

#define LEN 4

extern t_params params[];
extern double pol2[][2], pol4[][4], pol8[][8], pol16[][16], pol32[][32];

t_pol gen_message(int n)
{
	t_pol m;
	m.len = n;
	m.coeffs = malloc(sizeof(double) * m.len);
	srand(time(NULL));
	for (int i = 0; i < n; i++)
		m.coeffs[i] = rand() % Q;
	return (m);
}

int main(int argc, char **argv)
{
	int dim = 2, n = 0;
	t_sk key;

	if (argc > 1)
		dim = atoi(argv[1]);
	while (dim >>= 1)
		n++;
	dim = params[n - 1].n;
	
	t_pol f = {.len = dim}, g = {.len = dim}, F = {.len = dim}, G = {.len = dim};
	if (dim == 2)
		f.coeffs = pol2[0], g.coeffs = pol2[1], F.coeffs = pol2[2], G.coeffs = pol2[3];
	else if (dim == 4)
		f.coeffs = pol4[0], g.coeffs = pol4[1], F.coeffs = pol4[2], G.coeffs = pol4[3];
	else if (dim == 8)
		f.coeffs = pol8[0], g.coeffs = pol8[1], F.coeffs = pol8[2], G.coeffs = pol8[3];
	else if (dim == 16)
		f.coeffs = pol16[0], g.coeffs = pol16[1], F.coeffs = pol16[2], G.coeffs = pol16[3];
	else if (dim == 32)
		f.coeffs = pol32[0], g.coeffs = pol32[1], F.coeffs = pol32[2], G.coeffs = pol32[3];
	else
	{
		printf("Dimension not supported\n");
		printf("Supported dimensions are: 2, 4, 8, 16, 32\n");
		return (1);
	}
	key = gen_sk(f, g, F, G, params[n - 1].sigma);

	t_pol message = gen_message(dim);
	printf("message = ");
	for (int i = 0; i < message.len; i++)
		printf("%f\t", message.coeffs[i]);

	t_pol sig = pseudo_sign(message, key, params[n - 1]);
	printf("\nsignature = ");
	for (int i = 0; i < sig.len; i++)
		printf("%f\t", sig.coeffs[i]);

	if (pseudo_verify(message, sig, key.h, params[n - 1]))
		printf("\nsignature accepted\n");
	else
		printf("\nsignature too large\n");

	return (0);
}