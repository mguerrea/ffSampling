#include "ffsampling.h"
#include <time.h>

#define LEN 4

extern t_params params[];
extern double pol2[][2], pol4[][4], pol8[][8], pol16[][16], pol32[][32],
	pol64[][64], pol128[][128];

int usage(char **argv)
{
	printf("Usage: %s [dim] [-s|-v] [message]\n", argv[0]);
	return (1);
}

int main(int argc, char **argv)
{
	int dim = 2, n = 0;
	t_sk key;

	if (argc < 4)
		return (usage(argv));
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
	else if (dim == 64)
		f.coeffs = pol64[0], g.coeffs = pol64[1], F.coeffs = pol64[2], G.coeffs = pol64[3];
	else if (dim == 128)
		f.coeffs = pol128[0], g.coeffs = pol128[1], F.coeffs = pol128[2], G.coeffs = pol128[3];
	else
	{
		printf("Dimension not supported\n");
		printf("Supported dimensions are: 2, 4, 8, 16, 32\n");
		return (1);
	}
	key = gen_sk(f, g, F, G, params[n - 1].sigma);

	if (strcmp(argv[2], "-s") == 0)
	{
		t_pol sig = pseudo_sign(argv[3], key, params[n - 1]);
		FILE *out = fopen("message.sig", "w");
		if (out == NULL)
			return(printf("Could not create message.sig\n"));
		for (int i = 0; i < sig.len; i++)
			fprintf(out, "%f\n", sig.coeffs[i]);
		printf("message: %s\n", argv[3]);
		printf("signature written to message.sig\n");
		fclose(out);
		free(sig.coeffs);
	}
	else if (strcmp(argv[2], "-v") == 0)
	{
		FILE *in = fopen("message.sig", "r");
		t_pol sig = {.len = dim};
		if (in == NULL)
			return(printf("Could not open message.sig\n"));
		sig.coeffs = malloc(sizeof(double) * sig.len);
		for (int i = 0; i < dim; i++)
			fscanf(in, "%lf", &(sig.coeffs[i]));
		if (pseudo_verify(argv[3], sig, key.h, params[n - 1]))
			printf("signature accepted\n");
		else
			printf("signature too large\n");
		fclose(in);
		free(sig.coeffs);
	}

	free_sk(key);
	return (0);
}