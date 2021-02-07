#include "ffsampling.h"
#include <time.h>

#define LEN 4

extern t_params params[];

double pol2[][2] = {{72, -37}, {73,-4}, {-80, 11}, {65, 49}};
double pol4[][4] = {{-80, 39, 19, -15}, {12, 16, -9, -44}, {-2, -121, -5, 20}, {-75, -17, -15, -14}};

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
	dim = 1 << n;

	if (dim == 2)
	{
		t_pol f = {2, pol2[0]}, g = {2, pol2[1]}, F = {2, pol2[2]}, G = {2, pol2[3]};
		key = gen_sk(f, g, F, G, params[n - 1].sigma);
	}
	else
	{
		t_pol f = {4, pol4[0]}, g = {4, pol4[1]}, F = {4, pol4[2]}, G = {4, pol4[3]};
		key = gen_sk(f, g, F, G, params[n - 1].sigma);
	}

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