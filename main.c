#include "ffsampling.h"
#include <time.h>

#define LEN 4

extern t_params params[];

int main()
{
	double tab[][2] = {{72, -37}, {73,-4}, {-80, 11}, {65, 49}};
	t_pol f = {2, tab[0]}, g = {2, tab[1]}, F = {2, tab[2]}, G = {2, tab[3]};
	
	t_sk key = gen_sk(f, g, F, G, params[0].sigma);
	double m[2];
	srand(time(NULL));
	m[0] = rand() % Q;
	m[1] = rand() % Q;
	t_pol message = {2, m}, sig;
	sig = pseudo_sign(message, key, params[0]);
	printf("message = %f\t%f\n", message.coeffs[0], message.coeffs[1]);
	printf("signature = %f\t%f\n", sig.coeffs[0], sig.coeffs[1]);

	if (pseudo_verify(message, sig, key.h, params[0]))
		printf("signature accepted\n");
	else
		printf("signature too large\n");
	return (0);
}