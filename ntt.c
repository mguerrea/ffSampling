#include "ffsampling.h"

extern int *roots_Zq[];

int inv_mod(int a, int b)
{
    int b0 = b, t, q;
	int x0 = 0, x1 = 1;
	if (b == 1)
        return 1;
	while (a > 1) {
		q = a / b;
		t = b;
        b = a % b;
        a = t;
		t = x0;
        x0 = x1 - q * x0;
        x1 = t;
	}
	if (x1 < 0)
        x1 += b0;
	return x1;
}

void split_ntt(t_pol f, t_pol *f0, t_pol *f1)
{
    int pow = 0;
    int n = f.len;
    while (n >>= 1)
        pow++;
    int *w = roots_Zq[pow - 1];
    f0->len = f.len / 2;
    f1->len = f.len / 2;
    f0->coeffs = malloc(sizeof(complex double)*f0->len);
    f1->coeffs = malloc(sizeof(complex double)*f1->len);
    for (int i = 0; i < f0->len; i++)
    {
        f0->coeffs[i] = (int)(I2 * (f.coeffs[2 * i] + f.coeffs[2 * i + 1])) % Q;
        f1->coeffs[i] = (int)(I2 * (f.coeffs[2 * i] - f.coeffs[2 * i + 1]) * inv_mod(w[2*i], Q)) % Q;
    }
}

t_pol merge_ntt(t_pol f0, t_pol f1)
{
    t_pol f;
    int pow = 0;
    int n = 2 * f0.len;
    while (n >>= 1)
        pow++;
    int *w = roots_Zq[pow - 1];
    f.len = 2 * f0.len;
    f.coeffs = malloc(sizeof(complex double) * f.len);
    for (int i = 0; i < f0.len; i++)
    {
        f.coeffs[2 * i] = (int)(f0.coeffs[i] + w[2 * i] * f1.coeffs[i]) % Q;
        f.coeffs[2 * i + 1] = (int)(f0.coeffs[i] - w[2 * i] * f1.coeffs[i]) % Q;
    }
    return (f);
}

t_pol ntt(t_pol f)
{
    t_pol f0, f1, f0_ntt, f1_ntt, f_ntt;
    if(f.len > 2)
    {
        split(f, &f0, &f1);
        f0_ntt = ntt(f0);
        free(f0.coeffs);
        f1_ntt = ntt(f1);
        free(f1.coeffs);
        f_ntt = merge_ntt(f0_ntt, f1_ntt);
        free(f0_ntt.coeffs);
        free(f1_ntt.coeffs);
    }
    else if (f.len == 2)
    {
        f_ntt.len = f.len;
        f_ntt.coeffs = malloc(sizeof(double) * f_ntt.len);
        f_ntt.coeffs[0] = ((int)f.coeffs[0] + SQR1 * (int)f.coeffs[1]) % Q;
        f_ntt.coeffs[0] = (f_ntt.coeffs[0] < 0) ? f_ntt.coeffs[0] + Q : f_ntt.coeffs[0];
        f_ntt.coeffs[1] = ((int)f.coeffs[0] - SQR1 * (int)f.coeffs[1]) % Q;
        f_ntt.coeffs[1] = (f_ntt.coeffs[1] < 0) ? f_ntt.coeffs[1] + Q : f_ntt.coeffs[1];
    }
    return (f_ntt);
}

t_pol intt(t_pol f_ntt)
{
    t_pol f0_ntt, f1_ntt, f, f0, f1; 
    if (f_ntt.len > 2)
    {
        split_ntt(f_ntt, &f0_ntt, &f1_ntt);
        f0 = intt(f0_ntt);
        free(f0_ntt.coeffs);
        f1 = intt(f1_ntt);
        free(f1_ntt.coeffs);
        f = merge(f0, f1);
        free(f0.coeffs);
        free(f1.coeffs);
    }
    else if (f_ntt.len == 2)
    {
        f.len = 2;
        f.coeffs = malloc(sizeof(double) * f.len);
        f.coeffs[0] = (I2 * (int)(f_ntt.coeffs[0] + f_ntt.coeffs[1])) % Q;
        f.coeffs[0] = (f.coeffs[0] < 0) ? f.coeffs[0] + Q : f.coeffs[0];
        f.coeffs[1] = ((I2 * ((int)f_ntt.coeffs[0] - (int)f_ntt.coeffs[1]) % Q) * inv_mod(SQR1, Q)) % Q;
        f.coeffs[1] = (f.coeffs[1] < 0) ? f.coeffs[1] + Q : f.coeffs[1];
    }
    return (f);
}