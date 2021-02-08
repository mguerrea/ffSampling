#include "ffsampling.h"

extern complex double *roots[];

void split_fft(t_pol_fft f, t_pol_fft *f0, t_pol_fft *f1)
{
    int pow = 0;
    int n = f.len;
    while (n >>= 1)
        pow++;
    double complex *w = roots[pow - 1];
    f0->len = f.len / 2;
    f1->len = f.len / 2;
    f0->coeffs = malloc(sizeof(complex double)*f0->len);
    f1->coeffs = malloc(sizeof(complex double)*f1->len);
    for (int i = 0; i < f0->len; i++)
    {
        complex double w_conj = creal(w[2*i]) - I*cimag(w[2*i]);
        f0->coeffs[i] = 0.5*(f.coeffs[2 * i] + f.coeffs[2 * i + 1]);
        f1->coeffs[i] = 0.5*(f.coeffs[2 * i] - f.coeffs[2 * i + 1]) * w_conj;
    }
}

// split a polynomial f into two polynomials f0 and f1
void split(t_pol f, t_pol *f0, t_pol *f1)
{
    f0->len = f.len / 2;
    f1->len = f.len / 2;
    f0->coeffs = malloc(sizeof(double)*f0->len);
    f1->coeffs = malloc(sizeof(double)*f1->len);
    for (int i = 0; i < f0->len; i++)
    {
        f0->coeffs[i] = f.coeffs[2 * i];
        f1->coeffs[i] = f.coeffs[2 * i + 1];
    }
}

// merge two polynomials in fft representation into another polynomial
t_pol_fft merge_fft(t_pol_fft f0, t_pol_fft f1)
{
    t_pol_fft f_fft;
    int pow = 0;
    int n = 2 * f0.len;
    while (n >>= 1)
        pow++;
    double complex *w = roots[pow - 1];
    f_fft.len = 2 * f0.len;
    f_fft.coeffs = malloc(sizeof(complex double) * f_fft.len);
    for (int i = 0; i < f_fft.len; i++)
        f_fft.coeffs[i] = 0;
    for (int i = 0; i < f0.len; i++)
    {
        f_fft.coeffs[2 * i] = f0.coeffs[i] + w[2 * i] * f1.coeffs[i];
        f_fft.coeffs[2 * i + 1] = f0.coeffs[i] - w[2 * i] * f1.coeffs[i];
    }
    return (f_fft);
}

t_pol_fft fft(t_pol f)
{
    t_pol f0, f1;
    t_pol_fft f0_fft, f1_fft, f_fft;
    int n = f.len;
    if (n > 2)
    {
        split(f, &f0, &f1);
        f0_fft = fft(f0);
        free(f0.coeffs);
        f1_fft = fft(f1);
        free(f1.coeffs);
        f_fft = merge_fft(f0_fft, f1_fft);
        free(f0_fft.coeffs);
        free(f1_fft.coeffs);
    }
    else if (n == 2)
    {
        f_fft.len = n;
        f_fft.coeffs = malloc(sizeof(complex double) * n);
        f_fft.coeffs[0] = f.coeffs[0] + I * f.coeffs[1];
        f_fft.coeffs[1] = f.coeffs[0] - I * f.coeffs[1];
    }
    return f_fft;
}

t_pol merge(t_pol f0, t_pol f1)
{
    t_pol f;
    f.len = 2 *  f0.len;
    f.coeffs = malloc(sizeof(double)*f.len);
    for (int i = 0; i < f0.len; i++)
    {
        f.coeffs[2 * i] = f0.coeffs[i];
        f.coeffs[2 * i + 1] = f1.coeffs[i];
    }
    return (f);
}

t_pol ifft(t_pol_fft f_fft)
{
    t_pol_fft f0_fft, f1_fft;
    t_pol f, f0, f1; 
    if (f_fft.len > 2)
    {
        split_fft(f_fft, &f0_fft, &f1_fft);
        f0 = ifft(f0_fft);
        free(f0_fft.coeffs);
        f1 = ifft(f1_fft);
        free(f1_fft.coeffs);
        f = merge(f0, f1);
        free(f0.coeffs);
        free(f1.coeffs);
    }
    else if (f_fft.len == 2)
    {
        f.len = 2;
        f.coeffs = malloc(sizeof(double) * f.len);
        f.coeffs[0] = creal(f_fft.coeffs[0]);
        f.coeffs[1] = cimag(f_fft.coeffs[0]);
    }
    return (f);
}