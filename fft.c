#include "ffsampling.h"

extern complex double *roots[];

// split a polynomial f into two polynomials f0 and f1
void split(t_pol f, t_pol *f0, t_pol *f1)
{
    f0->len = f.len / 2;
    f1->len = f.len / 2;
    f0->coeffs = malloc(sizeof(float)*f0->len);
    f1->coeffs = malloc(sizeof(float)*f1->len);
    for (int i = 0; i < f0->len; i++)
    {
        f0->coeffs[i] = f.coeffs[2 * i];
        f1->coeffs[i] = f.coeffs[2 * i + 1];
    }
}

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
    f_fft.len = 2 * f0.len;
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
        free(f.coeffs);
        f0_fft = fft(f0);
        f1_fft = fft(f1);
        f_fft = merge_fft(f0_fft, f1_fft);
        free(f0_fft.coeffs);
        free(f1_fft.coeffs);
    }
    else if (n == 2)
    {
        f_fft.len = n;
        f_fft.coeffs = malloc(sizeof(complex double) * n);
        for (int i = 0; i < n; i++)
            f_fft.coeffs[i] = 0;
        f_fft.coeffs[0] = f.coeffs[0] + I * f.coeffs[1];
        f_fft.coeffs[1] = f.coeffs[0] - I * f.coeffs[1];

    }
    return f_fft;
}