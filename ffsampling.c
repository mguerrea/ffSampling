#include "ffsampling.h"

extern char *seed;

t_pol_fft *ffSampling(t_pol_fft t[2], t_tree *T, t_params params)
{
    t_pol_fft *z1, *z0, *z = malloc(sizeof(t_pol_fft) * 2);
    if (t[0].len == 1)
    {
        double sigma = creal(T->value.coeffs[0]);
        printf("seed = %s\n", seed);
        z[0].coeffs = malloc(sizeof(complex double));
        z[0].len = 1;
        z[0].coeffs[0] = SamplerZ(creal(t[0].coeffs[0]), sigma, params.sigmin);
        z[1].coeffs = malloc(sizeof(complex double));
        z[1].len = 1;
        z[1].coeffs[0] = SamplerZ(creal(t[1].coeffs[0]), sigma, params.sigmin);
        printf("sampler z0 = %f\nsampler z1 = %f\n", creal(z[0].coeffs[0]), creal(z[1].coeffs[0]));
        return(z);
    }
    t_pol_fft l = T->value;
    t_tree *T0 = T->leftchild, *T1 = T->rightchild;
    t_pol_fft t1[2], t0[2];
    split_fft(t[1], &(t1[0]), &(t1[1]));
    z1 = ffSampling(t1, T1, params);
    free(t1[0].coeffs);
    free(t1[1].coeffs);
    z[1] = merge_fft(z1[0], z1[1]);
    t_pol_fft t_prime;
    t_prime.len = t[0].len;
    t_prime.coeffs = malloc(sizeof(complex double) * t_prime.len);
    sub_fft(&t_prime, t[1], z[1]);
    mul_fft(&t_prime, t_prime, l);
    add_fft(&t_prime, t_prime, t[0]);
    split_fft(t_prime, &(t0[0]), &(t0[1]));
    free(t_prime.coeffs);
    z0 = ffSampling(t0, T0, params);
    free(t0[0].coeffs);
    free(t0[1].coeffs);
    z[0] = merge_fft(z0[0], z0[1]);
    return (z);
}