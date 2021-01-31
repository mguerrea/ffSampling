#include "ffsampling.h"

t_pol_fft add(t_pol_fft *res, t_pol_fft f, t_pol_fft g)
{
    for (int i = 0; i < res->len; i++)
        res->coeffs[i] = f.coeffs[i] + g.coeffs[i];
    return (*res);
}

t_pol_fft mul(t_pol_fft *res, t_pol_fft f, t_pol_fft g)
{
    for (int i = 0; i < res->len; i++)
        res->coeffs[i] = f.coeffs[i] * g.coeffs[i];
    return (*res);
}

t_pol_fft adj(t_pol_fft *res, t_pol_fft f)
{
    for (int i = 0; i < res->len; i++)
        res->coeffs[i] = creal(f.coeffs[i]) - I*cimag(f.coeffs[i]);
    return(*res);
}