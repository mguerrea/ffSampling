#include "ffsampling.h"
#include "math.h"

t_pol_fft add_fft(t_pol_fft *res, t_pol_fft f, t_pol_fft g)
{
    for (int i = 0; i < res->len; i++)
        res->coeffs[i] = f.coeffs[i] + g.coeffs[i];
    return (*res);
}

t_pol_fft sub_fft(t_pol_fft *res, t_pol_fft f, t_pol_fft g)
{
    for (int i = 0; i < res->len; i++)
        res->coeffs[i] = f.coeffs[i] - g.coeffs[i];
    return (*res);
}

t_pol_fft mul_fft(t_pol_fft *res, t_pol_fft f, t_pol_fft g)
{
    for (int i = 0; i < res->len; i++)
        res->coeffs[i] = f.coeffs[i] * g.coeffs[i];
    return (*res);
}

t_pol_fft adj_fft(t_pol_fft *res, t_pol_fft f)
{
    for (int i = 0; i < res->len; i++)
        res->coeffs[i] = creal(f.coeffs[i]) - I*cimag(f.coeffs[i]);
    return(*res);
}

t_pol_fft div_fft(t_pol_fft *res, t_pol_fft f, t_pol_fft g)
{
    for (int i = 0; i < res->len; i++)
        res->coeffs[i] = f.coeffs[i] / g.coeffs[i];
    return (*res);
}

t_pol_fft sqrt_fft(t_pol_fft *res, t_pol_fft f)
{
    for (int i = 0; i < res->len; i++)
        res->coeffs[i] = csqrt(f.coeffs[i]);
    return (*res);
}