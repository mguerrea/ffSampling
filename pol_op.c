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

t_pol_ntt div_ntt(t_pol_ntt *res, t_pol_ntt f, t_pol_ntt g)
{
    for (int i = 0; i < f.len; i++)
    {
        res->coeffs[i] = (int)(f.coeffs[i] * inv_mod(g.coeffs[i], Q)) % Q;
        res->coeffs[i] = (res->coeffs[i] < 0) ? res->coeffs[i] + Q : res->coeffs[i];
    }
    return (*res);
}

t_pol_ntt mul_ntt(t_pol_ntt *res, t_pol_ntt f, t_pol_ntt g)
{
    for (int i = 0; i < f.len; i++)
    {
        res->coeffs[i] = (int)(f.coeffs[i] * g.coeffs[i]) % Q;
        res->coeffs[i] = (res->coeffs[i] < 0) ? res->coeffs[i] + Q : res->coeffs[i];
    }
    return (*res);
}

t_pol_ntt add_ntt(t_pol_ntt *res, t_pol_ntt f, t_pol_ntt g)
{
    for (int i = 0; i < f.len; i++)
    {
        res->coeffs[i] = (int)(f.coeffs[i] + g.coeffs[i]) % Q;
        res->coeffs[i] = (res->coeffs[i] < 0) ? res->coeffs[i] + Q : res->coeffs[i];
    }
    return (*res);
}

t_pol_ntt sub_ntt(t_pol_ntt *res, t_pol_ntt f, t_pol_ntt g)
{
    for (int i = 0; i < f.len; i++)
    {
        res->coeffs[i] = (int)(f.coeffs[i] - g.coeffs[i]) % Q;
        res->coeffs[i] = (res->coeffs[i] < 0) ? res->coeffs[i] + Q : res->coeffs[i];
    }
    return (*res);
}

t_pol div_zq(t_pol f, t_pol g)
{
    t_pol res;
    t_pol_ntt f_ntt, g_ntt;
    f_ntt = ntt(f);
    g_ntt = ntt(g);
    div_ntt(&f_ntt, f_ntt, g_ntt);
    res = intt(f_ntt);
    free(f_ntt.coeffs);
    free(g_ntt.coeffs);
    return (res);
}

t_pol mul_zq(t_pol f, t_pol g)
{
    t_pol res;
    t_pol_ntt f_ntt, g_ntt;
    f_ntt = ntt(f);
    g_ntt = ntt(g);
    mul_ntt(&f_ntt, f_ntt, g_ntt);
    res = intt(f_ntt);
    free(f_ntt.coeffs);
    free(g_ntt.coeffs);
    return (res);
}

t_pol add_zq(t_pol f, t_pol g)
{
    t_pol res;
    t_pol_ntt f_ntt, g_ntt;
    f_ntt = ntt(f);
    g_ntt = ntt(g);
    add_ntt(&f_ntt, f_ntt, g_ntt);
    res = intt(f_ntt);
    free(f_ntt.coeffs);
    free(g_ntt.coeffs);
    return (res);
}

t_pol sub_zq(t_pol f, t_pol g)
{
    t_pol res;
    t_pol_ntt f_ntt, g_ntt;
    f_ntt = ntt(f);
    g_ntt = ntt(g);
    sub_ntt(&f_ntt, f_ntt, g_ntt);
    res = intt(f_ntt);
    free(f_ntt.coeffs);
    free(g_ntt.coeffs);
    return (res);
}