#include "ffsampling.h"
#include "test.h"
#include <assert.h>

char *seed;
extern struct samplerZ_KAT samplerZ[];
extern struct sign_KAT sign[];
extern const t_params params[];
extern struct sign2 sign2[];

void test_ffsampling(t_pol h, t_sk key, t_params params)
{
    init_RCDT();
    t_pol_fft point = fft(h);
    t_pol_fft a = key.basis[0][0];
    t_pol_fft b = key.basis[0][1];
    t_pol_fft c = key.basis[1][0];
    t_pol_fft d = key.basis[1][1];
    t_pol_fft t0, t1, t[2];
    t0.len = point.len;
    t1.len = point.len;
    t0.coeffs = malloc(sizeof(complex double) * t0.len);
    t1.coeffs = malloc(sizeof(complex double) * t1.len);
    for (int i = 0; i < t0.len; i++)
    {
        t0.coeffs[i] = (point.coeffs[i] * d.coeffs[i]) / (Q);
        t1.coeffs[i] = (-point.coeffs[i] * b.coeffs[i]) / (Q);
    }
    t[0] = t0;
    t[1] = t1;
    t_pol_fft *z = ffSampling(t, key.T, params);
    print_fft(z[0]);
    printf("\n");
    print_fft(z[1]);
}

void test_samplerZ()
{
    double sigmin = 1.277833697;
	int z;
	init_RCDT();
    for (int i = 0; i < 8; i++)
    {
        seed = samplerZ[i].randombytes;
	    z = SamplerZ(samplerZ[i].mu, samplerZ[i].sigma, sigmin);
	    printf("z = %d\n", z);
        assert(samplerZ[i].z == z);
    }
}

int main()
{
    t_pol f, g, F, G, h;
    t_sk key;
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            seed = sign2[j].seed;
            f.len = 1 << (i + 1);
            f.coeffs = sign2[j].f;
            g.len = 1 << (i + 1);
            g.coeffs = sign2[j].g;
            F.len = 1 << (i + 1);
            F.coeffs = sign2[j].F;
            G.len = 1 << (i + 1);
            G.coeffs = sign2[j].G;
            h.len = 1 << (i + 1);
            h.coeffs = sign2[j].h;
            key = gen_sk(f, g, F, G, params[i].sigma);
            test_ffsampling(h, key, params[i]);
        }
    }    
	return (0);
}