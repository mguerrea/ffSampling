#include "ffsampling.h"
#include "test.h"
#include <assert.h>

char *seed;
extern struct samplerZ_KAT samplerZ[];
extern struct sign_KAT sign[];
extern const t_params params[];
extern struct sign2 sign2[];

static void test_ffsampling(t_pol point, t_sk key, t_params params)
{
    init_RCDT();
    t_pol_fft point_fft = fft(point);
    t_pol_fft a = key.basis[0][0];
    t_pol_fft b = key.basis[0][1];
    t_pol_fft c = key.basis[1][0];
    t_pol_fft d = key.basis[1][1];
    t_pol_fft t[2], v[2];
    t_pol s[2];

    t[0] = new_pol(point_fft.len);
    t[1] = new_pol(point_fft.len);
    v[0] = new_pol(point_fft.len);
    v[1] = new_pol(point_fft.len);

    // we compute t which is a preimage of point (but not a short one)
    for (int i = 0; i < t[0].len; i++)
    {
        t[0].coeffs[i] = (point_fft.coeffs[i] * d.coeffs[i]) / (Q);
        t[1].coeffs[i] = (-point_fft.coeffs[i] * b.coeffs[i]) / (Q);
    }

    // ffSampling gives us a random point z on our lattice close to t
    t_pol_fft *z = ffSampling(t, key.T, params);

    vect_mat_mul(v, z, key.basis);

    s[0] = ifft(v[0]);
    s[1] = ifft(v[1]);

    for (int i = 0; i < v[1].len; i++)
    {
        s[0].coeffs[i] = point.coeffs[i] - s[0].coeffs[i];
        s[1].coeffs[i] = -s[1].coeffs[i];
    }

    printf("s[0] = %f\t%f\n", s[0].coeffs[0], s[0].coeffs[1]);
    printf("s[1] = %f\t%f\n", s[1].coeffs[0], s[1].coeffs[1]);

    t_pol tmp = mul_zq(s[1], key.h);
    t_pol new_point = add_zq(s[0], tmp);

    printf("new = %f\t%f\n", new_point.coeffs[0], new_point.coeffs[1]);
    printf("point = %f\t%f\n", point.coeffs[0], point.coeffs[1]);
}

void test_samplerZ()
{
    int z;
    init_RCDT();
    for (int i = 0; i < 12; i++)
    {
        seed = samplerZ[i].randombytes;
        z = SamplerZ(samplerZ[i].mu, samplerZ[i].sigma, samplerZ[i].sigmin);
        printf("z = %d\n", z);
        assert(samplerZ[i].z == z);
    }
}

int main()
{
    // t_pol f, g, F, G, point;
    // t_sk key;
    // for (int i = 0; i < 1; i++)
    // {
    //     for (int j = 0; j < 1; j++)
    //     {
    //         seed = sign2[j].seed;
    //         f.len = 1 << (i + 1);
    //         f.coeffs = sign2[j].f;
    //         g.len = 1 << (i + 1);
    //         g.coeffs = sign2[j].g;
    //         F.len = 1 << (i + 1);
    //         F.coeffs = sign2[j].F;
    //         G.len = 1 << (i + 1);
    //         G.coeffs = sign2[j].G;
    //         point.len = 1 << (i + 1);
    //         point.coeffs = sign2[j].point;
    //         key = gen_sk(f, g, F, G, params[i].sigma);
    //         printf("h = %f\t%f\n", key.h.coeffs[0], key.h.coeffs[1]);
    //         test_ffsampling(point, key, params[i]);
    //     }
    // }
    test_samplerZ();
    return (0);
}