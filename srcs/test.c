#include "ffsampling.h"
#include "test.h"

char *seed;
extern struct samplerZ_KAT samplerZ[];

void test_samplerZ()
{
    int z;
    init_RCDT();
    printf(KBLU "Testing SamplerZ\n\n");
    printf("%10s\t|%10s\t|%10s\t|%8s\n", "mu", "sigma", "sigmin", "exp z");
    printf("---------------------------------");
    printf("--------------------------\n");
    for (int i = 0; i < 12; i++)
    {
        printf(KBLU "%12f\t|%12f\t|%12f\t|%5d", samplerZ[i].mu, samplerZ[i].sigma, samplerZ[i].sigmin, samplerZ[i].z);
        seed = samplerZ[i].randombytes;
        z = SamplerZ(samplerZ[i].mu, samplerZ[i].sigma, samplerZ[i].sigmin);
        assert(samplerZ[i].z == z);
        printf(KGRN "\tOK\n");
    }
    free_RCDT();
}

void test_fft()
{
    printf(KBLU "Testing FFT\n\n");
    srand(time(NULL));
    for (int i = 1; i < 8; i++)
    {
        int dim = 1 << i;
        printf(KBLU "Degree = %-5d", dim);
        t_pol f = {.len = dim}, g = {.len = dim};
        f.coeffs = malloc(sizeof(double) * f.len);
        g.coeffs = malloc(sizeof(double) * g.len);
        for (int j = 0; j < dim; j++)
        {
            f.coeffs[j] = rand() % 20 - 10;
            g.coeffs[j] = rand() % 20 - 10;
        }
        t_pol_fft f_fft = fft(f);
        t_pol_fft g_fft = fft(g);
        mul_fft(&g_fft, f_fft, g_fft);
        div_fft(&g_fft, g_fft, f_fft);
        t_pol h = ifft(g_fft);
        for (int j = 0; j < h.len; j++)
            assert(round(h.coeffs[j]) == round(g.coeffs[j]));
        free(f.coeffs);
        free(g.coeffs);
        free(f_fft.coeffs);
        free(g_fft.coeffs);
        free(h.coeffs);
        printf(KGRN "OK\n" KNRM);
    }
}

void test_ntt()
{
    printf(KBLU "Testing NTT\n\n");
    srand(time(NULL));
    for (int i = 1; i < 8; i++)
    {
        int dim = 1 << i;
        printf(KBLU "Degree = %-5d", dim);
        t_pol f = {.len = dim}, g = {.len = dim};
        f.coeffs = malloc(sizeof(double) * f.len);
        g.coeffs = malloc(sizeof(double) * g.len);
        for (int j = 0; j < dim; j++)
        {
            f.coeffs[j] = rand() % Q;
            g.coeffs[j] = rand() % Q;
        }
        t_pol h = mul_zq(f, g);
        t_pol k = div_zq(h, f);
        for (int j = 0; j < h.len; j++)
            assert(round(k.coeffs[j]) == round(g.coeffs[j]));
        free(f.coeffs);
        free(g.coeffs);
        free(h.coeffs);
        free(k.coeffs);
        printf(KGRN "OK\n" KNRM);
    }
}

int main()
{
    printf("\n");
    test_samplerZ();
    printf(KMAG "\n********************\n\n");
    test_fft();
    printf(KMAG "\n********************\n\n");
    test_ntt();
    return (0);
}