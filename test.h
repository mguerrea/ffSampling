#ifndef TEST_H
#define TEST_H

struct samplerZ_KAT {
    double mu;
    double sigma;
    char *randombytes;
    int z;
};

struct sign_KAT
{
    int n;
    double **f;
    double **g;
    double **F;
    double **G;
    double **s2;
    double **h;
    char **seed;
};

struct sign2 {
    double f[2];
    double g[2];
    double F[2];
    double G[2];
    double h[2];
    char *seed;
};


#endif