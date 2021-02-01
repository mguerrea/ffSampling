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
};


#endif