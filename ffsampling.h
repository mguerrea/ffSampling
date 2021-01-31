#ifndef FFSAMPLING_H
#define FFSAMPLING_H

#include <gmp.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <complex.h>

#define LN2 0.69314718056
#define SIGMAX 1.8205
#define INV_2SIGMA2 1 / (2 * (SIGMAX * SIGMAX))

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef struct polynomial {
    int len;
    float *coeffs;
}               t_pol;

typedef struct polynomial_fft {
    double complex *coeffs;
    int len;
}                   t_pol_fft;

typedef struct secret_key {
    t_pol_fft basis[2][2];
    t_pol_fft gram[2][2];
    int nb;
}               t_sk;

void init_RCDT();
int SamplerZ(double mu, double sigma, double sigmin);
void random_bytes(int nb, unsigned char *buff);

t_pol_fft fft(t_pol f);
void split_fft(t_pol_fft f, t_pol_fft *f0, t_pol_fft *f1);
t_pol_fft merge_fft(t_pol_fft f0, t_pol_fft f1);
t_sk gen_sk(t_pol f, t_pol g, t_pol F, t_pol G);
void    print_fft(t_pol_fft f);
void    print_mat(t_pol_fft mat[2][2]);

t_pol_fft add(t_pol_fft *res, t_pol_fft f, t_pol_fft g);
t_pol_fft mul(t_pol_fft *res, t_pol_fft f, t_pol_fft g);
t_pol_fft adj(t_pol_fft *res, t_pol_fft f);

#endif