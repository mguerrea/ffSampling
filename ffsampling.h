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
    float *coeffs;
    int len;
}               t_pol;

typedef struct polynomial_fft {
    double complex *coeffs;
    int len;
}                   t_pol_fft;

void init_RCDT();
int SamplerZ(double mu, double sigma, double sigmin);
void random_bytes(int nb, unsigned char *buff);

t_pol_fft fft(t_pol f);

#endif