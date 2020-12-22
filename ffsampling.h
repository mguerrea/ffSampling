#ifndef FFSAMPLING_H
#define FFSAMPLING_H

#include <gmp.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <string.h>

#define LN2 0.69314718056
#define SIGMAX 1.8205

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

void init_RCDT();
int SamplerZ(float mu, float sigma, float sigmin, unsigned char *stream);

#endif