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
#define Q 12289
#define SQR1 1479
#define I2 6145

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef struct polynomial {
	int len;
	double *coeffs;
}               t_pol;

typedef struct polynomial_fft {
	double complex *coeffs;
	int len;
}                   t_pol_fft;

typedef struct polynomial_ntt {
	int *coeffs;
	int len;
}                   t_pol_ntt;

typedef struct tree {
	t_pol_fft value;
	struct tree *leftchild;
	struct tree *rightchild;
}               t_tree;

typedef struct secret_key {
	t_pol_fft basis[2][2];
	t_tree *T;
	t_pol h;
}               t_sk;

typedef struct params {
	int n;
	double sigma;
	double sigmin;
	int bound;
}               t_params;

void init_RCDT();
int SamplerZ(double mu, double sigma, double sigmin);
void random_bytes(int nb, unsigned char *buff);

t_pol_fft   fft(t_pol f);
t_pol       ifft(t_pol_fft f_fft);
void        split_fft(t_pol_fft f, t_pol_fft *f0, t_pol_fft *f1);
t_pol_fft   merge_fft(t_pol_fft f0, t_pol_fft f1);
t_pol       intt(t_pol_ntt f_ntt);
t_pol_ntt	ntt(t_pol f);

t_sk gen_sk(t_pol f, t_pol g, t_pol F, t_pol G, double sigma);

void	init_matrix(t_pol_fft mat[2][2], int deg);
void	free_matrix(t_pol_fft mat[2][2]);
void    print_fft(t_pol_fft f);
void    print_mat(t_pol_fft mat[2][2]);
void	print_tree(t_tree *T);


t_tree		*new_node(t_pol_fft value);
t_pol_fft	dup_pol(t_pol_fft f);
t_pol_fft	new_pol(int len);

void    vect_mat_mul(t_pol_fft res[2], t_pol_fft vect[2], t_pol_fft mat[2][2]);
int inv_mod(int a, int b);
void split(t_pol f, t_pol *f0, t_pol *f1);
t_pol merge(t_pol f0, t_pol f1);


t_pol_fft *ffSampling(t_pol_fft t[2], t_tree *T, t_params params);
t_pol pseudo_sign(t_pol message, t_sk key, t_params params);
int pseudo_verify(t_pol message, t_pol sig, t_pol h, t_params params);

/*
** pol_op.c
*/

t_pol_fft	add_fft(t_pol_fft *res, t_pol_fft f, t_pol_fft g);
t_pol_fft	mul_fft(t_pol_fft *res, t_pol_fft f, t_pol_fft g);
t_pol_fft	adj_fft(t_pol_fft *res, t_pol_fft f);
t_pol_fft	div_fft(t_pol_fft *res, t_pol_fft f, t_pol_fft g);
t_pol_fft	sub_fft(t_pol_fft *res, t_pol_fft f, t_pol_fft g);
t_pol_fft	sqrt_fft(t_pol_fft *res, t_pol_fft f);
t_pol		div_zq(t_pol f, t_pol g);
t_pol		mul_zq(t_pol f, t_pol g);
t_pol		add_zq(t_pol f, t_pol g);
t_pol		sub_zq(t_pol f, t_pol g);

// int round(float nb);

#endif