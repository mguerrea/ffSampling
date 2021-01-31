#include "ffsampling.h"

static void init_matrix(t_pol_fft mat[2][2], int deg)
{
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            mat[i][j].len = deg;
            mat[i][j].coeffs = malloc(sizeof(complex double) * deg);
        }
}

static void free_matix(t_pol_fft mat[2][2])
{
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            mat[i][j].len = 0;
            free(mat[i][j].coeffs);
            mat[i][j].coeffs = NULL;
        }
}

void    print_fft(t_pol_fft f)
{
    printf("(");
    for (int i = 0; i < f.len; i++)
    {
        printf("%.5f%+.5fi", crealf(f.coeffs[i]), cimagf(f.coeffs[i]));
        if (i < f.len - 1)
            printf(", ");
    }
    printf(")");
}

void    print_mat(t_pol_fft mat[2][2])
{
    printf("[");
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            print_fft(mat[i][j]);
            if (j == 0)
                printf("\t");
        }
        if (i == 0)
            printf("\n");
    }
    printf("]\n");
}

static void gram(t_pol_fft B[2][2], t_pol_fft G[2][2])
{
    t_pol_fft tmp1, tmp2;
    t_pol_fft B_star[2][2];
    tmp1.len = B[0][0].len;
    tmp1.coeffs = malloc(sizeof(complex double) * tmp1.len);
    tmp2.len = B[0][0].len;
    tmp2.coeffs = malloc(sizeof(complex double) * tmp2.len);

    init_matrix(B_star, B[0][0].len);
    for(int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            adj(&(B_star[i][j]), B[j][i]);

    add(&(G[0][0]), mul(&tmp1, B[0][0], B_star[0][0]), mul(&tmp2, B[0][1], B_star[1][0]));
    add(&(G[0][1]), mul(&tmp1, B[0][0], B_star[0][1]), mul(&tmp2, B[0][1], B_star[1][1]));
    add(&(G[1][0]), mul(&tmp1, B[1][0], B_star[0][0]), mul(&tmp2, B[1][1], B_star[1][0]));
    add(&(G[1][1]), mul(&tmp1, B[1][0], B_star[0][1]), mul(&tmp2, B[1][1], B_star[1][1]));
    free(tmp1.coeffs);
    free(tmp2.coeffs);
    free_matix(B_star);
}



t_sk gen_sk(t_pol f, t_pol g, t_pol F, t_pol G)
{
    t_sk key;

    key.basis[0][0] = fft(g);
    key.basis[0][1] = fft(f);
    for (int i = 0; i < key.basis[0][1].len; i++)
        key.basis[0][1].coeffs[i] *= -1;
    key.basis[1][0] = fft(G);
    key.basis[1][1] = fft(F);
    for (int i = 0; i < key.basis[1][1].len; i++)
        key.basis[1][1].coeffs[i] *= -1;
    init_matrix(key.gram, f.len);
    gram(key.basis, key.gram);
    return (key);
}