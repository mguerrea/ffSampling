#include "ffsampling.h"

static void LDL(t_pol_fft G[2][2], t_pol_fft L[2][2], t_pol_fft D[2][2])
{
    t_pol_fft tmp1, tmp2;
    tmp1.len = G[0][0].len;
    tmp1.coeffs = malloc(sizeof(complex double) * tmp1.len);
    tmp2.len = G[0][0].len;
    tmp2.coeffs = malloc(sizeof(complex double) * tmp2.len);
    
    for (int i = 0; i < G[0][0].len; i++)
        D[0][0].coeffs[i] = G[0][0].coeffs[i];
    div_fft(&(L[1][0]), G[1][0], G[0][0]);
    mul_fft(&tmp2, adj_fft(&tmp1, L[1][0]), L[1][0]);
    mul_fft(&tmp1, tmp2, G[0][0]);
    sub_fft(&(D[1][1]), G[1][1], tmp1);
    free(tmp1.coeffs);
    free(tmp2.coeffs);
}

t_tree *ffLDL(t_pol_fft G[2][2])
{
    t_tree *T;
    t_pol_fft L[2][2], D[2][2], G0[2][2], G1[2][2];
    init_matrix(L, G[0][0].len);
    init_matrix(D, G[0][0].len);
    LDL(G, L, D);
    free_matrix(G);
    T = new_node(L[1][0]);
    if (D[0][0].len == 2)
    {
        T->leftchild = new_node(D[0][0]);
        T->rightchild = new_node(D[1][1]);
    }
    else
    {
        split_fft(D[0][0], &(G0[0][0]), &(G0[0][1]));
        split_fft(D[1][1], &(G1[0][0]), &(G1[0][1]));
        adj_fft(&(G0[1][0]), G0[0][1]);
        adj_fft(&(G1[1][0]), G1[0][1]);
        G0[1][1] = dup_pol(G0[0][0]);
        G1[1][1] = dup_pol(G1[0][0]);
        T->leftchild = ffLDL(G0);
        T->rightchild= ffLDL(G1);
    }
    free_matrix(L);
    free_matrix(D);
    return (T);
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
            adj_fft(&(B_star[i][j]), B[j][i]);

    add_fft(&(G[0][0]), mul_fft(&tmp1, B[0][0], B_star[0][0]), mul_fft(&tmp2, B[0][1], B_star[1][0]));
    add_fft(&(G[0][1]), mul_fft(&tmp1, B[0][0], B_star[0][1]), mul_fft(&tmp2, B[0][1], B_star[1][1]));
    add_fft(&(G[1][0]), mul_fft(&tmp1, B[1][0], B_star[0][0]), mul_fft(&tmp2, B[1][1], B_star[1][0]));
    add_fft(&(G[1][1]), mul_fft(&tmp1, B[1][0], B_star[0][1]), mul_fft(&tmp2, B[1][1], B_star[1][1]));
    free(tmp1.coeffs);
    free(tmp2.coeffs);
    free_matrix(B_star);
}

void normalize_tree(t_tree *T, double sigma)
{
    if (T->leftchild)
    {
        normalize_tree(T->leftchild, sigma);
        normalize_tree(T->rightchild, sigma);
    }
    else
    {
        sqrt_fft(&(T->value), T->value);
        T->value.coeffs[0] = sigma / T->value.coeffs[0];
        T->value.coeffs[1] = 0;
        // for (int i = 0; i < T->value.len; i++)
        //     T->value.coeffs[i] = sigma / T->value.coeffs[i];
    }
}

t_sk gen_sk(t_pol f, t_pol g, t_pol F, t_pol G, double sigma)
{
    t_sk key;
    t_pol_fft gram_basis[2][2];

    key.basis[0][0] = fft(g);
    key.basis[0][1] = fft(f);
    for (int i = 0; i < key.basis[0][1].len; i++)
        key.basis[0][1].coeffs[i] *= -1;
    key.basis[1][0] = fft(G);
    key.basis[1][1] = fft(F);
    for (int i = 0; i < key.basis[1][1].len; i++)
        key.basis[1][1].coeffs[i] *= -1;
    init_matrix(gram_basis, f.len);
    gram(key.basis, gram_basis);
    key.T = ffLDL(gram_basis);
    normalize_tree(key.T, sigma);
    key.h = div_zq(g, f);
    return (key);
}