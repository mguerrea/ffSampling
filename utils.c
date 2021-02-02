#include "ffsampling.h"

void print_tree(t_tree *T)
{
        print_fft(T->value);
        if (T->leftchild)
        {
        printf("\n[");
        print_tree(T->leftchild);
        print_tree(T->rightchild);
        printf("]");
        }
}

t_tree *new_node(t_pol_fft value)
{
    t_tree *new;
    new = malloc(sizeof(t_tree));
    new->value = dup_pol(value);
    new->leftchild = NULL;
    new->rightchild = NULL;
    return (new);
}

t_pol_fft dup_pol(t_pol_fft f)
{
    t_pol_fft new;
    new.len = f.len;
    new.coeffs = malloc(sizeof(complex double) * new.len);
    for (int i = 0; i < new.len; i++)
        new.coeffs[i] = f.coeffs[i];
    return (new);
}

void init_matrix(t_pol_fft mat[2][2], int deg)
{
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            mat[i][j].len = deg;
            mat[i][j].coeffs = malloc(sizeof(complex double) * deg);
        }
}

void free_matrix(t_pol_fft mat[2][2])
{
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            mat[i][j].len = 0;
            free(mat[i][j].coeffs);
            mat[i][j].coeffs = NULL;
        }
}

t_pol_fft new_pol(int len)
{
    t_pol_fft new;
    new.len = len;
    new.coeffs = malloc(sizeof(complex double) * len);
    for (int i = 0; i < len; i++)
        new.coeffs[i] = 0;
    return (new);
}

void    print_fft(t_pol_fft f)
{
    printf("(");
    for (int i = 0; i < f.len; i++)
    {
        if (creal(f.coeffs[i]))
            printf("%f", creal(f.coeffs[i]));
        if (cimag(f.coeffs[i]))
            printf("%+fi", cimag(f.coeffs[i]));
        if (f.coeffs[i] == 0)
            printf("0");
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

void    vect_mat_mul(t_pol_fft res[2], t_pol_fft vect[2], t_pol_fft mat[2][2])
{
    t_pol_fft new[2];
    t_pol_fft tmp1 = new_pol(vect[0].len);
    t_pol_fft tmp2 = new_pol(vect[0].len);
    new[0] = new_pol(vect[0].len);
    new[1] = new_pol(vect[0].len);
    add_fft(&(new[0]), mul_fft(&tmp1, vect[0], mat[0][0]), mul_fft(&tmp2, vect[1], mat[1][0]));
    add_fft(&(new[1]), mul_fft(&tmp1, vect[0], mat[0][1]), mul_fft(&tmp2, vect[1], mat[1][1]));
    free(tmp1.coeffs);
    free(tmp2.coeffs);
    for (int i = 0; i < new[0].len; i++)
    {
        res[0].coeffs[i] = new[0].coeffs[i];
        res[1].coeffs[i] = new[1].coeffs[i];
    }
    free(new[0].coeffs);
    free(new[1].coeffs);
}
