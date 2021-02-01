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

void    print_fft(t_pol_fft f)
{
    printf("(");
    for (int i = 0; i < f.len; i++)
    {
        printf("%f%+fi", crealf(f.coeffs[i]), cimagf(f.coeffs[i]));
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
