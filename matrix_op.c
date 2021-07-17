#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "matrix_op.h"


void free_matrix( double  ** A, int N)
{
    int i;
    for (i = 0; i < N; i++)
        free(A[i]);
    free(A);
}


void print_mat( double  ** mat, int N, int d)
{
    int i, j;
    for (i = 0; i < N; i++) {
        printf("[");
        for (j = 0; j < d; j++)
            printf("%.4f,", mat[i][j]);
        printf("]\n");
    }
}


 double **create_matrix(int N, int d)
 {
    double** mat = (double  **)calloc(N , sizeof(double*));
     assert(mat);
    int i;
    for (i = 0; i < N; i++) {
        mat[i] = (double  *) calloc(d,sizeof(double));
        assert(mat[i]);

    }
    return mat;

}

/*
 //  TODO implement
double** create_Id_matrix(int N);
double** mult_matrix(A double[][],B double[][],int N);
int find_max_ele_off_diag(A double[][],int N); // int size 2, (i,j)
int is_diagonal(A double[][],int N); //
double** renorm_matrix_cols(U, int N);

// Vector
double mult_vec(a double[],a double[][],int n);
double find_vec_norm_diff(a double [],b double [], int n);
double find_vec_norm(a double [], int n);
double sum_vector (a double[],int n);
*/