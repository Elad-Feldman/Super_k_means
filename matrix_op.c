#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include<math.h>
#include "matrix_op.h"

#define Arr_size(x)  (sizeof(x) / sizeof((x)[0]))


/*************** implementation Vector ******************/

double dot_mult_vector(double *a, double *b, int n) {
    assert(a);
    assert(b);
    assert(n > 0);
    double sum;
    int i;
    sum = 0;
    for (i = 0; i < n; i++)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

double find_vec_norm(double* a, int n) {
    int i;
    double sum;
    assert(a);
    assert(n > 0);
    sum = 0;
    for ( i = 0; i < n; i++)
    {
        sum += (a[i] * a[i]);
    }
    return sqrt(sum);
}

double find_vec_norm_diff(double* a, double* b, int n) {
     /* returns the euclidian distance bitween two vectors a, b */
    int i;
    double norm;
    double * c =  (double  *) calloc(n,sizeof(double));
    for (i = 0; i < n; i++)
        c[i]= a[i]-b[i];
    norm = find_vec_norm(c,n);
    free(c);
    return norm;
}

double sum_vector(double a[], int n) {
    int sum, i;
    sum = 0;
    for (i = 0; i < n; i++) {
        sum += a[i];
    }
    return sum;
}

void print_vector(double* a,int n){
    int i;
    printf("[");
    for (i = 0; i < n; i++)
        printf("%.0f,",a[i]);
    printf("]\n");
}


/*************** implementation MATRIX ******************/
double **create_matrix(int n, int d)
{
    double** mat = (double  **)calloc(n , sizeof(double*));
    assert(mat);
    int i;
    for (i = 0; i < n; i++) {
        mat[i] = (double  *) calloc(d,sizeof(double));
        assert(mat[i]);

    }
    return mat;

}

void free_matrix( double  ** A, int n)
{
    int i;
    for (i = 0; i < n; i++)
        free(A[i]);
    free(A);
}

double** transpose_mat(double** mat, int n, int D)
{
    assert(mat);
    assert(n > 0);
    assert(D > 0);
    double** mat_T = create_matrix(D, n);
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < D; j++)
        {
            mat_T[j][i] = mat[i][j];
        }
    }
    return mat_T;
}

void mult_matrix(double** A, double** B, double ** C ,int n) {
    assert(A);
    assert(B);
    assert(C);
    assert(n > 0);

    double** B_T = transpose_mat(B, n, n);
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            C[i][j] = dot_mult_vector(A[i], B_T[j], n);
        }
    }
    free_matrix(B_T,n);
}

void sub_matrix(double** A, double** B, double** C,int n) {
    int i, j;
    for ( i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
    }

}

void print_mat( double  ** mat, int n, int d)
{
    int i, j;
    printf("==============================\n");
    for (i = 0; i < n; i++) {
        printf("[");
        for (j = 0; j < d; j++)
            printf("%.4f,", mat[i][j]);
        printf("]\n");
    }
    printf("==============================\n");
}

double** create_Id_matrix(int n) {
     assert(n > 0);
     double** mat = create_matrix(n,n);
    int i;
     for (i = 0; i < n; i++) {
         mat[i][i] = 1; /*  for i=\=j ,calloc  allocated memory block to zero */

     }
     return mat;
 }


/***************  SPK functions ******************/




 

/*
   TODO implement
    1.to optimize adj_weighted_mat
    implement jacobi algorithim


*/

