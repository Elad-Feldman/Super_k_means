#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include<math.h>
#include "matrix_op.h"

/*23.7 */
void free_matrix( double  ** A, int n)
{
    int i;
    for (i = 0; i < n; i++)
        free(A[i]);
    free(A);
}


void print_mat( double  ** mat, int n, int d)
{
    int i, j;
    for (i = 0; i < n; i++) {
        printf("[");
        for (j = 0; j < d; j++)
            printf("%.4f,", mat[i][j]);
        printf("]\n");
    }
}


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
 //creates a matrix such that every element on the diagon is 1 and all other elements are 0
 double** create_Id_matrix(int n) {
     assert(n > 0);
     double** mat = create_matrix(n,n);
    int i;
     for (i = 0; i < n; i++) {
         mat[i][i] = 1;
     }
     return mat;
 }
 //transpose matrix
 double** transpose_mat(double** mat, int n, int D) {
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
 //scalar multiplication
 double dot_vector(double a[], double b[], int n) {
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
 
 //vectorial multiplication
 double** mult_matrix(double** A, double** B,int n) {
     assert(A);
     assert(B);
     assert(n > 0);
     double** mat = create_matrix(n, n);
     double** B_T = transpose_mat(B, n, n);
     int i, j;
     for (i = 0; i < n; i++)
     {
         for (j = 0; j < n; j++)
         {
             mat[i][j] = dot_vector(A[i], B_T[j],n);
         }
     }
     return mat;
 }
 //finding the biggest element that is not on the diagon
 double find_max_ele_off_diag(double** A, int n) {
     assert(A);
     assert(n > 0);
     int i, j;
     double max = 0;

     for (i = 0; i < n; i++)
     {

         for (j = 0; j < n; j++)
         {
             if (i != j && max < A[i][j]) {
                 max = A[i][j];
             }
         }
     }
     return max;
 }
 //checks if matrix is diagonal
 int is_diagonal(double** A, int n) {
     assert(A);
     assert(n > 0);
     int i, j;
     for (i = 0; i < n; i++)
     {
         for (j = 0; j < n; j++)
         {
             if (i != j && A[i][j] != 0) {
                 return 0;
             }
         }
     }
     return 1;
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
//returns a normlized vector
 double* renorm_vec(double* a, int n) {
     assert(a);
     assert(n > 0);
     int i;
     double norm = find_vec_norm(a, n);
     double* norm_vec = calloc(n, sizeof(double));
     if (norm == 0) {
         return norm_vec;
     }
    
     for (i = 0; i < n; i++) {
         norm_vec[i] = (a[i] / norm);
     }
     return norm_vec;
 }
 //normalizes a matrix's columns
 double** renorm_matrix_cols(double** U, int n) {
     double** renorm_mat = create_matrix(n, n);
     int i,j;
     double* norm_vec;
     for (i = 0; i < n; i++)
     {
         norm_vec = renorm_vec(U[i], n);
         for (j = 0; j < n; j++) {
             renorm_mat[i][j] = norm_vec[j];
         }
     }
     free(norm_vec);
     return renorm_mat;
 
 }
//returns the euclidian distance bitween two vectors a, b
 double find_vec_norm_diff(double* a, double* b, int n) {
 
 }
 int main(int argc, char* argv[]) {
     int i, j;
     int n = 6;
     int D = 6;
     double** mat = create_matrix(n, D);
     double** id_mat = create_Id_matrix(n);
     int count;
     count = 1;
     for ( i = 0; i < n; i++)
     {

         for (j = 0; j < D; j++)
         {
             mat[i][j] = count;
             count++;
         }
     }
     double** mat_renormed = renorm_matrix_cols(mat,n);
     print_mat(mat, n, n);
     printf("\n");
     print_mat(mat_renormed, n, n);

     return 0;
 }
 
/*
 //  TODO implement

// Vector
double find_vec_norm_diff(a double [],b double [], int n);
double sum_vector (a double[],int n);
*/