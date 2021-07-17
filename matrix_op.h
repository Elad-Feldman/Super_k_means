#include <string.h>

#define Arr_size(x)  (sizeof(x) / sizeof((x)[0]))

typedef  struct tuple2{
    int i;
    int j;
} tuple2;

//Matrix
double** create_matrix(int N,int d);
void free_matrix(double** A,int N);
double** transpose_mat(double** mat, int N, int D);
double dot_vector(double a[], double b[], int n);
double** mult_matrix(double** A, double** B, int N);
double find_max_ele_off_diag(A double[][], int N);
int is_diagonal( double** A, int N);
double find_vec_norm(a double[], int n);
double* renorm_vec(double* a, int n);
double** renorm_matrix_cols(double** U, int N);
ouble find_vec_norm_diff(double* a, double* b, int n);

/*
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

