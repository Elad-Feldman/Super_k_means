#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include<math.h>


#define Arr_size(x)  (sizeof(x) / sizeof((x)[0]))

typedef  struct tuple2{
    int i;
    int j;
} tuple2;

/*************** Matrix functions ******************/
double** create_matrix(int N,int d);
void free_matrix(double** A,int N);
double** transpose_mat(double** mat, int N, int D);
double** mult_matrix(double** A, double** B, int N);
int* find_max_ele_off_diag(double** A , int N);
int is_diagonal( double** A, int N);
double** create_renorm_matrix_rows(double** U, int N);
double** build_adj_mat(double** dot_list, int N, int D);
/*************** Vectors functions ******************/
double* renormlized_vector(double* a, int n);
double find_vec_norm(double *a, int n);
double dot_mult_vector(double *a, double *b, int n);
double find_vec_norm_diff(double* a, double* b, int n);

/*************** implementation ******************/
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
    printf("==============================\n");
    for (i = 0; i < N; i++) {
        printf("[");
        for (j = 0; j < d; j++)
            printf("%.4f,", mat[i][j]);
        printf("]\n");
    }
    printf("==============================\n");
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

 double** create_Id_matrix(int N) {
     assert(N > 0);
     double** mat = create_matrix(N,N);
    int i;
     for (i = 0; i < N; i++) {
         mat[i][i] = 1; /*  for i=\=j ,calloc  allocated memory block to zero */

     }
     return mat;
 }

 double** transpose_mat(double** mat, int N, int D)
 {
     assert(mat);
     assert(N > 0);
     assert(D > 0);
     double** mat_T = create_matrix(D, N);
     int i, j;
     for (i = 0; i < N; i++)
     {
         for (j = 0; j < D; j++)
         {
             mat_T[j][i] = mat[i][j];
         }
     }
     return mat_T;
 }

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
 

 double** mult_matrix(double** A, double** B,int N) {
     assert(A);
     assert(B);
     assert(N > 0);
     double** mat = create_matrix(N, N);
     double** B_T = transpose_mat(B, N, N);
     int i, j;
     for (i = 0; i < N; i++)
     {
         for (j = 0; j < N; j++)
         {
             mat[i][j] = dot_mult_vector(A[i], B_T[j], N);
         }
     }
     free(B_T);
     return mat;
 }

 int* find_ind_max_ele_off_diag(double** A, int N) { // TODO we need to change name, we find the indexs, not the max value
     assert(A);
     assert(N > 0);
     int i, j;
     int* indices= calloc(2,sizeof(int));
     double max = A[0][1]; 
     for (i = 0; i < N; i++)
     {

         for (j = 0; j < N; j++)
         {
             if (i != j && max < A[i][j]) {
                 max = A[i][j];
                 indices[0] = i;
                 indices[1] = j;

             }
         }
     }
     return indices;
 }

 int is_diagonal(double** A, int N)
 {
     assert(A);
     assert(N > 0);
     int i, j;
     for (i = 0; i < N; i++) 
     {
         for (j = 0; j < N; j++)
         {
             if (i != j && A[i][j] != 0) {
                 return 0;
             }
         }
     }
     return 1;
 }



 double* renormlized_vector(double* a, int n) {
     assert(a);
     assert(n > 0);
     int i;
     double norm = find_vec_norm(a, n);
         double* norm_vec = calloc(n, sizeof(double));
     if (norm == 0) {
         return norm_vec; //TODO good idea to check, but should use  assert(norm==0); this should not happen;
     }
     for (i = 0; i < n; i++) {
         norm_vec[i] = (a[i] / norm);
     }
     return norm_vec;
 }

 double** create_renorm_matrix_rows(double** U, int N) {
     double** renorm_mat = create_matrix(N, N);
     int i,j;
     for (i = 0; i < N; i++)
       renorm_mat[i]= renormlized_vector(U[i], N);

     return renorm_mat;
 
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

//returns the euclidian distance bitween two vectors a, b
 double find_vec_norm_diff(double* a, double* b, int n) {
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

 double** build_adj_mat(double** dot_list, int N, int D){
     assert(N > 0);
     assert(D > 0);
     assert(dot_list);
     int i, j;
     double** mat = create_matrix(N, N);
     for (i = 0; i < N; i++)
     {
         for (j = 0; j < N; j++)
         {
             if(i!=j){//elements on the diagon are 0 because of calloc used in create matrix
                 mat[i][j] = exp((-find_vec_norm_diff(dot_list[i], dot_list[j], D)) / 2);
             }
         }
     }
     return mat;
 }

 //ns stands for negtive squre root this function returns D^(-1/2)
 double** create_diagonal_degree_mat_ns(double** adj_mat,int N) {
     assert(adj_mat);
     assert(N>0);
     double** mat = create_matrix(N,N);
     int i,j;
     double sum_of_weights;//sum of weights per row
     sum_of_weights = 0;
     for (i = 0; i < N; i++)
     {
         for (j = 0; j < length; j++)
         {
             sum_of_weights+= adj_mat[i][j];
         }
         assert(sum_of_weights > 0);
         mat[i][i] = sqrt(1 / sum_of_weights);
         sum_of_weights = 0;
     }
     return mat;
 }
 // subtracts B from A and returns the result
 double** matrix_subtraction(double** A, double** B,int N) {
     double** result_mat = create_matrix(N, N);
     int i, j;
     for ( i = 0; i < N; i++)
     {

         for (j = 0; j < N; j++)
         {
             result_mat[i][j] = A[i][j] - B[i][j];
         }
     }
     return result_mat;
 }
 
 calculate_L_norm(double** D,double** W,int N) {
     double** id_mat = create_Id_matrix(N);
     double** mat = mult_matrix(D, W, N);
     mat = mult_matrix(mat, D);
     mat = matrix_subtraction(id_mat, mat, N);
     free(id_mat);
     return mat;
 }
 int sign(double x) {
     if (x<0)
     {
         return -1;
     }
     return 1;
 }
 double calc_teta(double jj ,double ii,double ij) {
     if (ij == 0) {
         return 0;
     }
     return (jj - ii) / (2 * ij);
 }
 double calc_c(double teta) {
    sign(teta)/(abs(teta)-sqrt((teta*teta)+1)
 }

 double** build_rotation_mat(double** A,int N) {
     double** id_mat = create_Id_matrix(N);
     int* max_el_ind = find_ind_max_ele_off_diag(A, N);
     int i = max_el_ind[0];
     int j = max_el_ind[1];
     double teta = calc_teta(A[j][j],A[i][i],A[i][j]);
     double t = calc_t(teta);
 
 }
 int test_mat_op() {
     int i, j;
     int N = 6;
     int D = 6;
     double** mat = create_matrix(N, D);
     //double** id_mat = create_Id_matrix(N);
    // print_mat(id_mat, N, N);

     int count;
     count = 1;
     for ( i = 0; i < N; i++)
     {

         for (j = 0; j < D; j++)
         {
             mat[i][j] = count;
             count++;
         }
     }

     //double** mat_renormed = create_renorm_matrix_rows(mat, N);
     print_mat(mat, N, N);
     printf("\n");
     double** adj_mat = build_adj_mat(mat, N, N);
     printf("%f \n",adj_mat[0][0]);
     print_mat(adj_mat, N, N);
     free(mat);
     free(adj_mat);
     //free(id_mat);
     //free(mat_renormed);

     return 0;
 }
 
/*
 //  TODO implement
    implement jacobi algorithim
// Vector

*/

