#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include<math.h>


#define Arr_size(x)  (sizeof(x) / sizeof((x)[0]))
typedef  struct tuple2{
    int i;
    int j;
}tuple2;


typedef  struct eigen_obejct{
    double ** mat;
    double* eigen_values;
} e_o;

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
 double abs_d(double x){
    if(x<0){
        return -x;
    }
     return x;
}

 int* find_ind_max_ele_off_diag(double** A, int N) { // TODO we need to change name, we find the indexs, not the max value
     assert(A);
     assert(N > 0);
     int i, j;
     int* indices= calloc(2,sizeof(int));
     double max_abs = abs_d(A[0][1]);
     for (i = 0; i < N; i++)
     {

         for (j = 0; j < N; j++)
         {
             if (i != j && max_abs < abs_d(A[i][j])) {
                 max_abs = abs_d(A[i][j]);
                 indices[0] = i;
                 indices[1] = j;

             }
         }
     }
     print_mat(A,N,N);
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
     int i;
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

 double** create_adj_mat(double** dot_list, int N, int D){
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
         sum_of_weights= sum_vector(adj_mat[i],N);
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
 
double** calculate_L_norm(double** D,double** W,int N) {
     double** id_mat = create_Id_matrix(N);
     double** mat = mult_matrix(D, W, N);
     mat = mult_matrix(mat, D,N);
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
 double calc_t(double teta) {
   return sign(teta)/(abs(teta)+sqrt((teta*teta)+1));
 }
double calc_c(double t) {
    return 1 / sqrt((t * t) + 1);
}
 double** create_rotation_mat(double** A,int N) {
     double** id_mat = create_Id_matrix(N);
     int* max_el_ind = find_ind_max_ele_off_diag(A, N);
     int i = max_el_ind[0];
     int j = max_el_ind[1];
     free(max_el_ind);
     double teta = calc_teta(A[j][j],A[i][i],A[i][j]);
     double t = calc_t(teta);
     double c = calc_c(t);
     id_mat[i][i] = c;
     id_mat[j][j] = c;
     id_mat[i][j] = c*t;
     id_mat[j][i] = -1*c*t;
      return id_mat;
 }
 double* extract_eigen_values_from_mat(double** mat,int N){
    int i;
    double * eigen_values = calloc(N,sizeof(double));
    for(i=0;i<N;i++){
        eigen_values[i]=mat[i][i];
    }
    return eigen_values;
}

 e_o find_eigen_vectors(double** Lnorm,int N){//
     assert(Lnorm);
     assert(N>0);
     double** A = Lnorm;
     double** V= create_Id_matrix(N);
     int count=0;
     while(!is_diagonal(A,N)&&count<100){
         count++;
         print_mat(A,N,N);
         printf("\n");
         double** P =create_rotation_mat(A,N);
         double** P_T = transpose_mat(P,N,N);
         A= mult_matrix(P_T,A,N);
         A= mult_matrix(A,P,N);/*P^T*A*P*/
         V =mult_matrix(V,P,N);
    }
     e_o eigen_values_and_vectors;
     eigen_values_and_vectors.eigen_values = extract_eigen_values_from_mat(A,N);
     eigen_values_and_vectors.mat = V;
     return eigen_values_and_vectors;
}
 int test_mat_op() {
     int i, j;
     int N = 3;
     int D = 3;
     double** mat = create_matrix(N, D);
     int count;
     count = 1;
     for ( i = 0; i < N; i++)
     {

         for (j = 0; j < D; j++)
         {
             if(mat[i][j]==0){
                 mat[i][j] = count;
                 mat[j][i] = count;
                 count++;
             }

         }
     }

     print_mat(mat, N, N);
     printf("\n");
     e_o eo = find_eigen_vectors(mat,N);
     printf("worked");
     print_mat(eo.mat, N, N);
     for (int k = 0; k < N; k++) {
         printf("%f \n",eo.eigen_values[k]);
     }
     free(mat);
     //free(id_mat);
     //free(mat_renormed);

     return 0;
 }
 
/*
 //  TODO implement
    1.to build a struct which contains matrix and array for eigen values and eigen vectors
    2.to optimize adj_weighted_mat
    3.to change diagonal degree to only one loop

    implement jacobi algorithim
// Vector

*/

