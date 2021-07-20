#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include<math.h>


#define Arr_size(x)  (sizeof(x) / sizeof((x)[0]))
typedef  struct Tuple2{
    int i;
    int j;
}Tuple2;


typedef  struct Eigen_Obejct{
    double ** vectors;
    int mat_size;
    double* values;
} Eigen;

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
     free(B_T); // TODO use free_matrix !!
     return mat;
 }
 double abs_d(double x){
    if(x<0){
        return -x;
    }
     return x;
}

Tuple2 find_ind_max_ele_off_diag(double** A, int N)
     { // TODO we need to change name, we find the indexs, not the max value
     assert(A);
     assert(N > 0);
     int i, j;
     Tuple2 max_indeces;
     double max_abs = abs_d(A[0][1]);
     for (i = 0; i < N; i++)
     {

         for (j = 0; j < N; j++)
         {
             if (i != j && max_abs < abs_d(A[i][j])) {
                 max_abs = abs_d(A[i][j]);
                 max_indeces.i = i;
                 max_indeces.j = j;

             }
         }
     }
     print_mat(A,N,N);
     return max_indeces;
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
         return norm_vec;
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

 double** create_adj_mat(double** observations, int N, int D)
 {
     assert(N > 0);
     assert(D > 0);
     assert(observations);
     double norm ;
     int i, j;
     double** W = create_matrix(N, N);
     for (i = 0; i < N; i++)
     {
         for (j = 0; j < i; j++) /* matrix is symtric, W[j][i] = W[i][j];   */
         {
             norm = find_vec_norm_diff(observations[i], observations[j], D);
             W[i][j] = exp((-norm) / 2);
             W[j][i] = W[i][j];
         }
     }
     return W;
 }

 //ns stands for negtive squre root this function returns D^(-1/2)
 double** create_diagonal_degree_mat_ns(double** adj_mat,int N) {
     assert(adj_mat);
     assert(N>0);
     double** mat = create_matrix(N,N);
     int i;
     double sum_of_weights;//sum of weights per row
     sum_of_weights = 0;
     for (i = 0; i < N; i++)
     {
         sum_of_weights= sum_vector(adj_mat[i],N);
         assert(sum_of_weights > 0);
         mat[i][i] = 1 / sqrt( sum_of_weights);

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
 
double** calculate_L_norm(double** D,double** W,int N) { /* TODO what happen to memory here? not clear */
    double** after_sub_mat;
    double** id_mat = create_Id_matrix(N);
    double** mat = mult_matrix(D, W, N);
    mat = mult_matrix(mat, D,N); // TODO  you didn't free the first mat !
    after_sub_mat = matrix_subtraction(id_mat, mat, N);
    free(mat); // TODO use free_matrix !!
    free(id_mat); // TODO use free_matrix !!
    return after_sub_mat;
 }
 int sign(double x) {
     if (x<0)
     {
         return -1;
     }
     return 1;
 }
 double calc_theta(double A_jj ,double A_ii,double A_ij) {
     if (A_ij == 0) {
         return 0;
     }
     return (A_jj - A_ii) / (2 * A_ij);
 }
 double calc_t(double theta) {
   return sign(theta)/(abs_d(theta)+sqrt((theta*theta)+1));
 }
double calc_c(double t) {
    return 1 / sqrt((t * t) + 1);
}
 double** create_rotation_mat(double** A,int N)
    {
     int i,j;
     double c,t,s, theta;
     double** P = create_Id_matrix(N);
     Tuple2 max_indces = find_ind_max_ele_off_diag(A, N);
     i = max_indces.i;
     j = max_indces.j;
     theta = calc_theta(A[j][j],A[i][i],A[i][j]);
     t = calc_t(theta);
     c = calc_c(t);
     s = t * c;
        printf("c: %f s: %f t: %f theta: %f \n",c,s,t,theta);
     P[i][i] = c;
     P[j][j] = c;
     P[i][j] = s;
     P[j][i] = -1*s;
     print_mat(P,N,N);
     return P;
    }
 double* extract_eigen_values_from_mat(double** mat,int N){
    int i;
    double * eigen_values = calloc(N,sizeof(double));
    for(i=0;i<N;i++){
        eigen_values[i]=mat[i][i];
    }
    return eigen_values;
}
double sum_square_elements_off_diag(double ** A,int N){
    int i,j;
    double sum=0;
    for (i = 0; i <N ;i++) {
        for (j = 0; j <N ;j++) {
            if(i!=j){
                sum+=A[i][j]*A[i][j];
            }

        }
    }
    return sum;
}
int check_convergence(double** A,double** A1,int N){
    double sum_A =sum_square_elements_off_diag(A, N);
    double sum_A1 = sum_square_elements_off_diag(A1,N);
    if(sum_A-sum_A1<=0.001){
        return 1;
    }else{
        return 0;
    }
}

Eigen find_eigen_vectors(double** A, int N){
     /* Start with A=L_norm */
     assert(A);
     assert(N>0);
     double** V = create_Id_matrix(N);
     double **P,** P_T;
     double** A1;
     double** A2; // TODO maybe A_tmp?
     double** V1; // TODO maybe V_tmp?
     int convergence = 0;
     Eigen eigen;
     while( !convergence){ // TODO I split he code
         /*print_mat(A,N,N);*/
         printf("\n");
         P = create_rotation_mat(A,N);
         P_T = transpose_mat(P,N,N);

         A1 = mult_matrix(P_T,A,N);
         A2 = mult_matrix(A1,P,N);
         free(A1); // TODO use free_matrix !!

         convergence= check_convergence(A,A2,N); // TODO we always use the same A, but we still re-calcualting off(A) each time ! maybe we an pass A2 the off(A) ?
         free(A); // Todo this is the value we are passing
         A=A2;
         print_mat(A,N,N);
         A2=NULL; // TODO why?

         V1 = mult_matrix(V,P,N);
         free(V); // TODO use free_matrix !!
         V=V1;
         V1=NULL; // TODO why?

         free_matrix(P,N);
         free_matrix(P_T,N);
    }
     eigen.values = extract_eigen_values_from_mat(A, N);
     eigen.vectors = V;
     eigen.mat_size = N;
     return eigen;
}
void sort_eigen_values(Eigen eign_obj){
    double** mat_T = transpose_mat(eign_obj.vectors, eign_obj.mat_size);

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
     Eigen eo = find_eigen_vectors(mat, N);
     printf("eigen vectors: \n");
     print_mat(eo.vectors, N, N);
     printf("eigen values: \n");
     for (int k = 0; k < N; k++) {
         printf("%f \n", eo.values[k]);
     }
     free(mat);
     /* free(id_mat); */
     /* (mat_renormed); */

     return 0;
 }
 
/*
   TODO implement
    1.to optimize adj_weighted_mat
    implement jacobi algorithim


*/

