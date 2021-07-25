#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include<math.h>

typedef  struct Tuple2{
    int i;
    int j;
}Tuple2;


typedef  struct Eigen_Obejct{
    double ** vectors;
    double* values;
    int *ranks;
    int mat_size;
} Eigen;

/*************** Vectors functions ******************/

double dot_mult_vector(double *a, double *b, int n);
double find_vec_norm(double *a, int n);
double find_vec_norm_diff(double* a, double* b, int n);
double* renormlized_vector(double* a, int n);
double sum_vector(double* a, int n);
void print_vector(double* a,int n);
void swap_int(int *a,int *b);
void swap_double(double *a, double *b);




/*************** Matrix functions ******************/
double** create_matrix(int n,int d);
void free_matrix(double** A,int n);
double** transpose_mat(double** mat, int n, int D);
void mult_matrix(double** A, double** B, double ** C ,int n);
void copy_matrix(double** A, double** B ,int n);
void sub_matrix(double** A, double** B, double ** C, int n);
void print_mat( double  ** mat, int n, int d);
double** create_Id_matrix(int n);

void re_order_matrix_by_indces(double** A,int* indces, int n);

