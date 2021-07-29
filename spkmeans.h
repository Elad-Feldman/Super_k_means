#include "matrix_op.h"
#include "kmeans.h"
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



#define  is_goal(string) strcmp(goal,string) == 0
typedef  struct spk_results{
    int k;
    double **mat; /*TODO better name */
    Eigen eigen;
}spk_results;

void print_verbose(char* string);

void my_assert(int  cond);

/***************  SPK functions ******************/
double abs_d(double x);
void find_ind_max_ele_off_diag(double** A, int N,int* I, int* J);
double* renormlized_vector(double* a, int n);
void renorm_matrix_rows(double** U, int n, double** T);
void create_adj_mat(double** observations, int N, int D,double** W );
void create_diagonal_degree_mat_ns(double** adj_mat,int N,double** W);
spk_results activate_flag(char* goal,double** observations , int k, int n, int d);
/*TODO   add decleration for all functions here */

#ifndef FINAL_PROJECT_SPKMEANS_H
#define FINAL_PROJECT_SPKMEANS_H




#endif //FINAL_PROJECT_SPKMEANS_H
