
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#define  is_goal(string) strcmp(goal,string) == 0

#define verbose 0    \
/* TODO  to zero before submitting */

/****** Strcuts START *******/
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
typedef  struct spk_results{
    int k;
    double **mat; /*TODO better name */
    Eigen eigen;
}spk_results;
/****** Strcuts END *******/

/****** small function START *******/
void print_verbose(char* string);
void my_assert(int  cond);
int assert_goal(char* goal);
/****** small function START *******/



/*************** Vectors functions START ******************/
double dot_mult_vector(double *a, double *b, int n);
double find_vec_norm(double *a, int n);
double find_vec_norm_diff(double* a, double* b, int n);
double* renormlized_vector(double* a, int n);
double sum_vector(double* a, int n);
void print_vector(double* a,int n);
void swap_int(int *a,int *b);
void swap_double(double *a, double *b);
double* renormlized_vector(double* a, int n);
/*************** Vectors functions END ******************/


/*************** Matrix functions  START ******************/
double** create_matrix(int n,int d);
void free_matrix(double** A,int n);
double** transpose_mat(double** mat, int n, int D);
void mult_matrix(double** A, double** B, double ** C ,int n);
void copy_matrix(double** A, double** B ,int n);
void sub_matrix(double** A, double** B, double ** C, int n);
void print_mat( double  ** mat, int n, int d);
double** create_Id_matrix(int n);
void re_order_matrix_by_indces(double** A,int* indces, int n);
void renorm_matrix_rows(double** U, int n, double** T);
/*************** Matrix functions  END ******************/

/*************** kmean  START ******************/
static double find_distance(double *dot, double *center, int d);
static int get_index_of_closest_cluster(double* dot, double** cluster_list, int d, int k );
static void update_cluster_center(double* dot, double * center,int cluster_size,int d,int sign);
double** get_init_clusters_list(double** T,int k);
int * init_clusters_indexes(int k);
void simple_kmean (double ** T_mat, double ** T_cluster_list, int * cluster_index_list,double ** observations, int n, int k, int d);
/*************** kmean  END ******************/



/************    Lnorm  START    ********/
void create_adj_mat(double** observations, int n, int d,double** W);

void create_diagonal_degree_mat_ns(double** adj_mat,int n,double** D);

void create_L_norm(double** D, double** W, int n, double** L );
/************    Lnorm  END    ********/


/************    Jacobi  START    ********/
double abs_d(double x);
void find_ind_max_ele_off_diag(double** A, int n,int* I, int* J);
int sign(double x) ;
double calc_theta(double A_jj ,double A_ii,double A_ij);
double calc_t(double theta) ;
double calc_c(double t) ;
double** create_rotation_mat(double** A,int n);
double sum_square_elements_off_diag(double ** A,int n);
int check_convergence(double** A,double** A1,int n);
/************    Jacobi  END    ********/


/**********   Eigen_values START **********/
double* extract_eigen_values_from_mat(double** mat,int n);

int partition (double* e_values,int* ranks, int low, int high);

void Qsort_eigen_values(double* e_values,int* ranks,int low, int high);

Eigen find_eigen_vectors_and_values(double** L, int n);

int  eigengap_huristic(Eigen eigen);

void free_eigen(Eigen eigen);
/**********   Eigen_values END **********/

/***************  SPK START ******************/
spk_results activate_flag(char* goal,double** observations , int k, int n, int d);
/***************  SPK END ******************/

/******** C Interface ******/
void load_string(char** str,char* cpy);

int string_to_doubles(char *row,double* arr);
Tuple2 load_observations_from_file(double** observations, char* file_name);

int main(int argc, char* argv[]);