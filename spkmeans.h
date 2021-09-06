
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>


#define  is_goal(string) strcmp(goal,string) == 0
#define Ver 0     /* TODO  to zero before submitting */
#define print_verbose(x) if(Ver && printf(x)){}
#define my_assert(cond) assert( (cond) && "An Error Has Occured" )
#define isNull(x) ( x == NULL )
#define assert_not_null(x) super_assert( !isNull(x) ) /* if x is a null - this is an error */
#define assert_positive(x) super_assert( (x > 0) ) /* if x is a null - this is an error */

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
    double **T;
    int T_size;
}spk_results;
/****** Strcuts END *******/

/****** small function START *******/

double fix_neg_zero(double num);
void super_assert(int cond);
int assert_goal(char* goal);

/****** small function  END *******/


/*************** Vectors  START ******************/
void print_vector_int(int* a,int n);
void print_vector(double* a,int n);
double dot_mult_vector(double *a, double *b, int n);
double find_vec_norm(double* a, int n) ;
double find_vec_norm_diff(double* a, double* b, int n);
double sum_vector(double* a, int n) ;
void swap_int(int *a,int *b);
void swap_double(double *a, double *b);
void swap_double_pointers(double **a, double **b);
double* renormlized_vector(double* a, int n) ;
/*************** Vectors  END ******************/

/*************** Matrix  START ******************/
double **create_matrix(int rows, int cols);
void free_matrix( double  ** A, int rows);
void inplace_transpose_mat(double** mat, int rows, int cols);
double** transpose_mat(double** mat, int rows, int cols);
void mult_matrix(double** A, double** B, double ** C ,int n);
void copy_matrix(double** A, double** B ,int n,int m);
void sub_matrix(double** A, double** B, double** C,int n) ;
void print_mat( double  ** mat, int n, int d);
double** create_Id_matrix(int n);
void re_order_matrix_by_indces(double** A,int* indces, int n);
void renorm_matrix_rows(double** U, int n,int k, double** T);
/*************** Matrix END ******************/


/*************** kmean  START ******************/
double find_distance(double *dot, double *center, int d);
int get_index_of_closest_cluster(double* dot, double** cluster_list, int d, int k );
void update_cluster_center(double* dot, double * center,int cluster_size,int d,int sign) ;
double** get_init_clusters_list(double** T,int k);
int * init_clusters_indexes(int k);
void simple_kmean (double ** T_mat, double ** T_cluster_list, int* cluster_index_list, int n, int k, int d,int is_Py_call);
/*************** kmean  END ******************/

/************    Lnorm  START    ********/
void create_adj_mat(double** observations, int n, int d,double** W);
void create_diagonal_degree_mat(double** adj_mat,int n,double** D);
void D_sqrt(double** D,int n) ;
void create_L_norm(double** D, double** W, int n, double** L );
/************    Lnorm  END    ********/


/************    Jacobi  START    ********/
double abs_d(double x);
void find_ind_max_ele_off_diag(double** A, int n,int* I, int* J);
double sign(double x) ;
double calc_theta(double a_jj ,double a_ii,double a_ij) ;
double calc_t(double theta);
double calc_c(double t) ;
void update_V(int i, int j, int  n, double c, double s,double** V);
void update_A_f(int i, int j, int n, double c, double s, double** A, double ** A1);
double sum_square_elements_off_diag(double ** A,int n);
int check_convergence(double** A,double** A1,int n);
Eigen find_eigen_vectors_and_values(double** L, int n);
/************    Jacobi  END    ********/


/**** Mergesort start ****/
void merge(double* arr,int* arr_i, int l, int m, int r);
void mergeSort(double* arr,int* arr_i, int l, int r);
/**** Mergesort end ****/

/**********   Eigen_values START **********/
double* extract_eigen_values_from_mat(double** mat,int n);
int eigengap_huristic(Eigen eigen,int k);
void eigen_to_matrix(Eigen eigen, double** E,int n);
void free_eigen(Eigen eigen);
void print_eigen(Eigen eigen);
/**********   Eigen_values END **********/

/***************  SPK START ******************/

void start_wam(double** observations , int n, int d, double*** W);
void start_ddg(double** observations , int n, int d, double*** D);
void start_lnorm(double** observations , int n, int d, double*** L);
void start_jacobi(double** observations , int n, double*** E);
spk_results activate_flag(char* goal,double** observations , int k, int n, int d);
/***************  SPK END ******************/

/******** C Interface ******/
void load_string(char** str,char* cpy);
int string_to_doubles(char *row,double* arr);
Tuple2 load_observations_from_file(double** observations, char* file_name);
int main(int argc, char* argv[]);