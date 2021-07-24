
typedef  struct spk_results{
    int k;
    double **mat; /*TODO better name */
    Eigen eigen;
}spk_results;

/***************  SPK functions ******************/
double abs_d(double x);
Tuple2 find_ind_max_ele_off_diag(double** A, int N);
double* renormlized_vector(double* a, int n);
void renorm_matrix_rows(double** U, int n, double** T)
void create_adj_mat(double** observations, int N, int D,double** W );
void create_diagonal_degree_mat_ns(double** adj_mat,int N,double** W);

/*TODO   add decleration for all functions here */

#ifndef FINAL_PROJECT_SPKMEANS_H
#define FINAL_PROJECT_SPKMEANS_H




#endif //FINAL_PROJECT_SPKMEANS_H
