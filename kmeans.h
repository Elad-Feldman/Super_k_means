#include "matrix_op.h"
 double find_distance(double *dot, double *center, int d);
 int get_index_of_closest_cluster(double* dot, double** cluster_list, int d, int k );
 void update_cluster_center(double* dot, double * center,int cluster_size,int d,int sign);
 void simple_kmean (double ** T_mat, double ** T_cluster_list, int * cluster_index_list,double ** observations, int n, int k, int d);