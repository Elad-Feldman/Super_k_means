#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "matrix_op.h"
static double find_distance(double *dot, double *center, int d);
static int get_index_of_closest_cluster(double* dot, double** cluster_list, int d, int k );
static void update_cluster_center(double* dot, double * center,int cluster_size,int d,int sign);
static  void simple_kmean (double ** T_mat, double ** T_cluster_list, double * cluster_index_list, int n, int k, int d);