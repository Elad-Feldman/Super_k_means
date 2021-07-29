#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kmeans.h"


 double find_distance(double *dot, double *center, int d){
    double dis;
    int i;
    dis = 0;

    for ( i = 0; i < d; i++)
        dis += (dot[i] - center[i]) * (dot[i] - center[i]);
    return  dis;

}
 int get_index_of_closest_cluster(double* dot, double** cluster_list, int d, int k )
{
    int j;
    int i;
    double min_dis;
    double tmp_dis;
    j = 0;
    min_dis = find_distance(dot, cluster_list[0], d);

    for (i = 1; i < k; i++)
    {
        tmp_dis = find_distance(dot, cluster_list[i], d);
        if (tmp_dis <= min_dis)
        {
            min_dis = tmp_dis;
            j = i;
        }
    }
    return j;
}
 void update_cluster_center(double* dot, double * center,int cluster_size,int d,int sign) {
    double* center_temp;
    int i;
    if (cluster_size+sign==0)
        printf("error \n ");

    center_temp  = (double *) calloc(d,sizeof(double ));
    for (i = 0; i < d; i++)
        center_temp[i] = (center[i] * (cluster_size));

    for (i = 0; i < d; i++){
        center_temp[i] += (dot[i]*sign);
        center[i] = center_temp[i] / (cluster_size+sign);

    }
    free(center_temp);
}


void simple_kmean (double ** T_mat, double ** T_cluster_list, int * cluster_index_list,double ** observations, int n, int k, int d) {
    int i,j;
    int *T_at;
    int *move_T_to;
    int *T_cluster_size;
    int max_iter;
    max_iter = 300;

    int is_a_cluster_changed;
    int count_iter;

    T_at = (int*) calloc(n, sizeof (int));
    move_T_to = (int*) calloc(n, sizeof (int));
    T_cluster_size = (int*) calloc(k, sizeof (int));


    for (i = 0; i < n; i++) {
        T_at[i] = -1;
        move_T_to[i] = 0;
    }

    for (i = 0; i < k; i++) { /* set inial  locations */
        j = (int) cluster_index_list[i]; /*[ 44,56,73 ] */
        T_at[j] = i;
        T_cluster_size[i] = 1;

    }


    is_a_cluster_changed = 1;
    count_iter = 0;
    while (count_iter < max_iter && is_a_cluster_changed) {
        int i, j;
        is_a_cluster_changed = 0;
        count_iter++;

        for (i = 0; i < n; i++) /*find nearest clusters */
            move_T_to[i] = get_index_of_closest_cluster(T_mat[i], T_cluster_list, d, k);

        for (j = 0; j < n; j++) {/* update clusters*/
            if (T_at[j] == -1) {
                T_at[j] = move_T_to[j];
                update_cluster_center(T_mat[j], T_cluster_list[move_T_to[j]], T_cluster_size[move_T_to[j]], d, 1); /*add dot to center*/
                T_cluster_size[move_T_to[j]]++;
                is_a_cluster_changed = 1;
            } else {
                if (T_at[j] != move_T_to[j]) {
                    update_cluster_center(T_mat[j], T_cluster_list[T_at[j]], T_cluster_size[T_at[j]], d,
                                          -1); /*remove dot from center */
                    update_cluster_center(T_mat[j], T_cluster_list[move_T_to[j]], T_cluster_size[move_T_to[j]], d,
                                          1); /*add dot to center */
                    T_cluster_size[T_at[j]]--;
                    T_cluster_size[move_T_to[j]]++;
                    T_at[j] = move_T_to[j];
                    is_a_cluster_changed = 1;
                }
            }
        }
    }

    /**** create the cluster for observations ****/
    double **Ob_clusters = create_matrix(n,d);
    int *Ob_cluster_size = (int*) calloc(k, sizeof (int));

    for (i = 0; i < n; i++) {/* update clusters*/
        j = T_at[i]; /*[ 44,56,73 ] */
        update_cluster_center(observations[i], Ob_clusters[j], Ob_cluster_size[j], d, 1); /*add dot to center*/
        T_cluster_size[i]++;
    }
    print_vector((double*)cluster_index_list,k);
    print_mat(Ob_clusters,n,d);
    
    free_matrix(Ob_clusters,n);
    free(Ob_cluster_size);
    free(T_at);
    free(move_T_to);
    free(T_cluster_size);

}