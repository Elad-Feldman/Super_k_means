#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

static double find_distance(double *dot, double *center, int d){
    double dis;
    int i;
    dis = 0;

    for ( i = 0; i < d; i++)
        dis += (dot[i] - center[i]) * (dot[i] - center[i]);
    return  dis;

}
static int get_index_of_closest_cluster(double* dot, double** cluster_list, int d, int k )
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
static void update_cluster_center(double* dot, double * center,int cluster_size,int d,int sign) {
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


static  void simple_kmean (double ** dot_list, double ** cluster_list, double * cluster_index_list, int n,int k,int d) {
    int i,j;
    int *dot_at;
    int *move_dot_to;
    int *cluster_size;
    int max_iter;
    max_iter = 300;

    int is_a_cluster_changed;
    int count_iter;

    dot_at = (int*) calloc(n,sizeof (int));
    move_dot_to = (int*) calloc(n,sizeof (int));
    cluster_size = (int*) calloc(k,sizeof (int));


    for (i = 0; i < n; i++) {
        dot_at[i] = -1;
        move_dot_to[i] = 0;
    }

    for (i = 0; i < k; i++) {
        j = (int) cluster_index_list[i]; /*[ 44,56,73 ] */
        dot_at[j] = i;
        cluster_size[i] = 1;

    }


    is_a_cluster_changed = 1;
    count_iter = 0;
    while (count_iter < max_iter && is_a_cluster_changed) {
        int i, j;
        is_a_cluster_changed = 0;
        count_iter++;

        for (i = 0; i < n; i++) /*find nearest clusters */
            move_dot_to[i] = get_index_of_closest_cluster(dot_list[i], cluster_list, d, k);

        for (j = 0; j < n; j++) {/* update clusters*/
            if (dot_at[j] == -1) {
                dot_at[j] = move_dot_to[j];
                update_cluster_center(dot_list[j], cluster_list[move_dot_to[j]], cluster_size[move_dot_to[j]], d,1); /*add dot to center*/
                cluster_size[move_dot_to[j]]++;
                is_a_cluster_changed = 1;
            } else {
                if (dot_at[j] != move_dot_to[j]) {
                    update_cluster_center(dot_list[j], cluster_list[dot_at[j]], cluster_size[dot_at[j]], d,
                                          -1); /*remove dot from center */
                    update_cluster_center(dot_list[j], cluster_list[move_dot_to[j]], cluster_size[move_dot_to[j]], d,
                                          1); /*add dot to center */
                    cluster_size[dot_at[j]]--;
                    cluster_size[move_dot_to[j]]++;
                    dot_at[j] = move_dot_to[j];
                    is_a_cluster_changed = 1;
                }
            }
        }
    }

    /* print_matrix(cluster_list, k, d); */
    free(dot_at);
    free(move_dot_to);
    free(cluster_size);

}
