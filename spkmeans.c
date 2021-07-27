#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "spkmeans.h"
// gcc spkmeans.c && gcc  -o spkmeans spkmeans.c && spkmeans  5 wam  in1.txt

#define Arr_size(x)  (sizeof(x) / sizeof((x)[0]))
#define  is_goal(string) strcmp(goal,string) == 0


/************* implementation  SPK -  matrix  ***********/



void renorm_matrix_rows(double** U, int n, double** T)
{
    int i;
    for (i = 0; i < n; i++)
        T[i]= renormlized_vector(U[i], n);
}

void create_adj_mat(double** observations, int n, int d,double** W)
{
    assert(n > 0);
    assert(d > 0);
    assert(observations);
    double norm ;
    int i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++) /* matrix is symtric, W[j][i] = W[i][j];   */
        {
            norm = find_vec_norm_diff(observations[i], observations[j], d);
            W[i][j] = exp((-norm) / 2);
            W[j][i] = W[i][j];
        }
    }

}

void create_diagonal_degree_mat_ns(double** adj_mat,int n,double** D) {
    /* ns stands for negtive squre root this function returns D^(-1/2) */
    assert(adj_mat);
    assert(n>0);
    int i;
    double sum_of_weights;//sum of weights per row
    sum_of_weights = 0;
    for (i = 0; i < n; i++)
    {
        sum_of_weights = sum_vector(adj_mat[i],n);
        if (sum_of_weights==0)
            print_vector(adj_mat[i],n);
       // assert(sum_of_weights > 0);
        D[i][i] = 1 / sqrt( sum_of_weights);

    }

}

void create_L_norm(double** D, double** W, int n, double** L ) {
    double** id_mat = create_Id_matrix(n);
    double** DW = create_matrix(n,n);
    double** DWD = create_matrix(n,n);

    mult_matrix(D, W, DW, n);
    mult_matrix(DW, D,DWD,n);
    sub_matrix(id_mat, DWD,L, n);

    free_matrix(DW,n);
    free_matrix(DWD,n);
    free_matrix(id_mat,n);

}

/************  implementation SPK - Rotation matrix    ********/
double abs_d(double x){
    if(x<0){
        return -x;
    }
    return x;
}

void find_ind_max_ele_off_diag(double** A, int n,int* I, int* J)
{
    assert(A);
    assert(n > 0);
    int i, j;
    double max_abs = -1;
    double tmp;

    for (i = 0; i < n; i++)
    {

        for (j = 0; j < n; j++)
        {
            tmp = abs_d(A[i][j]);
            if (( i != j ) && (max_abs < tmp )) {
                         max_abs = tmp;
                *I = i;
                *J = j;

            }
        }
    }
}

int sign(double x) {
    if (x<0)
    {
        return -1;
    }
    return 1;
}
double calc_theta(double A_jj ,double A_ii,double A_ij) {
    //printf("A_jj:%.4f, A_ii:%.4f, A_jj:%.4f, A_ij:%.4f",A_jj , A_ii, A_ij);
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
double** create_rotation_mat(double** A,int n)
{
    assert(A);
    int i,j;
    double c,t,s, theta;
    double** P = create_Id_matrix(n);
    find_ind_max_ele_off_diag(A, n,&i,&j);
    theta = calc_theta(A[j][j],A[i][i],A[i][j]);
    t = calc_t(theta);
    c = calc_c(t);
    s = t * c;
    //printf("c: %f s: %f t: %f theta: %f \n",c,s,t,theta);
    P[i][i] = c;
    P[j][j] = c;
    P[i][j] = s;
    P[j][i] = -1*s;
    return P;
}

double sum_square_elements_off_diag(double ** A,int n){
    int i,j;
    double sum=0;
    for (i = 0; i <n ;i++) {
        for (j = 0; j <n ;j++) {
            if(i!=j){
                sum+=A[i][j]*A[i][j];
            }

        }
    }
    return sum;
}
int check_convergence(double** A,double** A1,int n){
    double sum_A =sum_square_elements_off_diag(A, n);
    double sum_A1 = sum_square_elements_off_diag(A1,n);
    if(sum_A-sum_A1<=0.001){
        return 1;
    }else{
        return 0;
    }
}


/**********   implementation eigen_values **********/

double* extract_eigen_values_from_mat(double** mat,int n){
    int i;
    double * eigen_values = calloc(n,sizeof(double));
    assert(eigen_values);
    for(i=0;i<n;i++){
        eigen_values[i]=mat[i][i];
    }
    return eigen_values;
}
int partition (double* e_values,int* ranks, int low, int high)
{
     double pivot = e_values[high];
    int i = (low - 1);
    int j;
    for (j = low; j <= high- 1; j++)
    {
        if (e_values[j] <= pivot)
        {
            i++;
            swap_double(&e_values[i], &e_values[j] );
            swap_int(&ranks[i], &ranks[j] );
        }
    }
    swap_double(&e_values[i+1], &e_values[high] );
    swap_int(&ranks[i+1], &ranks[high]);
    return (i + 1);
}



void Qsort_eigen_values(double* e_values,int* ranks,int low, int high){
    /* based on https://hackr.io/blog/quick-sort-in-c */
   // print_vector(e_values,10);
    if (low < high) {
            int pi = partition(e_values, ranks, low, high);
            Qsort_eigen_values(e_values, ranks, low, pi - 1);
            Qsort_eigen_values(e_values,ranks, pi + 1, high);
        }


}

Eigen find_eigen_vectors_and_values(double** L, int n){
    /* Start with A = L_norm */
    double ** A = create_matrix(n,n);
    copy_matrix(A,L,n);
    print_verbose("start: find eigen vactors");
    assert(A);
    assert(n>0);
    int i;
    double** V = create_Id_matrix(n);
    double** V_tmp =  create_matrix(n,n);

    double **P,**P_T;
    double** A_tmp = create_matrix(n,n);
    double** A_f = create_matrix(n,n);

    int max_iter,convergence;
    max_iter = 100;
    convergence = 0;
    i = 0;
    Eigen eigen;

    while( !convergence ||  i<max_iter){
        i++;

        P = create_rotation_mat(A,n);
        P_T = transpose_mat(P,n,n);


        mult_matrix(P_T,A,A_tmp,n); /* A_tmp = pT*A  */
        mult_matrix(A_tmp,P,A_f,n); /* A_f = pT*A*p  */

        convergence = check_convergence(A,A_f,n);
        copy_matrix(A,A_f,n);
        mult_matrix(V,P,V_tmp,n);
        copy_matrix(V,V_tmp,n);
        free_matrix(P,n);
        free_matrix(P_T,n);
    }
    free_matrix(A_tmp,n);
    free_matrix(A,n);
    //free_matrix(V,n);
    print_verbose("\nfound vectors!\n");

    eigen.vectors = V;
    eigen.mat_size = n;
    eigen.values = extract_eigen_values_from_mat(A, n);
    eigen.ranks =  (int*) calloc(n , sizeof (int));
    assert(eigen.ranks);

    for (i = 0; i < n; i++) /* after sorting, in [i]=j, j would the be the rank of the i vector */
        eigen.ranks[i] = i;

    Qsort_eigen_values(eigen.values,eigen.ranks,0,n-1);
    re_order_matrix_by_indces(eigen.vectors, eigen.ranks, n); /* TODO are the rows the eigenvectors or the colmuns ? */
    return eigen;
}

int  eigengap_huristic(Eigen eigen){
    int k,m,i;
    m =  (int ) eigen.mat_size / 2; /* it said that k<n/2 */
    double  delta_i ;
    double max = - 1 ;
    k = 0;
    for (i=0 ; i <= m ; i++){
        delta_i = abs_d (eigen.values[i] - eigen.values[i+1]) ;
        if (delta_i > max)
            k=i+1;

    }
    assert(k>0);
    return k;
}


double* renormlized_vector(double* a, int n) {
    assert(a);
    assert(n > 0);
    int i;
    double norm = find_vec_norm(a, n);
    double* norm_vec = calloc(n, sizeof(double));
    assert(norm_vec);
    if (norm == 0) {
        return norm_vec;
    }
    for (i = 0; i < n; i++) {
        norm_vec[i] = (a[i] / norm);
    }
    return norm_vec;
}

void free_eigen(Eigen eigen){
    free_matrix(eigen.vectors,eigen.mat_size);
    free(eigen.ranks);
    free(eigen.values);
}


void test_sort()

{
    int i,n;
    n= 10;
    double * val  = (double *) malloc(n * sizeof (double ));
    int* ranks  =(int*) malloc(n * sizeof (int));
    for (i = 0; i < n; i++)
    {
        ranks[i] = i;
        val[i] = n/2 - i;
    }
    print_vector(val,n);
    Qsort_eigen_values(val,ranks,0,n-1);
    print_vector(val,n);


}



int test_mat_op() {
    int i, j;
    int n = 3;
    int D = 3;
    double** mat = create_matrix(n, D);
    int count;
    count = 1;
    for ( i = 0; i < n; i++)
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

    print_mat(mat, n, n);
    printf("\n");
    Eigen eo = find_eigen_vectors_and_values(mat, n);
    printf("eigen vectors: \n");
    print_mat(eo.vectors, n, n);
    printf("eigen values: \n");
    for (int k = 0; k < n; k++) {
        printf("%f \n", eo.values[k]);
    }
    free(mat);
    /* free(id_mat); */
    /* (mat_renormed); */

    return 0;
}


 /******** C Interface ******/

int assert_goal(char* goal)
{
   if (strcmp(goal,"spk")==0)
       return 1;
    if (strcmp(goal,"wam")==0)
        return 1;
    if (strcmp(goal,"ddg")==0)
        return 1;
    if (strcmp(goal,"lnorm")==0)
        return 1;
    if (strcmp(goal,"jacobi")==0)
        return 1;
    assert(!"goal is not define"); /* TODO invalid input */
    return 0;



}

void load_string(char** str,char* cpy)
{
    int len;
    len = (int)strlen(cpy);
    *str = (char  *) malloc(len * sizeof(char));
    assert(*str);
    strcpy(*str, cpy);
}

int string_to_doubles(char *row,double* arr)
{
    int i;
    char* ptr;
    i = 0;
    ptr = strtok(row, ",");
    while(ptr != NULL)
    {
        arr[i] = atof(ptr);
        ptr = 	strtok(NULL, ",");
        i++;
    //printf("we got %d features",i);

    }
    return i;


}
Tuple2 load_observations_from_file(double** observations, char* file_name)
{
    int d,n,i;
    FILE *fp;
    fp = fopen(file_name,"r");
    char* row = (char* )calloc(1000, sizeof(char));
    assert(row);
    i =0;
    while(fscanf(fp,"%s",row)==1){ // load data
        d = string_to_doubles(row, observations[i]);
        i++;
    }
    n=i;

    // change the size of observations to match the file
    observations = (double **) realloc(observations,n* sizeof(double *));
    assert(observations);
    for (i = 0; i < n; i++)
    {
        observations[i] = (double *) realloc(observations[i],d* sizeof(double));
        assert(observations[i]);
    }
    fclose(fp);
    free(row);
    Tuple2 sizes;
    sizes.i = n;
    sizes.j = d;
    return sizes;

}

void test_swap(){
    int n = 10;
    double *a = malloc(n * sizeof (double));
    assert(a);
    int i;
    for (i = 0; i < n; i++)
        a[i] =(double) i ;


    print_vector(a,n);

    free(a);
}

spk_results activate_flag(char* goal,double** observations , int k, int n, int d)
{
    /* run all the flags, that are not spk */
    spk_results Res;
    Res.k = k;
    Res.eigen.values = NULL;
    Res.eigen.vectors = NULL;
    Res.eigen.mat_size = 0;


    double** W = create_matrix(n, n);
    create_adj_mat(observations,n,d,W);

    if (is_goal("wam")){
        printf("wam:\n");
        Res.mat = W;
        print_mat(Res.mat,n,n);
        free_matrix(W,n);
        printf("\n");
        return Res;
    }

    double** D = create_matrix(n, n);
    create_diagonal_degree_mat_ns(W,n,D);
    if (is_goal("ddg"))
    {
        printf("ddg:\n");
        Res.mat = D;
        free_matrix(W,n);
        print_mat(Res.mat,n,n);
        free_matrix(D,n);//
        printf("\n");
        return Res;
    }


    double** L = create_matrix(n, n);
    create_L_norm(D,W,n,L);
    if (is_goal("lnorm"))
    {
        printf("lnorm:\n");
        Res.mat = L;
        free_matrix(W,n);
        free_matrix(D,n);
        print_mat(Res.mat,n,n);
        free_matrix(L,n);
        printf("\n");
        return Res;
    }

    //print_verbose("start jacobi\n");
    Res.eigen = find_eigen_vectors_and_values(L, n);
    if (k==0)
        k = eigengap_huristic(Res.eigen);

    if (is_goal("jacobi"))
    {
        printf("jacobi:\n");
        printf("eigen vectors:\n");
        print_mat(Res.eigen.vectors,n,k);
        printf("eigen values:\n");
        print_vector(Res.eigen.values,n);
        free_matrix(W,n);
        free_matrix(D,n);
        free_matrix(L,n);
        free_eigen(Res.eigen);
        print_verbose("finish jacobi");
        return Res;
    }


    int i;
    double** U  = (double**)  malloc(n* sizeof (double *)) ;
    assert(U);
    double** T = create_matrix(n,k);
    for (i=0; i<n; i++)
        U[i] =Res.eigen.vectors[i];
    renorm_matrix_rows(U, n, T);



    free_matrix(W,n);
    free_matrix(D,n);
    free_matrix(L,n);
    free(U); /* does not hold new vectors, just point to eigen vectors */
    Res.mat = T;
    Res.k = k;
    //print_mat(Res.eigen.vectors,Res.eigen.mat_size,Res.eigen.mat_size);

   free_eigen(Res.eigen);
    return Res;
}


double** init_clusters_list(double** T,int n,int k){
    int i;
    double** cluster_list = create_matrix(k,k);
    copy_matrix(cluster_list,T,k);
    return cluster_list;
}
init_clusters_indexes(int k){
    int * clusters_indexes = calloc(k,sizeof(int));
    int i;
    for(i=0;i<k;i++){
      clusters_indexes[i]=i;
    }
    return clusters_indexes;
}

int main(int argc, char* argv[])
{
    int k,n,d;
    n = 1000;
    d = 10;
    char*  goal;
    char*  file_name;
    double** observations = create_matrix(n,d);
    assert(argc==4);
      //TODO assert k is an integer
        k = atoi(argv[1]);
        load_string(&goal,argv[2]);
        load_string(&file_name,argv[3]);


    assert(k>=0);

    assert_goal(goal);
    printf("k: %d\ngoal %s\nfile_name: %s\n",k,goal,file_name);
    Tuple2 sizes = load_observations_from_file(observations, file_name);
    n=sizes.i;
    d=sizes.j;


    spk_results Res;
    Res = activate_flag( goal, observations , k,  n, d);
    print_verbose("\nfinish activate_flag\n");
    if(!is_goal("spk")){
        return;
    }
    k=Res.k;
    printf("found k: %d",Res.k);
    printf("create full spk here");
    double** T_clusters_list = init_clusters_list(Res.mat,n,k);
    int * T_clusters_indexes = init_clusters_indexes(k);
    simple_kmean(Res.mat, T_clusters_list, T_clusters_indexes,observations,n,k,d);

    //Free all
    free(goal);
    free(file_name);
    free_matrix(observations, n);
    printf("\n FREE all !");

}

/*TODO:
check why free eigen dosnt work
lines 551,529
*/