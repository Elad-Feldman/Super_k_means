#include "matrix_op.h"
#include "spkmeans.h"
#include "kmeans.h"
// gcc spkmeans.c && gcc  -o spkmeans spkmeans.c && spkmeans  5 spk  in1.txt

#define Arr_size(x)  (sizeof(x) / sizeof((x)[0]))



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
        sum_of_weights= sum_vector(adj_mat[i],n);
        assert(sum_of_weights > 0);
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

Tuple2 find_ind_max_ele_off_diag(double** A, int n)
{ // TODO we need to change name, we find the indexs, not the max value
    assert(A);
    assert(n > 0);
    int i, j;
    Tuple2 max_indeces;
    double max_abs = abs_d(A[0][1]);
    for (i = 0; i < n; i++)
    {

        for (j = 0; j < n; j++)
        {
            if (i != j && max_abs < abs_d(A[i][j])) {
                max_abs = abs_d(A[i][j]);
                max_indeces.i = i;
                max_indeces.j = j;

            }
        }
    }
    print_mat(A,n,n);
    return max_indeces;
}

int sign(double x) {
    if (x<0)
    {
        return -1;
    }
    return 1;
}
double calc_theta(double A_jj ,double A_ii,double A_ij) {
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
    int i,j;
    double c,t,s, theta;
    double** P = create_Id_matrix(n);
    Tuple2 max_indces = find_ind_max_ele_off_diag(A, n);
    i = max_indces.i;
    j = max_indces.j;
    theta = calc_theta(A[j][j],A[i][i],A[i][j]);
    t = calc_t(theta);
    c = calc_c(t);
    s = t * c;
    printf("c: %f s: %f t: %f theta: %f \n",c,s,t,theta);
    P[i][i] = c;
    P[j][j] = c;
    P[i][j] = s;
    P[j][i] = -1*s;
    print_mat(P,n,n);
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
    for(i=0;i<n;i++){
        eigen_values[i]=mat[i][i];
    }
    return eigen_values;
}

void swap_array_elements(void * a, void * b, size_t len)
{
    /*  https://stackoverflow.com/questions/29596151/swap-function-using-void-pointers*/
    assert(a);
    assert(b);
    assert( a != b);
    unsigned char * p = a, * q = b, tmp;
    for (size_t i = 0; i != len; ++i)
    {
        tmp = p[i];
        p[i] = q[i];
        q[i] = tmp;
    }
}

void sort_eigen_values(Eigen eigen_obj,int first,int last){
    /* based on https://hackr.io/blog/quick-sort-in-c */
    double* values = eigen_obj.values;
    int * ranks = eigen_obj.values_index_sorted;
    int i, j, pivot;
    if( first < last) {
        pivot = first;
        i = first;
        j = last;

        while( i < j )
        {
            while( values[i] <= values[pivot] && i < last )
                i++;
            while( values[j] > values[pivot] )
                j--;
            if(i<j)
            {
                swap_array_elements(&values[i], &values[j], sizeof(values[i]) );
                swap_array_elements(&ranks[i], &ranks[j], sizeof(ranks[i]) );
            }
        }

        swap_array_elements(&values[pivot], &values[j], sizeof(values[i]) );
        swap_array_elements(&ranks[pivot], &ranks[j], sizeof(ranks[i]) );

        sort_eigen_values(eigen_obj,first,j-1);
        sort_eigen_values(eigen_obj,j+1,last);
    }


}

Eigen find_eigen_vectors_and_values(double** A, int n){
    /* Start with A = L_norm */
    assert(A);
    assert(n>0);
    int i;
    double** V = create_Id_matrix(n);
    double** V_tmp;

    double **P,**P_T;
    double** A_tmp;
    double** A_f;

    int convergence = 0;
    Eigen eigen;

    while( !convergence){
        P = create_rotation_mat(A,n);
        P_T = transpose_mat(P,n,n);

        /* A_f = P_T* A * P */
        mult_matrix(P_T,A,A_tmp,n);
        mult_matrix(A_tmp,P,A_f,n);

        convergence = check_convergence(A,A_f,n);
        free_matrix(A_tmp,n);
        free_matrix(A,n);

        A = A_f;
        A_f = NULL;

        mult_matrix(V,P,V_tmp,n);
        free_matrix(V,n);
        V = V_tmp;
        V_tmp = NULL;

        free_matrix(P,n);
        free_matrix(P_T,n);
    }
    eigen.values = extract_eigen_values_from_mat(A, n);
    eigen.values_index_sorted = malloc(n* sizeof (int));
    eigen.vectors = V;
    eigen.mat_size = n;

    for (i = 0; i < n; i++) /* after sorting, in [i]=j, j would the be the rank of the i vector */
        eigen.values_index_sorted[i] = i;
    sort_eigen_values(eigen,0,eigen.mat_size);
    /* TODO sort vectors by values */
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
    if (norm == 0) {
        return norm_vec;
    }
    for (i = 0; i < n; i++) {
        norm_vec[i] = (a[i] / norm);
    }
    return norm_vec;
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
    assert(!"goal is not define");
    return 0;



}

void load_string(char** str,char* cpy)
{
    int len;
    len = strlen(cpy);
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
    int i;
    for (i = 0; i < n; i++)
        a[i] =(double) i ;
   swap_array_elements(&a[3],&a[7],sizeof(double));
    swap_array_elements(&a[1],&a[5],sizeof(double));

    print_vector(a,n);

    free(a);
}

int main(int argc, char* argv[])
{
    test_swap();

     int k,n,d;
    n = 1000;
    d = 10;
    char*  goal;
    char*  path_init_centroids;
    char*  file_name;
    double** observations = create_matrix(n,d);
   // assert(argc==4);
      //TODO assert k is an integer
        k = atoi(argv[1]);
        load_string(&goal,argv[2]);
        load_string(&path_init_centroids,"");
        load_string(&file_name,argv[3]);


    assert(k>=0);

    assert_goal(goal);
    printf("k: %d\n goal %s\n init: %s\n file_name:%s\n ",k,goal,path_init_centroids,file_name);
    Tuple2 sizes = load_observations_from_file(observations, file_name);
    n=sizes.i;
    d=sizes.j;

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


    print_mat(observations,n,d);


    /* Free all */

    free(goal);
    free(path_init_centroids);
    free(file_name);
    free_matrix(observations, n);




}

spk_results run_kmean(double** observations,double** centroids,double * init_centroids_index , int k, int n, int d)
{
    double** U = create_matrix(n, k); /* TODO get first k eigenvectors  */
    double** T = create_matrix(n, k);
    renorm_matrix_rows(U, n, T);

    simple_kmean(T, centroids, init_centroids_index, n, k,k);
    /*  TODO find cluster for each T[i], T_At[i]=j
     * place the observations i at centriod j,
     * calculate  new centroies, return them */

}

spk_results activate_flag(char* goal,double** observations,double** centroids,double * init_centroids_index , int k, int n, int d)
{
    spk_results Res;
    Res.k = k;
    Res.eigen.values = NULL;
    Res.eigen.vectors = NULL;
    Res.eigen.mat_size = 0;


    double** W = create_matrix(n, n);
    create_adj_mat(observations,n,d,W);

    if (strcmp(goal,"wam")==0){
        Res.mat = W;
        return Res;
    }

    double** D = create_matrix(n, n);
    create_diagonal_degree_mat_ns(W,n,D);
    if (strcmp(goal,"ddg")==0)
    {
        Res.mat = D;
        free_matrix(W,n);
        return Res;
    }


    double** L = create_matrix(n, n);
    create_L_norm(D,W,L,n);
    if (strcmp(goal,"lnorm")==0)
    {
        Res.mat = L;
        free_matrix(W,n);
        free_matrix(D,n);
        return Res;
    }

    Res.eigen = find_eigen_vectors_and_values(L, n);
    if (strcmp(goal,"jacobi")==0)
    {
        free_matrix(W,n);
        free_matrix(D,n);
        free_matrix(L,n);
        return Res;
    }

    if (k==0)
        k = eigengap_huristic(Res.eigen);





        return Res;
}


