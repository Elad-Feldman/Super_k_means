#include "spkmeans.h"

/* gcc spkmeans.c && gcc  -o spkmeans spkmeans.c && spkmeans  5 spk  dots_10.txt  */

#define Ver 0    \
/* TODO  to zero before submitting */
#define print_verbose(x) if(Ver && printf(x)){}
#define my_assert(cond) assert( (cond) && "An Error Has Occured" )
/****** small function START *******/
void my_assert2(int  cond)
{
    if (!cond){
        printf("hello");
        assert(0);
    }

}

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
     printf("Invalid Input");
     assert(0);
     return 0;



}
/****** small function START *******/

double fix_neg_zero(double num){
    double EPS = 0.00000001;
    if ((num < EPS) && (-1* num > -1*EPS))
        return 0;
    return num;
}


/*************** Vectors  START ******************/
double dot_mult_vector(double *a, double *b, int n) {
    double sum;
    int i;
    assert(a);
    assert(b);
    assert(n > 0);

    sum = 0;
    for (i = 0; i < n; i++)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

double find_vec_norm(double* a, int n) {
    int j;
    double sum;
    double tmp;
    assert(a);
    assert(n > 0);
    sum = 0;
    for ( j = 0; j < n; j++)
    {
        tmp = pow(a[j],2);
        sum += tmp;
    }

    assert(sum>0);
    return sqrt(sum);
}

double find_vec_norm_diff(double* a, double* b, double* c, int n) {
    /* returns the euclidian distance bitween two vectors a, b */
    int i;
    double norm;
    for (i = 0; i < n; i++)
        c[i]= a[i]-b[i];
    norm = find_vec_norm(c,n);
    return norm;
}

double sum_vector(double* a, int n) {
    double sum;
    int i;
    sum = 0;
    for (i = 0; i < n; i++) {
        sum += a[i];
    }
    return sum;
}
void print_vector_int(int* a,int n){
    int i;
    for (i = 0; i < n; i++){
        printf("%d",a[i]);
        if (i<n-1)
            printf(",");

    }

    printf("\n");

}
void print_vector(double* a,int n){
    int i;
    int M = 1;
    for (i = 0; i < n; i++){
        printf("%.4f",a[i]*M);
        if (i<n-1)
            printf(",");

    }

    printf("\n");

}
void swap_int(int *a,int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

void swap_double(double *a, double *b)
{
    double t = *a;

    *a = *b;
    *b = t;

}

void swap_double_pointers(double **a, double **b)
{
    double* t = *a;
    *a = *b;
    *b = t;

}

double* renormlized_vector(double* a, int n) {
    int j;
    double norm;
    double* norm_vec;
    my_assert(a != NULL);
    my_assert(n > 0);
    norm = find_vec_norm(a, n);
    norm_vec = calloc(n, sizeof(double));


    my_assert(norm_vec != NULL);
    if (norm == 0) {
        return norm_vec;
    }
    for (j = 0; j < n; j++) {
        norm_vec[j] = (a[j] / norm);
    }
    return norm_vec;
}
/*************** Vectors  END ******************/


/*************** Matrix  START ******************/
double **create_matrix(int rows, int cols)
{
    int i;
    double** mat = (double  **)calloc(rows , sizeof(double*));
    assert(mat);
    for (i = 0; i < rows; i++) {
        mat[i] = (double  *) calloc(cols,sizeof(double));
        assert(mat[i]);

    }
    return mat;

}

void free_matrix( double  ** A, int rows)
{
    int i;
    for (i = 0; i < rows; i++){
        free(A[i]);
        A[i] = NULL;
    }
    free(A);
    A = NULL;
}

void inplace_transpose_mat(double** mat, int rows, int cols){
    /* works for ros=cols only */
   double** mat_T;
   mat_T = transpose_mat(mat,rows,cols);
   copy_matrix(mat,mat_T,rows,cols);
   free_matrix(mat_T,rows);


}

double** transpose_mat(double** mat, int n, int d)
{
    int i, j;
    double** mat_T;
    assert(mat);
    assert(n > 0);
    assert(d > 0);
    mat_T = create_matrix(d, n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            mat_T[j][i] = mat[i][j];
        }
    }
    return mat_T;
}

void mult_matrix(double** A, double** B, double ** C ,int n, int is_diag) {
    /* is_diag: 0= full, 1= A diag , 2= B diag */

    int i, j;
    double **B_T;
    assert(A);
    assert(B);
    assert(C);
    assert(n > 0);
    assert(is_diag >= 0 && is_diag <= 2);
        B_T = transpose_mat(B, n, n);
    if (is_diag == 0) { /* full matrix multilation  */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                C[i][j] = dot_mult_vector(A[i], B_T[j], n);
        }
        free_matrix(B_T, n);
    }

    if (is_diag == 1) { /* A is diga  */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                C[i][j] = B_T[i][j] * A[i][i];

        }
    }
    if (is_diag == 2) { /* B is diga  */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                C[i][j] = A[i][j] * B[i][i];

        }


    }

}

void copy_matrix(double** A, double** B ,int n,int m){
    /* copy B into A, B override A */
    int i, j;
    assert(A);
    assert(B);
    assert(n > 0);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
            A[i][j] =  B[i][j];

    }
}

void sub_matrix(double** A, double** B, double** C,int n) {
    int i, j;
    for ( i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
    }

}

void print_mat( double  ** mat, int n, int d)
{
    int i;
    print_verbose("==============================\n");
    for (i = 0; i < n; i++) {
        print_vector(mat[i],d);
    }
    print_verbose("==============================\n");
}

double** create_Id_matrix(int n) {
    int i;
    double** mat;
    assert(n > 0);
     mat = create_matrix(n,n);
    for (i = 0; i < n; i++) {
        mat[i][i] = 1; /*  for i=\=j ,calloc  allocated memory block to zero */

    }
    return mat;
}

void re_order_matrix_by_indces(double** A,int* indces, int n)
{
    int i,j;
    for (i=0;i<n;i++){
        j = indces[i];
        swap_double_pointers(&A[i],&A[j]);
    }
}

void renorm_matrix_rows(double** U, int n,int k, double** T)
{
    int i;

    for (i = 0; i < n; i++)
        T[i]= renormlized_vector(U[i], k);
}
/*************** Matrix END ******************/



/*************** kmean  START ******************/
double find_distance(double *dot, double *center, int d){
    double dis;
    int i;
    dis = 0;

    for ( i = 0; i < d; i++)
        dis += (dot[i] - center[i]) * (dot[i] - center[i]);
    return  dis;

}
int get_index_of_closest_cluster(double* dot, double** cluster_list, int d, int k ){
    int j;
    int i;
    double min_dis;
    double tmp_dis;
    j = 0;
    min_dis = find_distance(dot, cluster_list[0], d);

    for (i = 1; i < k; i++)
    {
        tmp_dis = find_distance(dot, cluster_list[i], d);
        if (tmp_dis < min_dis)
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
        printf("remove from  an empty cluster \n ");

    center_temp  = (double *) calloc(d,sizeof(double ));
    for (i = 0; i < d; i++)
        center_temp[i] = (center[i] * (cluster_size));

    for (i = 0; i < d; i++){
        center_temp[i] += (dot[i]*sign);
        center[i] = center_temp[i] / (cluster_size+sign);

    }
    free(center_temp);
}
double** get_init_clusters_list(double** T,int k){
    /* get the first k rows of T */
    double** cluster_list = create_matrix(k,k);
    copy_matrix(cluster_list,T,k,k);
    return cluster_list;
}
int * init_clusters_indexes(int k){
    int * clusters_indexes = calloc(k,sizeof(int));
    int i;
    for(i=0;i<k;i++){
        clusters_indexes[i]=i;
    }
    return clusters_indexes;
}
void simple_kmean (double ** T_mat, double ** T_cluster_list, int* cluster_index_list, int n, int k, int d) {
    int i, j ,max_iter ;
    int is_a_cluster_changed , count_iter;
    int *T_at,   *move_T_to, *T_cluster_size ;

    max_iter = 300;
    T_at = (int*) calloc(n, sizeof (int));
    my_assert2(T_at != NULL);
    move_T_to = (int*) calloc(n, sizeof (int));
    my_assert2(move_T_to != NULL);
    T_cluster_size = (int*) calloc(k, sizeof (int));
    my_assert2(T_cluster_size != NULL);

    for (i = 0; i < n; i++) { /*  dot are  */
        T_at[i] = -1;
        move_T_to[i] = -1;
    }

    for (i = 0; i < k; i++) { /* set initial   locations */
        j =(int) cluster_index_list[i]; /*[ 44,56,73 ] */
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
            if (T_at[j] == -1) { /* don't is not in any center */
                T_at[j] = move_T_to[j];
                update_cluster_center(T_mat[j], T_cluster_list[move_T_to[j]], T_cluster_size[move_T_to[j]], d, 1); /*add dot to center*/
                T_cluster_size[move_T_to[j]]++;
                is_a_cluster_changed = 1;
            } else {
                if (T_at[j] != move_T_to[j]) {
                    update_cluster_center(T_mat[j], T_cluster_list[T_at[j]], T_cluster_size[T_at[j]], d,-1); /*remove dot from center */
                    update_cluster_center(T_mat[j], T_cluster_list[move_T_to[j]], T_cluster_size[move_T_to[j]], d,1); /*add dot to center */
                    T_cluster_size[T_at[j]]--;
                    T_cluster_size[move_T_to[j]]++;
                    T_at[j] = move_T_to[j];
                    is_a_cluster_changed = 1;
                }
            }
        }
    }


    print_vector_int(cluster_index_list,k);
    print_mat(T_cluster_list,k,k);
    free(T_at);
    free(move_T_to);
    free(T_cluster_size);

}
/*************** kmean  END ******************/


/************    Lnorm  START    ********/
void create_adj_mat(double** observations, int n, int d,double** W)
{
    double norm ;
    int i, j;
    my_assert(n > 0);
    my_assert(d > 0);
    my_assert(observations != NULL);
    assert(n > 0);
    assert(d > 0);
    assert(observations);
    double * temp_diff =  (double  *) calloc(n,sizeof(double));

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++) /* matrix is symtric, W[j][i] = W[i][j];   */
        {
            norm = find_vec_norm_diff(observations[i], observations[j],temp_diff , d);
            W[i][j] = exp((-norm) / 2);
            W[j][i] = W[i][j];
        }
    }
    free(temp_diff );
}

void create_diagonal_degree_mat(double** adj_mat,int n,double** D) {
    /* ns stands for negtive squre root this function returns D^(-1/2) */
    int i;
    my_assert(adj_mat != NULL);
    my_assert(n>0);
    for (i = 0; i < n; i++)
        D[i][i] = sum_vector(adj_mat[i],n);


}

void D_sqrt(double** D,int n) {
    /* ns stands for negtive squre root this function returns D^(-1/2) */
    int i;
    my_assert(D != NULL);
    my_assert(n>0);
    for (i = 0; i < n; i++)
    {
        if (D[i][i] != 0)
                D[i][i] = 1/ sqrt(D[i][i]);

    }

}

void create_L_norm(double** D, double** W, int n, double** L ) {
    double** id_mat = create_Id_matrix(n);
    double** DW = create_matrix(n,n);
    double** DWD = create_matrix(n,n);

    mult_matrix(D, W, DW, n,0);
    mult_matrix(DW, D,DWD,n,0);
    sub_matrix(id_mat, DWD,L, n);

    free_matrix(DW,n);
    free_matrix(DWD,n);
    free_matrix(id_mat,n);

}
/************    Lnorm  END    ********/


/************    Jacobi  START    ********/
double abs_d(double x){
    if(x<0){
        return -x;
    }
    return x;
}

void find_ind_max_ele_off_diag(double** A, int n,int* I, int* J)
{
    int i, j;
    double tmp;
    double max_abs;
    my_assert(A != NULL);
    my_assert(n > 0);
    max_abs = -1;

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
    if (A_ij == 0) {
        return 0;
    }
    return (A_jj - A_ii) / (2 * A_ij);
}
double calc_t(double theta) {    return sign(theta)/(abs_d(theta)+sqrt((theta*theta)+1));}
double calc_c(double t) { return 1 / sqrt((t * t) + 1); }
void update_V(int i, int j, int  n, double c, double s,double** V){
    int r;
    double v1;
    double v2;
    for(r = 0; r < n; r++){
        v1 = c * V[r][i] - s * V[r][j];
        v2 = s * V[r][i] + c * V[r][j];
        V[r][i] = v1;
        V[r][j] = v2;
    }
}
double** find_new_A(int i, int j, int n, double c, double s, double** A, double ** A1){
    int r;
    for(r = 0; r < n; r++){
        if(r!=i&&r!=j){
         A1[r][i] = c*A[r][i] - s*A[r][j];
         A1[r][j] = c*A[r][j] + s*A[r][i];
         A1[i][r] = A1[r][i];
         A1[j][r] = A1[r][j];
        }

    }
    A1[i][i] = c*c*A[i][i] + s*s*A[j][j] - 2*s*c*A[i][j];
    A1[j][j] = s*s*A[i][i] + c*c*A[j][j] + 2*s*c*A[i][j];
    A1[i][j] = 0;
    if(i!=j)
    A1[j][i] = 0;
    return A1;
}
double** create_rotation_mat(double** A, double** A1, double** V, int n){
    int i,j;
    double c,t,s, theta;
    my_assert(A != NULL );
    find_ind_max_ele_off_diag(A, n,&i,&j);
    theta = calc_theta(A[j][j],A[i][i],A[i][j]);
    t = calc_t(theta);
    c = calc_c(t);
    s = t * c;
    update_V(i, j, n, c, s, V);
    return find_new_A(i, j, n, c, s, A, A1);
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
    double EPSILON =  1.0*exp(-15);
    double sum_A = sum_square_elements_off_diag( A, n );
    double sum_A1 = sum_square_elements_off_diag( A1, n );
    double diff =  sum_A-sum_A1;
    if( diff <= EPSILON || sum_A1 == 0 )
        return 1;
    else
        return 0;
}
/************    Jacobi  END    ********/
/**** Mergesort start ****/
void merge(double* arr,int* arr_i, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    double* L = (double*) malloc( n1 * sizeof (double));
    double* R = (double*) malloc( n2 * sizeof (double));
    int* L_i = (int*) malloc( n1 * sizeof (int));
    int* R_i = (int*) malloc( n2 * sizeof (int));

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){
        L[i] = arr[l + i];
        L_i[i] = arr_i[l + i];
    }

    for (j = 0; j < n2; j++){
        R[j] = arr[m + 1 + j];
        R_i[j] = arr_i[m + 1 + j];
    }

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            arr_i[k] = L_i[i];
            i++;
        }
        else {
            arr[k] = R[j];
            arr_i[k] = R_i[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1) {
        arr[k] = L[i];
        arr_i[k] = L_i[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2) {
        arr[k] = R[j];
        arr_i[k] = R_i[j];
        j++;
        k++;
    }

    free(L);
    free(R);
    free(L_i);
    free(R_i);
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort(double* arr,int* arr_i, int l, int r)
{
    /* based on https://www.geeksforgeeks.org/merge-sort/ */
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSort(arr,arr_i, l, m);
        mergeSort(arr,arr_i, m + 1, r);

        merge(arr,arr_i, l, m, r);
    }
}

/**** Mergesort end ****/

/**********   Eigen_values START **********/
double* extract_eigen_values_from_mat(double** mat,int n){
    int i;
    double * eigen_values = calloc(n,sizeof(double));
    my_assert(eigen_values != NULL );
    for(i=0;i<n;i++){
        eigen_values[i]=mat[i][i];
    }
    return eigen_values;
}

int check_is_diag(double** A,int n,int k){
    // TODO can be remove, use sum_A in convarge is  zero
    int i, j;
    for(i = 0; i < n; i++){
        for(j = 0; j < k; k++){
            if(i!=j){
                if(A[i][j]!=0){
                    return 0;
                }
            }
        }
    }
    return 1;
}
Eigen find_eigen_vectors_and_values(double** L, int n){
    /* Start with A = L_norm */
    int i;
    int max_iter,convergence;
    int is_diag =0;
    Eigen eigen;
    double **V,  **A, **A_f;
    my_assert(L != NULL);
    my_assert(n>0);
    A = create_matrix(n,n);
    copy_matrix(A,L,n,n);
    print_verbose("start: find eigen vectors");
    V = create_Id_matrix(n);
    A_f = create_matrix(n,n);
    copy_matrix(A_f,A,n,n);
    max_iter = 100;
    convergence = 0;
    i = 0;

    while( !convergence &&  i < max_iter && !is_diag){
        i++;
        A_f = create_rotation_mat(A,A_f,V,n);
        convergence = check_convergence(A,A_f,n);
        copy_matrix(A,A_f,n,n);
          }
    eigen.vectors = V;
    eigen.mat_size = n;
    eigen.values = extract_eigen_values_from_mat(A, n);
    eigen.ranks =  (int*) calloc(n , sizeof (int));
    my_assert(eigen.ranks != NULL);
    free_matrix(A,n);
    for (i = 0; i < n; i++) /* after sorting, in [i]=j, j would the be the rank of the i vector */
        eigen.ranks[i] = i;
    printf("\n eigen vectors 2:\n");
    print_mat(V,n,n);
    printf("\n ----------");
     return eigen;
}

int  eigengap_huristic(Eigen eigen){
    int k,m,i;
    double  delta_i, max;
    m =  (int ) eigen.mat_size / 2; /* it said that k<n/2 */
    max = - 1 ;
    k = 0;
    for (i = 0 ; i < m ; i++){
        delta_i = abs_d (eigen.values[i] - eigen.values[i+1]) ;
        if (delta_i > max)
        {
            max = delta_i ;
            k=i+1;
        }

    }
    my_assert(k>0);
    return k;
}

void eigen_to_matrix(Eigen eigen, double** E,int n)
{

    int i,j;
    for(j=0; j<n;j++){
        E[0][j] = eigen.values[j];
    }
    for(i=1; i<=n;i++)
    {
        for(j=0; j<n;j++)
            E[i][j] = eigen.vectors[i-1][j];
    }

}

void free_eigen(Eigen eigen){
    free_matrix(eigen.vectors,eigen.mat_size);
    free(eigen.ranks);
    free(eigen.values);
    eigen.ranks = NULL;
    eigen.values = NULL;
}
/**********   Eigen_values END **********/

void test_stable_sort(){
    int perm = 6; /* 0 7 8 */
    int n = 8;
    int i,j;
    double d0[] = {5,9,3,4,5,2,5,1};
    double** double_mat  = (double**) malloc(perm * sizeof(double*));
    int** ind_mat  = (int**) malloc(perm * sizeof(int*));
    for (i=0;i<perm;i++)
    {
        ind_mat[i] = (int*)  malloc(n * sizeof(int));
        double_mat[i] = (double*)  malloc(n * sizeof(double));

        for (j=0;j<n;j++)
        {
            double_mat[i][j] =d0[j]; /* create 6 version of d0 */
            ind_mat[i][j] = 6;

        }

    }

    ind_mat[0][0] = 0;ind_mat[0][4] = 7;ind_mat[0][6] = 8;
    ind_mat[1][0] = 0;ind_mat[1][4] = 8;ind_mat[1][6] = 7;
    ind_mat[2][0] = 7;ind_mat[2][4] = 0;ind_mat[2][6] = 8;
    ind_mat[3][0] = 7;ind_mat[3][4] = 8;ind_mat[3][6] = 0;
    ind_mat[4][0] = 8;ind_mat[4][4] = 0;ind_mat[4][6] = 7;
    ind_mat[5][0] = 8;ind_mat[5][4] = 7;ind_mat[5][6] = 0;

    for (i=0;i<perm;i++)
    {
        printf("===PERM ID %d===\nBEFORE:",i+1);
        printf("perm: %d,%d,%d\n",ind_mat[i][0],ind_mat[i][4],ind_mat[i][6]);
        mergeSort(double_mat[i],ind_mat[i],0,n-1);
      /*   for (j=0;j<n;j++)
             printf("%.0f,",double_mat[i][j]); */
        printf("AFTER:");
         printf("perm: %d,%d,%d\n",ind_mat[i][4],ind_mat[i][5],ind_mat[i][6]);
        for (j=0;j<n;j++) {
             printf("%.0f,",double_mat[i][j]);
         }
        printf("\n ");
        free(double_mat[i]);
        free(ind_mat[i]);

        }
    free(double_mat);
    free(ind_mat);
    }


/***************  SPK START ******************/

void start_wam(double** observations , int n, int d, double*** W)
{
    *W = create_matrix(n, n);
    create_adj_mat(observations,n,d,*W);

}
void start_ddg(double** observations , int n, int d, double*** D)
{
    double** W = create_matrix(n, n);
    *D = create_matrix(n, n);
    create_adj_mat(observations,n,d,W);
    create_diagonal_degree_mat(W,n,*D);
    free_matrix(W,n);
}

void start_lnorm(double** observations , int n, int d, double*** L)
{
    double** W = create_matrix(n, n);
    double** D = create_matrix(n, n);
    *L = create_matrix(n, n);
    create_adj_mat(observations,n,d,W);
    create_diagonal_degree_mat(W,n,D);
    D_sqrt(D,n);
    create_L_norm(D,W,n,*L);
    free_matrix(W,n);
    free_matrix(D,n);
}

void start_jacobi(double** observations , int n, double*** E)
{
    Eigen eigen = find_eigen_vectors_and_values(observations, n);
    inplace_transpose_mat(eigen.vectors,n,n); /*  now each row is a vector */
    printf("this is second jacobi:\n");
    print_mat(eigen.vectors,n,n);
    printf("\nsecond values:\n");
    print_vector(eigen.values,n);
   *E  = create_matrix(n+1, n);
    eigen_to_matrix(eigen,*E,n);
}


spk_results activate_flag(char* goal,double** observations , int k, int n, int d)
{
    /* run all the flags, that are not spk */

    int i,j;
    spk_results res;
    Eigen eigen;
    double **W ,**D ,**L, **E, **U, **T;
    res.T_size = n;
    res.k = n;



    if (is_goal("wam")){
        start_wam(observations,n,d,&W);
        res.T = W;
        print_mat(res.T,res.T_size,n);
        return res;
    }
    if (is_goal("ddg"))
    {
        start_ddg(observations,n,d,&D);
        res.T = D;
        print_mat(res.T,res.T_size,n);
        return res;
    }
    if (is_goal("lnorm"))
    {
        start_lnorm(observations,n,d,&L);
        res.T = L;
        print_mat(res.T,res.T_size,n);
        return res;
    }
    if (is_goal("jacobi"))
    {
       // start_jacobi(observations,n,&E);
        start_jacobi(observations,n,&E);
        printf("\n");
        res.T_size = n + 1;
        res.T = E;
        //print_mat(res.T,res.T_size,n);
        return res;
    }
    /* else goal = full  spk  */


    W = create_matrix(n, n);
    create_adj_mat(observations,n,d,W);

    D = create_matrix(n, n);
    create_diagonal_degree_mat(W,n,D);
    D_sqrt(D,n);

    L = create_matrix(n, n);
    create_L_norm(D,W,n,L);
    eigen = find_eigen_vectors_and_values(L, n);

    mergeSort(eigen.values,eigen.ranks,0,n-1);
    inplace_transpose_mat(eigen.vectors,n,n); /*  now each row is a vector */
    re_order_matrix_by_indces(eigen.vectors, eigen.ranks, n);
     if (k==0) //TODO what about k<0 ?
        k = eigengap_huristic(eigen);


    U  = create_matrix(n, k);
    for (i=0; i<n; i++)
        {
            for (j=0; j<k; j++)
                U[i][j] = eigen.vectors[j][i]; /* each vector is a column */
        }
    T = create_matrix(n,k);
    renorm_matrix_rows(U, n,k, T);
  /*  printf("C: n=%d,  k=%d\n",n,k);
    print_mat(T,n,k); */
    free_matrix(W,n);
    free_matrix(D,n);
    free_matrix(L,n);
    free_matrix(U,n);;
    free_eigen(eigen);
    res.T = T;
    res.k = k;
    return res;
}
/***************  SPK END ******************/

/******** C Interface ******/
void load_string(char** str,char* cpy)
{
    int len;
    len = (int) strlen (cpy);
    *str = (char  *) malloc(len * sizeof(char));
    my_assert(*str != NULL);
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
    char* row;
    Tuple2 sizes;

    fp = fopen(file_name,"r");
    row = (char* )calloc(1000, sizeof(char));
    my_assert(row != NULL);
    i =0;
    while(fscanf(fp,"%s",row)==1)
    { /* load data */
        d = string_to_doubles(row, observations[i]);
        i++;
    }
    n = i;
    /* change the size of observations to match the file */
    observations = (double **) realloc(observations,n* sizeof(double *));
    my_assert(observations != NULL);
    for (i = 0; i < n; i++)
    {
        observations[i] = (double *) realloc(observations[i],d* sizeof(double));
        my_assert(observations[i] != NULL);
    }
    fclose(fp);
    free(row);
    sizes.i = n;
    sizes.j = d;
    return sizes;

}



/***** TEST QSORT ******/


int main(int argc, char* argv[])
{

   // test_stable_sort();

    int k,n,d;
    char*  goal;
    char*  file_name;
    double** T_clusters_list;
    int * T_clusters_indexes;
    double** observations;
    Tuple2 sizes;
    spk_results Res;
    n = 1000;  d = 10;
    observations = create_matrix(n,d);
    my_assert(argc==4);
    // TODO my_assert k is an integer
    k = atoi(argv[1]);
    load_string(&goal,argv[2]);
    load_string(&file_name,argv[3]);


    my_assert(k>=0); // TODO change massage

    assert_goal(goal);
   // printf("=========\nk:%d\ngoal %s\nfile_name: %s\n==========\n",k,goal,file_name);
    sizes = load_observations_from_file(observations, file_name);
    n=sizes.i;    d=sizes.j;

    Res = activate_flag( goal, observations , k,  n, d);
    if(is_goal("spk")){
        k = Res.k;
      //  printf("found k: %d \n",Res.k);
        T_clusters_list = get_init_clusters_list(Res.T,k);
        T_clusters_indexes = init_clusters_indexes(k);
        simple_kmean(Res.T, T_clusters_list, T_clusters_indexes,n,k,d);

        free_matrix(T_clusters_list,k);
        free(T_clusters_indexes);
    }



    //Free all
    free(goal);
    free(file_name);
    free_matrix(observations, n);
    free_matrix(Res.T,Res.T_size );
    printf("\n  C done !");
    return 0;
}






