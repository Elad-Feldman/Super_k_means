#include "spkmeans.h"

/* gcc spkmeans.c && gcc  -o spkmeans spkmeans.c && spkmeans  5 spk  dots_10.txt*/

/****** small function START *******/
void my_assert(int  cond)
{
    if (!cond){
        printf("An Error Has Occured");
        assert(0);
    }

}
void print_verbose(char* string){
    if (verbose)
        printf("%s",string);
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
    assert(!&&"Invalid Input"); //TODO invalid input
    return 0;



}
/****** small function START *******/



/*************** Vectors  START ******************/
double dot_mult_vector(double *a, double *b, int n) {
    assert(a);
    assert(b);
    assert(n > 0);
    double sum;
    int i;
    sum = 0;
    for (i = 0; i < n; i++)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

double find_vec_norm(double* a, int n) {
    int i;
    double sum;
    assert(a);
    assert(n > 0);
    sum = 0;
    for ( i = 0; i < n; i++)
    {
        sum += (a[i] * a[i]);
    }
    return sqrt(sum);
}

double find_vec_norm_diff(double* a, double* b, int n) {
    /* returns the euclidian distance bitween two vectors a, b */
    int i;
    double norm;
    double * c =  (double  *) calloc(n,sizeof(double));
    for (i = 0; i < n; i++)
        c[i]= a[i]-b[i];
    norm = find_vec_norm(c,n);
    free(c);
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
    my_assert(a != NULL);
    my_assert(n > 0);
    int i;
    double norm = find_vec_norm(a, n);
    double* norm_vec = calloc(n, sizeof(double));
    my_assert(norm_vec != NULL);
    if (norm == 0) {
        return norm_vec;
    }
    for (i = 0; i < n; i++) {
        norm_vec[i] = (a[i] / norm);
    }
    return norm_vec;
}
/*************** Vectors  END ******************/


/*************** Matrix  START ******************/
double **create_matrix(int n, int d)
{
    double** mat = (double  **)calloc(n , sizeof(double*));
    assert(mat);
    int i;
    for (i = 0; i < n; i++) {
        mat[i] = (double  *) calloc(d,sizeof(double));
        assert(mat[i]);

    }
    return mat;

}

void free_matrix( double  ** A, int n)
{
    int i;
    for (i = 0; i < n; i++){
        free(A[i]);
        A[i] = NULL;
    }
    free(A);
    A = NULL;
}

double** transpose_mat(double** mat, int n, int D)
{
    assert(mat);
    assert(n > 0);
    assert(D > 0);
    double** mat_T = create_matrix(D, n);
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < D; j++)
        {
            mat_T[j][i] = mat[i][j];
        }
    }
    return mat_T;
}

void mult_matrix(double** A, double** B, double ** C ,int n) {
    assert(A);
    assert(B);
    assert(C);
    assert(n > 0);

    double** B_T = transpose_mat(B, n, n);
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            C[i][j] = dot_mult_vector(A[i], B_T[j], n);
        }
    }
    free_matrix(B_T,n);
}

void copy_matrix(double** A, double** B ,int n){
    assert(A);
    assert(B);
    assert(n > 0);


    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] =  B[i][j];
        }
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
    assert(n > 0);
    double** mat = create_matrix(n,n);
    int i;
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

void renorm_matrix_rows(double** U, int n, double** T)
{
    int i;
    for (i = 0; i < n; i++)
        T[i]= renormlized_vector(U[i], n);
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
double** get_init_clusters_list(double** T,int k){
    /* get the first k rows of T */
    double** cluster_list = create_matrix(k,k);
    copy_matrix(cluster_list,T,k);
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
/*************** kmean  END ******************/


/************    Lnorm  START    ********/
void create_adj_mat(double** observations, int n, int d,double** W)
{
    my_assert(n > 0);
    my_assert(d > 0);
    my_assert(observations != NULL);
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
    my_assert(adj_mat != NULL);
    my_assert(n>0);
    int i;
    double sum_of_weights;//sum of weights per row
    sum_of_weights = 0;
    for (i = 0; i < n; i++)
    {
        sum_of_weights = sum_vector(adj_mat[i],n);
        if (sum_of_weights==0)
            print_vector(adj_mat[i],n);
        // my_assert(sum_of_weights > 0);
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
    my_assert(A != NULL);
    my_assert(n > 0);
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
double calc_t(double theta) {    return sign(theta)/(abs_d(theta)+sqrt((theta*theta)+1));}
double calc_c(double t) { return 1 / sqrt((t * t) + 1); }
double** create_rotation_mat(double** A,int n){
    my_assert(A != NULL );
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
/************    Jacobi  END    ********/







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
    my_assert(L != NULL);
    my_assert(n>0);
    double ** A = create_matrix(n,n);
    copy_matrix(A,L,n);
    print_verbose("start: find eigen vectors");

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
    print_verbose("\nfound vectors!\n");

    eigen.vectors = V;
    eigen.mat_size = n;
    eigen.values = extract_eigen_values_from_mat(A, n);
    eigen.ranks =  (int*) calloc(n , sizeof (int));
    my_assert(eigen.ranks != NULL);
    free_matrix(A,n);
    for (i = 0; i < n; i++) /* after sorting, in [i]=j, j would the be the rank of the i vector */
        eigen.ranks[i] = i;
    //TODO sort only for spk, not for jacobi
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
    my_assert(k>0);
    return k;
}

void free_eigen(Eigen eigen){
    free_matrix(eigen.vectors,eigen.mat_size);
    free(eigen.ranks);
    free(eigen.values);
    eigen.ranks = NULL;
    eigen.values = NULL;
}
/**********   Eigen_values END **********/



/***************  SPK START ******************/
spk_results activate_flag(char* goal,double** observations , int k, int n, int d)
{
    // TODO make it shorter - split
    /* run all the flags, that are not spk */
    spk_results Res;
    Res.k = k;
    Res.eigen.values = NULL;
    Res.eigen.vectors = NULL;
    Res.eigen.mat_size = 0;
    double** W = create_matrix(n, n);
    create_adj_mat(observations,n,d,W);
    if (is_goal("wam")){
        print_verbose("wam:\n");
        print_mat(W,n,n);
        free_matrix(W,n);
        return Res;
    }

    double** D = create_matrix(n, n);
    create_diagonal_degree_mat_ns(W,n,D);
    if (is_goal("ddg"))
    {
        print_verbose("ddg:\n");
        print_mat(D,n,n);
        free_matrix(W,n);
        free_matrix(D,n);//
        return Res;
    }


    double** L = create_matrix(n, n);
    create_L_norm(D,W,n,L);
    if (is_goal("lnorm"))
    {
        print_verbose("lnorm:\n");
        print_mat(L,n,n);
        free_matrix(L,n);
        free_matrix(W,n);
        free_matrix(D,n);
        return Res;
    }
    Res.eigen = find_eigen_vectors_and_values(L, n);
    if (is_goal("jacobi")){
        print_verbose("jacobi:\n");
        print_verbose("eigen vectors:\n");
        print_mat(Res.eigen.vectors,n,n);
        print_verbose("eigen values:\n");
        print_vector(Res.eigen.values,n); //TODO transpose
        free_matrix(W,n);
        free_matrix(D,n);
        free_matrix(L,n);
        print_verbose("finish jacobi");
        return Res;
    }

    if (k==0)
        k = eigengap_huristic(Res.eigen);
    int i;
    double** U  = (double**)  malloc(n* sizeof (double *)) ;
    my_assert(U != NULL);
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
    free_eigen(Res.eigen);


    return Res;
}
/***************  SPK END ******************/

/******** C Interface ******/
void load_string(char** str,char* cpy)
{
    int len;
    len = (int)strlen(cpy);
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
    fp = fopen(file_name,"r");
    char* row = (char* )calloc(1000, sizeof(char));
    my_assert(row != NULL);
    i =0;
    while(fscanf(fp,"%s",row)==1){ // load data
        d = string_to_doubles(row, observations[i]);
        i++;
    }
    n=i;

    // change the size of observations to match the file
    observations = (double **) realloc(observations,n* sizeof(double *));
    my_assert(observations != NULL);
    for (i = 0; i < n; i++)
    {
        observations[i] = (double *) realloc(observations[i],d* sizeof(double));
        my_assert(observations[i] != NULL);
    }
    fclose(fp);
    free(row);
    Tuple2 sizes;
    sizes.i = n;
    sizes.j = d;
    return sizes;

}

int main(int argc, char* argv[])
{
    int k,n,d;
    n = 1000;  d = 10;
    char*  goal;
    char*  file_name;
    double** observations = create_matrix(n,d);
    my_assert(argc==4);
    //TODO my_assert k is an integer
    k = atoi(argv[1]);
    load_string(&goal,argv[2]);
    load_string(&file_name,argv[3]);


    my_assert(k>=0);// TODO change massage

    assert_goal(goal);
    printf("k: %d\ngoal %s\nfile_name: %s\n",k,goal,file_name);
    Tuple2 sizes = load_observations_from_file(observations, file_name);
    n=sizes.i;    d=sizes.j;

    spk_results Res;
    Res = activate_flag( goal, observations , k,  n, d);
    print_verbose("\nfinish activate_flag\n");
    if(!is_goal("spk")){
        return 0;
    }
    k=Res.k;
    printf("found k: %d",Res.k);
    printf("create full spk here");
    double** T_clusters_list = get_init_clusters_list(Res.mat,k);
    int * T_clusters_indexes = init_clusters_indexes(k);
    simple_kmean(Res.mat, T_clusters_list, T_clusters_indexes,observations,n,k,d);

    //Free all
    free(goal);
    free(file_name);
    free_matrix(T_clusters_list,k);
    free(T_clusters_indexes);
    free_matrix(observations, n);
    printf("\n  C done !");

}






