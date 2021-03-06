#include "spkmeans.h"
#define MAXSIZE 50
/* gcc spkmeans.c && gcc  -o spkmeans spkmeans.c && spkmeans  5  jacobi   tests/test_data/jacobi_tests/test7.csv  */


/****** small function START *******/

double fix_neg_zero(double num){ /* round  numbers and fix the proble, of -0.0000 to 0.0000 */
    if (-0.00005 < num && num < 0)
        num = 0;
    return num;
}
void super_assert(int cond){ /* if  conditions not apply  print the required message and exit */
    if (cond) {
        return;
    }
    printf("An Error Has Occured");
    exit(0); /* https://stackoverflow.com/questions/15495004/memory-leaks-occur-if-exitexitcode-in-c  */


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

/****** small function  END *******/


/*************** Vectors  START ******************/


void print_vector_int(int* a,int n){ /* print an array of integers in a row */
    int i;
    for (i = 0; i < n; i++){
        printf("%d",a[i]);
        if (i<n-1)
            printf(",");

    }
    printf("\n");

}
void print_vector(double* a,int n)
{ /* print an array of doubles, with 4 digits, in a row */
    int i;
    for (i = 0; i < n; i++){
       printf("%.4f", fix_neg_zero(a[i]) ); /* TOOD fix negative zero */
       if (i<n-1)
           printf(",");
         }
    }


double dot_mult_vector(double *a, double *b, int n) { /* scalar multiplied  of 1d  vectors */
    double sum;
    int i;
    assert_not_null(a);
    assert_not_null(b);
    assert_positive(n);

    sum = 0;
    for (i = 0; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}
double find_vec_norm(double* a, int n) { /* find the norm of  a vector  */
    int j;
    double sum,tmp;
    assert_not_null(a);
    assert_positive(n);
    sum = 0;
    for ( j = 0; j < n; j++)
    {
        tmp = pow(a[j],2); /* fix rounding error */
        sum +=  tmp;
    }


    return sqrt(sum);
}
double find_vec_norm_diff(double* a, double* b, int n) {
    /* returns the euclidean distance between two vectors a, b */
    int i;
    double norm;
    double* c ;
    c =  (double*) calloc(n, sizeof(double) );
    assert_not_null(c);
    for (i = 0; i < n; i++) {
        c[i]= a[i]-b[i];
    }
    norm = find_vec_norm(c,n);
    free(c);
    return norm;
}

double sum_vector(double* a, int n) { /* sum the elements of a vector */
    double sum;
    int i;
    sum = 0;
    for (i = 0; i < n; i++)
        sum += a[i];
    return sum;
}


double* renormlized_vector(double* a, int n) { /* find the norm of the vector and normalize the vectors */
    int j;
    double norm;
    double* norm_vec;
    assert_not_null(a);
    assert_positive(n);
    norm = find_vec_norm(a, n);

    norm_vec = calloc(n, sizeof(double));
    assert_not_null(norm_vec);

    if (norm == 0) /* SHOULD NOT HAPPEN AS SAID IN https://moodle.tau.ac.il/mod/forum/discuss.php?d=163019 _ THIS IS a ZERO VECTOR*/
    {
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
{ /* allcote memory of a row X cols matrix and return it */
    int i;
    double** mat;
    mat = (double  **) calloc(rows , sizeof( double*) );
    assert_not_null(mat);
    for (i = 0; i < rows; i++) {
        mat[i] = (double  *) calloc( cols, sizeof( double) );
        assert_not_null(mat[i]);

    }
    return mat;

}

void free_matrix( double  ** A, int rows)
{ /* free matrix A, free each row and then A */
    int i;
    for (i = 0; i < rows; i++){
        free(A[i]);
    }
    free(A);

}

void inplace_transpose_mat(double** mat, int rows, int cols){
    /* works for ros=cols only */
   double** mat_T;
   mat_T = transpose_mat(mat,rows,cols);
   copy_matrix(mat,mat_T,rows,cols);
   free_matrix(mat_T,rows);


}

double** transpose_mat(double** mat, int rows, int cols) /* transpose the matrix and return a copy of the T_mat */
{
    int i, j;
    double** mat_T;
    assert_not_null(mat);
    assert_positive(rows );
    assert_positive(cols );
    mat_T = create_matrix(rows, cols);
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            mat_T[j][i] = mat[i][j];
        }
    }
    return mat_T;
}

void mult_matrix(double** A, double** B, double ** C ,int n) {
    /* multiple A and B, save result into C */

    int i, j;
    double **B_T;
    assert_not_null(A);
    assert_not_null(B);
    assert_not_null(C);
    assert_positive(n);
    B_T = transpose_mat(B, n, n); /* would be easier to multplie this way */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                C[i][j] = dot_mult_vector(A[i], B_T[j], n);
        }
        free_matrix(B_T, n);
}


void copy_matrix(double** A, double** B ,int n,int m){
    /* copy B into A, B override A, assume dim(A)=dim(B) */
    int i, j;
    assert_not_null(A);
    assert_not_null(B);
    assert_positive(n );
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
            A[i][j] =  B[i][j];

    }
}

void sub_matrix(double** A, double** B, double** C,int n) {
    /*  C = A - B */
    int i, j;
    for ( i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
    }

}

void print_mat( double  ** mat, int n, int d)
{ /* print  matrix with n rows and d cols */
    int i;
    for (i = 0; i < n; i++) {
        print_vector(mat[i],d);
        if (i<n-1)
            printf("\n");
    }

}

double** create_Id_matrix(int n) { /*  create a sqaure matrix n X n , with 1 on diagonal, else zeros */
    int i;
    double** mat;
    assert_positive(n);
    mat = create_matrix(n,n);
    for (i = 0; i < n; i++) {
        mat[i][i] = 1; /*  for i != j ,calloc  allocated memory block to zero */

    }
    return mat;
}

void re_order_matrix_by_indces(double** A,int* indces, int n)
{ /* get a matrix A and switch the rows of it by indeces order  */

    double** hold_row = create_matrix(n,n);
    int i,j,P;
    for (i = 0; i<n; i++) {
        P = indces[i];
        for (j = 0; j < n; j++) { /*coping A into hold row */
            hold_row[i][j] = A[P][j];
        }
    }
    for (i=0;i<n;i++) {
        for (j = 0; j < n; j++) { /* override A orginial values */
            A[i][j] =  hold_row[i][j];
        }
    }
    free_matrix(hold_row,n);
}



void renorm_matrix_rows(double** U, int n,int k, double** T)
{ /* Create Vector and allocate memory for each row */
    int i;
    for (i = 0; i < n; i++)
        T[i] = renormlized_vector(U[i], k);
}
/*************** Matrix END ******************/



/*************** kmean  START ******************/
double find_distance(double *dot, double *center, int d){
    /* find the distande of dot from the center  */
    double dis,dif;
    int i;
    dis = 0;
        for ( i = 0; i < d; i++){
            dif = dot[i] - center[i];
            dis += pow(dif,2);
        }
    return  dis;

}
int get_index_of_closest_cluster(double* dot, double** cluster_list, int d, int k ){
    /* find the index of the cluster with the minimal distance  */
    int j;
    int i;
    double min_dis;
    double tmp_dis;
    j = 0; /* start from the lowest index, increase only for shorter distance */
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
    /* remove or add the dot the cluster and update the center value  */
    double* center_temp;
    int i;
    if (cluster_size + sign == 0){
               return;
               /* remove dot from cluster, but don't change cluster position  SHOULD NOT HAPPEN */
    }

    if (cluster_size==0){
       /* SHOULD NOT HAPPAN */
        return;
    }

    center_temp  = (double *) calloc(d, sizeof(double ));
    assert_not_null(center_temp);
    for (i = 0; i < d; i++)
        center_temp[i] = (center[i] * (cluster_size));

    for (i = 0; i < d; i++){
        center_temp[i] += ( dot[i]* sign );
        center[i] = center_temp[i] / (cluster_size+sign);

    }
    free(center_temp);
}
double** get_init_clusters_list(double** T,int k){
    /* get the first k rows of T */
    double** cluster_list ;
    cluster_list = create_matrix(k,k);
    copy_matrix(cluster_list,T,k,k);
    return cluster_list;
}
int * init_clusters_indexes(int k){
    /*  get an array size k where arr[i] = i*/
    int i;
    int * clusters_indexes ;
    clusters_indexes =(int *) calloc(k, sizeof(int) );
    assert_not_null(clusters_indexes);
    for(i=0;i<k;i++){
        clusters_indexes[i] = i;
    }
    return clusters_indexes;
}
void simple_kmean (double ** T_mat, double ** T_cluster_list, int* cluster_index_list, int n, int k, int d,int is_Py_call) {
    /*  preform kmean on T, using the cluster list and indexes*/
    int i, j ,max_iter ;
    int is_a_cluster_changed , count_iter;
    int *dot_at,   *move_dot_to, *T_cluster_size ;
    max_iter = 100;

    dot_at = (int*) calloc(n, sizeof(int) );
    assert_not_null(dot_at);

    move_dot_to = (int*) calloc(n, sizeof(int) );
    assert_not_null(move_dot_to);

    T_cluster_size = (int*) calloc(k, sizeof(int) );
    assert_not_null(T_cluster_size);

    for (i = 0; i < n; i++) { /*  set clusters */
        dot_at[i] = -1;
        move_dot_to[i] = -1;
    }

    for (i = 0; i < k; i++) { /* set initial   locations */
        j =(int) cluster_index_list[i]; /* this is the j dot*/
        dot_at[j] = i;
        T_cluster_size[i] = 1;

    }

    is_a_cluster_changed = 1;
    count_iter = 0;

    while (count_iter < max_iter && is_a_cluster_changed) {
        is_a_cluster_changed = 0;
        count_iter++;

        for (i = 0; i < n; i++){
            j = get_index_of_closest_cluster(T_mat[i], T_cluster_list, d, k);
            move_dot_to[i] = j;

        } /*find nearest cluster for each dot */
        for (j = 0; j < n; j++) {/* update clusters*/
            if (dot_at[j] == -1) { /* dot is outside of all cluster */
                dot_at[j] = move_dot_to[j]; /* place dot at new place */
                update_cluster_center(T_mat[j], T_cluster_list[ move_dot_to[j] ], T_cluster_size[ move_dot_to[j] ], d, 1); /*add dot to center*/
                T_cluster_size[ move_dot_to[j] ]++;
                is_a_cluster_changed = 1;
            }
            else
                {
                if (dot_at[j] != move_dot_to[j])
                {
                    update_cluster_center(T_mat[j], T_cluster_list[ dot_at[j] ], T_cluster_size[ dot_at[j] ], d,-1); /*remove dot from center */
                    update_cluster_center(T_mat[j], T_cluster_list[ move_dot_to[j] ], T_cluster_size[ move_dot_to[j] ], d,1); /*add dot to center */
                    T_cluster_size[ dot_at[j] ]--;
                    T_cluster_size[ move_dot_to[j] ]++;
                    dot_at[j] = move_dot_to[j];
                    is_a_cluster_changed = 1;
                }
            }
        }
    }

    if ( is_Py_call) /* print indexes only for kmeans++ */
        print_vector_int(cluster_index_list,k);
    print_mat(T_cluster_list,k,k);

    free(dot_at);
    free(move_dot_to);
    free(T_cluster_size);



}
/*************** kmean  END ******************/


/************    Lnorm  START    ********/
void create_adj_mat(double** observations, int n, int d,double** W)
{   /*  create the adj matrix for observations */
    double norm ;
    int i, j;
    assert_not_null(observations);
    assert_positive(n);
    assert_positive(d);


    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++) /* matrix is symtric, W[j][i] = W[i][j];   */
        {
            norm = find_vec_norm_diff(observations[i], observations[j] , d);
            W[i][j] = exp( (-norm) / 2);
            W[j][i] = W[i][j];
        }
    }


}

void create_diagonal_degree_mat(double** adj_mat,int n,double** D) { /*create_diagonal_degree_mat from the adj matrix */
        int i;
    assert_not_null(adj_mat);
    assert_positive(n);
    for (i = 0; i < n; i++)
        D[i][i] = sum_vector(adj_mat[i],n);


}

void D_sqrt(double** D,int n) {
    /* get  diagonal_degree_mat and for each  diag element change it to 1/sqrt(element) */
    int i;
    assert_not_null(D);
    assert_positive(n);
    for (i = 0; i < n; i++)
    {
        if (D[i][i] != 0)
                D[i][i] = 1/ sqrt(D[i][i]);

    }

}

void create_L_norm(double** D, double** W, int n, double** L ) { /*  create L norm from diagonal_degree_mat and  adj matrix*/
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
double abs_d(double x){ /* get the abs( x ) */
    if(x<0)
        return -1.0 * x;
    return x;
}

void find_ind_max_ele_off_diag(double** A, int n,int* I, int* J)
{ /* find the indexes of the  biggest  absolute value  element of A that is off diag */
    int i, j;
    double tmp;
    double max_abs;
    assert_not_null(A);
    assert_positive(n);
    max_abs = -1.0;

    for (i = 0; i < n; i++)
    {
        for (j = i; j < n; j++) /* off diagonalm */
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

double sign(double x) { /* return sign of x  */
    if (x < 0)
        return -1.0;

    return 1.0;
}
double calc_theta(double a_jj ,double a_ii,double a_ij) { /*calc_theta */
    if (a_ij == 0) {
        return 0;
    }

    return (a_jj - a_ii) / (2.0 * a_ij);
}
double calc_t(double theta)
{ /*calc_t */
    double  mone,mechne;
    mone = sign(theta) ;
    mechne = abs_d( theta ) + sqrt( pow(theta,2) + 1.0 );
    return  mone / mechne ;
}
double calc_c(double t) { return 1.0 / sqrt( (pow(t,2)) + 1.0); } /* calc_c */


void update_V(int i, int j, int  n, double c, double s,double** V){
    /*update V matrix, taking into account the fact the P is almost diagonal to calculate only the relevant elements  */
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

void update_A_f(int i, int j, int n, double c, double s, double** A, double ** A1){ /* find A' using the instructions to A' definition */
    int r;
    for(r = 0; r < n; r++){
        if( r != i && r != j ){
            A1[r][i] = c * A[r][i] - s * A[r][j];
            A1[r][j] = c * A[r][j] + s * A[r][i];
            A1[i][r] = A1[r][i]; /* A1 is symatric */
            A1[j][r] = A1[r][j];
        }

    }
    A1[i][i] = pow(c,2) * A[i][i] +pow(s,2) * A[j][j] - 2.0 * s * c * A[i][j];
    A1[j][j] = pow(s,2) * A[i][i] + pow(c,2) * A[j][j] + 2.0 * s * c * A [i][j];
    A1[i][j] = 0;
    A1[j][i] = 0;


}

double sum_square_elements_off_diag(double ** A,int n){ /* sum all the elements off diagonal of A */
    int i,j;
    double sum;
    sum = 0;
    for (i = 0; i < n ;i++)
    {
        for (j = i; j < n ;j++) /* take only elements where i<j, A is stmetric so we take 2* A[i][j]^2 */
        {
            if(i!=j)
            {
                sum += 2 * pow(A[i][j],2);
            }

        }
    }
    return sum;
}
int check_convergence(double** A,double** A1,int n)
    { /* check the diff between A and A' and if A' is diagonal  */
     double eps, sum_A, sum_A1, diff;
     eps =   1.0 * pow(10,-15);
     sum_A = sum_square_elements_off_diag( A, n );
     sum_A1 = sum_square_elements_off_diag( A1, n );
     diff =  sum_A - sum_A1;
    if( diff <= eps  || sum_A1 == 0) /* sum_A1 == 0 A1 is diagonal */
        return 1;
    else
        return 0;
}

Eigen find_eigen_vectors_and_values(double** L, int n){ /* this is the jacobi algorithem */
    /* Start with A = L_norm */

    int iter_count, i, j;
    int max_iter,convergence;
    double **V, **A, **A_f;
    Eigen eigen;
    double c,t,s, theta;

    assert_not_null(L);
    assert_positive(n);
    A = create_matrix(n,n);
    A_f = create_matrix(n,n);
    V = create_Id_matrix(n);
    copy_matrix(A,L,n,n);
    copy_matrix(A_f,L,n,n);

    max_iter = 100;
    convergence = 0;
    iter_count = 0;

    while( !convergence &&   iter_count < max_iter  )
    {
        iter_count++;

        find_ind_max_ele_off_diag(A, n, &i, &j);
        theta = calc_theta( A[j][j], A[i][i], A[i][j] );
        t = calc_t(theta);
        c = calc_c(t);
        s = t * c;

        update_V(i, j, n, c, s, V); /* V = P(i) * P_(i+1) */
        update_A_f(i, j, n, c, s, A, A_f);
        convergence = check_convergence(A,A_f,n);
        copy_matrix(A,A_f,n,n);
    }

    eigen.vectors = V;
    eigen.mat_size = n;
    eigen.values = extract_eigen_values_from_mat(A, n);

    eigen.ranks =  (int*) calloc(n , sizeof (int));
    assert_not_null(eigen.ranks);



    for (i = 0; i < n; i++)
        eigen.ranks[i] = i;

    free_matrix(A,n);
    free_matrix(A_f,n);
    return eigen;
}

/************    Jacobi  END    ********/



/**** Mergesort start ****/
void merge(double* arr,int* arr_i, int l, int m, int r)
{ /*  merge two slices of the arr */
    double *L, *R;
    int *L_i, *R_i;
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    L = (double*) malloc( n1 * sizeof (double));
    assert_not_null(L);
    R = (double*) malloc( n2 * sizeof (double));
    assert_not_null(R);
    L_i = (int*) malloc( n1 * sizeof (int));
    assert_not_null(L_i);
    R_i = (int*) malloc( n2 * sizeof (int));
    assert_not_null(R_i);
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
    i = 0; /* Initial index of first subarray  */
    j = 0; /* Initial index of second subarray */
    k = l; /* Initial index of merged subarray */
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


void mergeSort(double* arr,int* arr_i, int l, int r)
{
    /* based on www.geeksforgeeks.org/merge-sort/ */
    /* l is for left index and r is right index of the
   sub-array of arr to be sorted */
    if (l < r) {
        /* Same as (l+r)/2, but avoids overflow for large l and h */
        int m = l + (r - l) / 2;

        /* Sort first and second halves */
        mergeSort(arr,arr_i, l, m);
        mergeSort(arr,arr_i, m + 1, r);

        merge(arr,arr_i, l, m, r);
    }
}

/**** Mergesort end ****/

/**********   Eigen_values START **********/
double* extract_eigen_values_from_mat(double** mat,int n){ /* get the eigen values of the diagonal of mat */
    int i;
    double* eigen_values;
    eigen_values = calloc(n,sizeof(double) );
    assert_not_null( eigen_values);
    for(i=0 ; i<n ; i++){
        eigen_values[i]= mat [i][i];
    }
    return eigen_values;
}


int  eigengap_huristic(Eigen eigen,int k){ /*find k by the gap between the sorted eigen values */
    int m,i;
    double  delta_i, max;
    if (k > 0)
        return k;

    m =  (int ) eigen.mat_size / 2; /* it said that k<n/2 */
    max = - 1 ; /* must be positive , as the values are sorted */
    k = 0;
    for (i = 0 ; i < m ; i++){
        delta_i = abs_d (eigen.values[i] - eigen.values[i+1]) ;
        if (delta_i > max)
        {
            max = delta_i ;
            k=i+1;
        }

    }
    return k;
}

void eigen_to_matrix(Eigen eigen, double** E,int n)
{ /* set the egien values and vectors into one matrix */

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

void free_eigen(Eigen eigen){ /*free the memory of the eigen strct */
    free_matrix(eigen.vectors,eigen.mat_size);
    free(eigen.ranks);
    free(eigen.values);
    eigen.ranks = NULL;
    eigen.values = NULL;

}

/**********   Eigen_values END **********/




/***************  SPK START ******************/

void start_wam(double** observations , int n, int d, double*** W)
{ /*start_wam flag, find  create_adj_mat */

    /*passing a pinter to the W */
    *W = create_matrix(n, n);
    create_adj_mat(observations,n,d,*W);

}
void start_ddg(double** observations , int n, int d, double*** D)
{ /*start_ddg flag, find  diagonal_degree_mat */
    double** W = create_matrix(n, n);
    *D = create_matrix(n, n);
    create_adj_mat(observations,n,d,W);
    create_diagonal_degree_mat(W,n,*D);
    free_matrix(W,n);
}

void start_lnorm(double** observations , int n, int d, double*** L)
{ /*l norm flag, get lnorm matrix */
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
{   /*get a symtric matrix, find eigen values and vbectors using the jacobi algoritem  */
    Eigen eigen = find_eigen_vectors_and_values(observations, n);
    inplace_transpose_mat(eigen.vectors,n,n); /*  now each row is a vector */
    *E  = create_matrix(n+1, n);
    eigen_to_matrix(eigen,*E,n);
    free_eigen(eigen);
}

double** Create_T_from(Eigen eigen,int n, int k){ /* create the T matrix form eigen struct */
    int i,j;
    double **U, **T;
    U  = create_matrix(n, k);
    for (i=0; i<n; i++){
        for (j=0; j<k; j++)  /* each  eigen vector is a column in U */
            U[i][j] = eigen.vectors[j][i];
    }

    T = (double  ** ) calloc(n , sizeof( double * ) );
    assert_not_null(T);
    renorm_matrix_rows(U, n,k, T);
    free_matrix(U,n);
    return T;
}

spk_results activate_flag(char* goal,double** observations , int k, int n, int d)
{   /* run all the flags, that are not spk */

    spk_results res;
    Eigen eigen;
    double **W ,**D ,**L, **E, **T;
    res.T_size = n;
    res.k = n;
    if (is_goal("wam")){
        start_wam(observations,n,d,&W);
        res.T = W;
      }
    if (is_goal("ddg")){
        start_ddg(observations,n,d,&D);
        res.T = D;
    }
    if (is_goal("lnorm")) {
        start_lnorm(observations,n,d,&L);
        res.T = L;
    }
    if (is_goal("jacobi")) {
        start_jacobi(observations,n,&E);
        res.T_size = n + 1; res.T = E;
    }
    if ( ! is_goal("spk")){
        print_mat(res.T, res.T_size ,n);
        return res;
    }
    /* implicit  goal == spk */
    W = create_matrix(n, n);
    create_adj_mat(observations,n,d,W);

    D = create_matrix(n, n);
    create_diagonal_degree_mat(W,n,D);
    D_sqrt(D,n);

    L = create_matrix(n, n);
    create_L_norm(D,W,n,L);
    eigen = find_eigen_vectors_and_values(L, n);
    inplace_transpose_mat(eigen.vectors,n,n); /*  now each row is an eigen  vector, easily reordered */
    mergeSort(eigen.values,eigen.ranks,0,n-1);
    re_order_matrix_by_indces(eigen.vectors, eigen.ranks, n);

    k = eigengap_huristic(eigen,k);
    T = Create_T_from(eigen,n,k);

    free_matrix(W,n);
    free_matrix(D,n);
    free_matrix(L,n);
    free_eigen(eigen);
    res.T = T;  res.k = k;
    return res;
}
/***************  SPK END ******************/
/***************  SPK END ******************/

/******** C Interface ******/
void load_string(char** str,char* cpy)
{ /* load string from args to var*/
    int len;
    len = ( (int) strlen (cpy)) + 1;
    *str = (char  *) malloc(  len * sizeof(char) );
    assert_not_null(*str);
    strcpy(*str, cpy);
}

int string_to_doubles(char *row,double* arr)
{ /* convert a string of doubles into array*/
    int i;
    char* ptr;
    i = 0;
    ptr = strtok(row, ",");
    while(ptr != NULL)
    {
        arr[i] = atof(ptr);
        ptr = 	strtok(NULL, ",");
        i++;
    }
    return i;


}

Tuple2 load_observations_from_file(double** observations, char* file_name)
{    /* read file and turn it into a 2d matrix of doubles*/
    int d,n,i;
    FILE *fp;
    char* row;
    Tuple2 sizes;
    d = 50;

    fp = fopen(file_name,"r");
    row = (char* ) calloc(1000, sizeof(char));
    assert_not_null(row);
    i =0;

    while(fscanf(fp,"%s",row)==1)
    { /* load data */
        d = string_to_doubles(row, observations[i]);
        i++;

    }
    n = i;
    fclose(fp);
    free(row);
    sizes.i = n;
    sizes.j = d;
    return sizes;

}






int main(int argc, char* argv[])
{   /* c interface of code */ 
    int k,n,d;
    char*  goal;
    char*  file_name;
    double** T_clusters_list;
    int * T_clusters_indexes;
    double** observations, **observations_load;
    Tuple2 sizes;
    spk_results Res;
    n = MAXSIZE ;  d = MAXSIZE;/*for jacobi d=50*/
    observations_load =  create_matrix(n, d);
    super_assert((argc==4) );
    super_assert( (atof(argv[1]) == atoi(argv[1])) ); /* is k an integer */

    k =  atoi(argv[1]);
    load_string(&goal,argv[2]);
    load_string(&file_name,argv[3]);

    assert_goal(goal);

    sizes = load_observations_from_file(observations_load, file_name);

    n = sizes.i;    d = sizes.j;
    observations =  create_matrix(n, d);
    copy_matrix(observations,observations_load,n,d);
    Res = activate_flag( goal, observations , k,  n, d);
    if(is_goal("spk")){
        k = Res.k;
        T_clusters_list = get_init_clusters_list(Res.T,k);
        T_clusters_indexes = init_clusters_indexes(k);

        simple_kmean(Res.T, T_clusters_list, T_clusters_indexes,n,k,k,0);
        free_matrix(T_clusters_list,k);
        free(T_clusters_indexes);

    }

    /* Free all */
    free_matrix(observations_load,MAXSIZE);
    free(goal);
    free(file_name);
    free_matrix(observations, n);
    free_matrix(Res.T,Res.T_size );
    return 0;
}






