#include <Python.h>
#include "matrix_op.c"

// gcc spkmeans.c && gcc  -o spkmeans spkmeans.c && spkmeans  5 spk  in1.txt

#define Arr_size(x)  (sizeof(x) / sizeof((x)[0]))

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
    int d,N,i;
    FILE *fp;
    fp = fopen(file_name,"r");
    char* row = (char* )calloc(1000, sizeof(char));
    i =0;
    while(fscanf(fp,"%s",row)==1){ // load data
        d = string_to_doubles(row, observations[i]);
        i++;
    }
    N=i;

    // change the size of observations to match the file
    observations = (double **) realloc(observations,N* sizeof(double *));
    assert(observations);
    for (i = 0; i < N; i++)
    {
        observations[i] = (double *) realloc(observations[i],d* sizeof(double));
        assert(observations[i]);
    }
    fclose(fp);
    free(row);
    Tuple2 sizes;
    sizes.i = N;
    sizes.j = d;
    return sizes;

}
//ns stands for negtive squre root this function returns D^(-1/2)
double** create_diagonal_degree_mat_ns(double** adj_mat,int N) {
    assert(adj_mat);
    assert(N>0);
    double** mat = create_matrix(N,N);
    int i;
    double sum_of_weights;//sum of weights per row
    sum_of_weights = 0;
    for (i = 0; i < N; i++)
    {
        sum_of_weights= sum_vector(adj_mat[i],N);
        assert(sum_of_weights > 0);
        mat[i][i] = 1 / sqrt( sum_of_weights);

    }
    return mat;
}
// subtracts B from A and returns the result
double** matrix_subtraction(double** A, double** B,int N) {
    double** result_mat = create_matrix(N, N);
    int i, j;
    for ( i = 0; i < N; i++)
    {

        for (j = 0; j < N; j++)
        {
            result_mat[i][j] = A[i][j] - B[i][j];
        }
    }
    return result_mat;
}

double** calculate_L_norm(double** D,double** W,int N) { /* TODO what happen to memory here? not clear */
    double** after_sub_mat;
    double** id_mat = create_Id_matrix(N);
    double** mat = mult_matrix(D, W, N);
    mat = mult_matrix(mat, D,N);
    after_sub_mat = matrix_subtraction(id_mat, mat, N);
    free(mat);
    free(id_mat);
    return after_sub_mat;
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
double** create_rotation_mat(double** A,int N)
{
    int i,j;
    double c,t,s, theta;
    double** P = create_Id_matrix(N);
    Tuple2 max_indces = find_ind_max_ele_off_diag(A, N);
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
    print_mat(P,N,N);
    return P;
}
double* extract_eigen_values_from_mat(double** mat,int N){
    int i;
    double * eigen_values = calloc(N,sizeof(double));
    for(i=0;i<N;i++){
        eigen_values[i]=mat[i][i];
    }
    return eigen_values;
}
double sum_square_elements_off_diag(double ** A,int N){
    int i,j;
    double sum=0;
    for (i = 0; i <N ;i++) {
        for (j = 0; j <N ;j++) {
            if(i!=j){
                sum+=A[i][j]*A[i][j];
            }

        }
    }
    return sum;
}
int check_convergence(double** A,double** A1,int N){
    double sum_A =sum_square_elements_off_diag(A, N);
    double sum_A1 = sum_square_elements_off_diag(A1,N);
    if(sum_A-sum_A1<=0.001){
        return 1;
    }else{
        return 0;
    }
}

Eigen find_eigen_vectors(double** A, int N){
    /* Start with A=L_norm */
    assert(A);
    assert(N>0);
    double** V = create_Id_matrix(N);
    double **P,** P_T;
    double** A1;
    double** A2;
    double** V1;
    int convergence =0;
    Eigen eigen;
    while( !convergence){
        /*print_mat(A,N,N);*/
        printf("\n");
        P = create_rotation_mat(A,N);
        P_T= transpose_mat(P,N,N);
        A1 = mult_matrix(P_T,A,N);
        A2 = mult_matrix(A1,P,N);

        free(A1);
        convergence= check_convergence(A,A2,N);
        free(A);
        A=A2;
        print_mat(A,N,N);
        A2=NULL;
        V1 = mult_matrix(V,P,N);
        free(V);
        V=V1;
        V1=NULL;
        free_matrix(P,N);
        free_matrix(P_T,N);
    }
    eigen.values = extract_eigen_values_from_mat(A, N);
    eigen.vectors = V;
    eigen.mat_size=N;
    return eigen;
}
void sort_eigen_values(Eigen eign_obj){
    double** mat_T = transpose_mat(eign_obj.vectors, eign_obj.mat_size);

}

int main(int argc, char* argv[])
{
    test_mat_op();
   /* int k,N,d;
    N = 1000;
    d = 10;
    char*  goal;
    char*  path_init_centroids;
    char*  file_name;
    double** observations = create_matrix(N,d);
    if (argc==5)
    {
        //TODO assert k is an integer
        k = atoi(argv[1]);
        load_string(&goal,argv[2]);
        load_string(&path_init_centroids,argv[3]);
        load_string(&file_name,argv[4]);

    }

    if (argc==4)
    {
        k = atoi(argv[1]);
        load_string(&goal,argv[2]);
        load_string(&path_init_centroids,"");
        load_string(&file_name,argv[3]);
    }
    // TODO if argc not 4 or 5
    //assert(k>0);

    assert_goal(goal);
    printf("k: %d\n goal %s\n init: %s\n file_name:%s\n ",k,goal,path_init_centroids,file_name);
    Tuple2 sizes = load_observations_from_file(observations, file_name);
    N=sizes.i;
    d=sizes.j;

    print_mat(observations,N,d);
    // Free all

    free(goal);
    free(path_init_centroids);
    free(file_name);
    free_matrix(observations, N);

    */


}