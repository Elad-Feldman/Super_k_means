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
tuple2 load_observations_from_file(double** observations,char* file_name)
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
    tuple2 sizes;
    sizes.i = N;
    sizes.j = d;
    return sizes;

}

int main(int argc, char* argv[])
{
    int k,N,d;
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
    assert(k>0);

    assert_goal(goal);
    printf("k: %d\n goal %s\n init: %s\n file_name:%s\n ",k,goal,path_init_centroids,file_name);
    tuple2 sizes = load_observations_from_file(observations, file_name);
    N=sizes.i;
    d=sizes.j;

    print_mat(observations,N,d);
    // Free all

    free(goal);
    free(path_init_centroids);
    free(file_name);
    free_matrix(observations, N);




}