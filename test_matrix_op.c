#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#import <time.h>

 double[][] create_matrix(int N) {
    double mat[N][N];
    int i, j;
    srand(time(NULL));
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
            mat[i][j] =i+j;
    }
    return mat;

}


int main() {
    int N

}
