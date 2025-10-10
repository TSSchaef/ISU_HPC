#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

#include "matrix.h"


int main(int argc, char *argv[]) {
    int n;
    printf("Enter the size of the square matrix: \n");
    scanf("%d", &n);

    matrix L = new_matrix(n, n);
    
    srand((unsigned int)time(NULL));

    int i, j;
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            if (i == j) {
                mget(L, i, j) = 1.0; // Diagonal elements are 1
            } else if (j > i) {
                mget(L, i, j) = ((double)rand()) / INT_MAX; // Upper triangular part random double 0-1
            } else {
                mget(L, i, j) = 0.0; // Lower triangular part is 0
            }
        }
    }
    print_matrix_full(&L, "L");

    matrix Lt = matrix_transpose(&L);
    print_matrix_full(&Lt, "L-transpose");

    //Generate symetric positive definite matrix A
    matrix A = matrix_mult(&Lt, &L);
    print_matrix_full(&A, "A = L^T * L");

    printf("Random vector b:\n");
    vector b = new_vector(n);

    for(i = 1; i <= n; i++) {
        vget(b, i) = ((double)rand()) / INT_MAX; // Random double 0-1
    }
    print_vector_full(&b, "b");

    // Solve Ax = b
    vector x = solve(&A, &b);
    printf("Solution vector x:\n");
    print_vector_full(&x, "x");

    // power iteration and inverse iteration
    printf("Power Iteration:\n");
    printf("%lf\n", power_iteration(&A, &b, 1e-6, 1000));

    printf("Inverse Iteration:\n");
    printf("%lf\n", inverse_iteration(&A, &b, 0.5, 1e-6, 1000));

    // Clean up data structures
    delete_matrix(&L);
    delete_matrix(&Lt);
    delete_matrix(&A);

    delete_vector(&b);
    delete_vector(&x);

    return 0;
}
