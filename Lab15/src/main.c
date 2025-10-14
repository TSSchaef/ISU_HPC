#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include "matrix.h"

/**
 * Main program: demonstrates matrix construction, solve, eigenvalue algorithms
 */
int main(int argc, char *argv[]) {
    int n;
    printf("Enter the size of the square matrix: \n");

    if(scanf("%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "Invalid input. Please enter a positive integer.\n");
        return 1;
    }

    // Create a random upper-triangular matrix L with 1 on diagonal
    matrix L = new_matrix(n, n);
    srand((unsigned int)time(NULL));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                mget(L, i, j) = 1.0;
            } else if (j > i) {
                mget(L, i, j) = ((double)rand()) / INT_MAX;
            } else {
                mget(L, i, j) = 0.0;
            }
        }
    }
    print_matrix_full(&L, "L");

    // Compute L-transpose
    matrix Lt = matrix_transpose(&L);
    print_matrix_full(&Lt, "L-transpose");

    // Compute A = L * L^T (symmetric positive definite)
    matrix A = matrix_mult(&L, &Lt);
    print_matrix_full(&A, "A = L * L^T");

    // Create and print random vector b
    /*printf("Random vector b:\n");
    vector b = new_vector(n);
    for (int i = 1; i <= n; i++)
        vget(b, i) = ((double)rand()) / INT_MAX;
    print_vector_full(&b, "b");

    // Solve Ax = b
    vector x = solve(&A, &b);
    printf("Solution vector x:\n");
    print_vector_full(&x, "x"); */

    /*
    // Power iteration: dominant eigenvalue
    printf("Power Iteration:\n");
    printf("%lf\n", power_iteration(&A, &b, 1e-6, 1000));

    // Inverse iteration: eigenvalue near shift mu
    printf("Inverse Iteration:\n");
    printf("%lf\n", inverse_iteration(&A, &b, 0.5, 1e-6, 1000));
    */

    // QR iteration: all eigenvalues
    printf("QR Iteration for Eigenvalues (all eigenvalues):\n");
    matrix A_qr = new_matrix(n, n);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            mget(A_qr, i, j) = mget(A, i, j);

    qr_eigenvalues(&A_qr, 1e-8, 100);

    printf("A after QR iterations (should be nearly diagonal):\n");
    print_matrix_full(&A_qr, "A_qr");

    printf("Eigenvalues (diagonal of A_qr):\n");
    for (int i = 1; i <= n; i++)
        printf("%.6e\n", mget(A_qr, i, i));

    // Clean up all memory
    delete_matrix(&L);
    delete_matrix(&Lt);
    delete_matrix(&A);
    delete_matrix(&A_qr);
    //delete_vector(&b);
    //delete_vector(&x);

    return 0;
}
