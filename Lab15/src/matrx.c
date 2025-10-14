#include "matrix.h"

/**
 * Allocates a zero matrix of size rows x cols
 */
matrix new_matrix(const int rows, const int cols) {
    matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    assert(rows > 0 && cols > 0);
    mat.data = (double*)malloc(sizeof(double) * rows * cols);
    for (int i = 0; i < rows * cols; i++)
        mat.data[i] = 0.0;
    return mat;
}

/**
 * Prints matrix contents with variable name
 */
void print_matrix_full(const matrix* mat, char* varname) {
    assert(mat->rows > 0 && mat->cols > 0);
    printf("\n %.100s =\n", varname);
    for (int i = 1; i <= mat->rows; i++) {
        printf(" | ");
        for (int j = 1; j <= mat->cols; j++) {
            printf(" %10.3e", mgetp(mat, i, j));
            if (j < mat->cols) printf(", ");
            else printf(" ");
        }
        printf("|\n");
    }
    printf("\n");
}

/**
 * Returns the sum of matrices A and B
 */
matrix matrix_add(const matrix* A, const matrix* B) {
    const int rows = A->rows, cols = A->cols;
    assert(rows == B->rows && cols == B->cols);
    matrix C = new_matrix(rows, cols);
    for (int i = 1; i <= rows; i++)
        for (int j = 1; j <= cols; j++)
            mget(C, i, j) = mgetp(A, i, j) + mgetp(B, i, j);
    return C;
}

/**
 * Returns the difference A - B
 */
matrix matrix_sub(const matrix* A, const matrix* B) {
    const int rows = A->rows, cols = A->cols;
    assert(rows == B->rows && cols == B->cols);
    matrix C = new_matrix(rows, cols);
    for (int i = 1; i <= rows; i++)
        for (int j = 1; j <= cols; j++)
            mget(C, i, j) = mgetp(A, i, j) - mgetp(B, i, j);
    return C;
}

/**
 * Returns A * B (matrix-matrix multiplication)
 */
matrix matrix_mult(const matrix* A, const matrix* B) {
    const int rowsA = A->rows, colsA = A->cols, rowsB = B->rows, colsB = B->cols;
    assert(colsA == rowsB);
    matrix C = new_matrix(rowsA, colsB);
    for (int i = 1; i <= rowsA; i++)
        for (int j = 1; j <= colsB; j++)
            for (int k = 1; k <= colsA; k++)
                mget(C, i, j) += mgetp(A, i, k) * mgetp(B, k, j);
    return C;
}

/**
 * Returns elementwise product (Hadamard product) of A and B
 */
matrix matrix_dot_mult(const matrix* A, const matrix* B) {
    const int rows = A->rows, cols = A->cols;
    assert(rows == B->rows && cols == B->cols);
    matrix C = new_matrix(rows, cols);
    for (int i = 1; i <= rows; i++)
        for (int j = 1; j <= cols; j++)
            mget(C, i, j) = mgetp(A, i, j) * mgetp(B, i, j);
    return C;
}

/**
 * Returns the transpose of A
 */
matrix matrix_transpose(const matrix* A) {
    const int rows = A->rows, cols = A->cols;
    matrix At = new_matrix(cols, rows);
    for (int i = 1; i <= rows; i++)
        for (int j = 1; j <= cols; j++)
            mget(At, j, i) = mgetp(A, i, j);
    return At;
}

/**
 * Frees a matrix's memory
 */
void delete_matrix(matrix* mat) {
    free(mat->data);
    mat->data = NULL;
    mat->rows = mat->cols = 0;
}

/**
 * Allocates a zero vector of given size
 */
vector new_vector(const int size) {
    vector vec;
    assert(size > 0);
    vec.size = size;
    vec.data = (double*)malloc(sizeof(double) * size);
    for (int i = 0; i < size; i++)
        vec.data[i] = 0.0;
    return vec;
}

/**
 * Prints a vector with variable name
 */
void print_vector_full(const vector* vec, char* varname) {
    assert(vec->size > 0);
    printf("\n %.100s =\n", varname);
    printf(" | ");
    for (int i = 1; i <= vec->size; i++) {
        printf(" %10.3e", vgetp(vec, i));
        if (i < vec->size) printf(", ");
    }
    printf(" |^T\n\n");
}

/**
 * Returns x + y (vector addition)
 */
vector vector_add(const vector* x, const vector* y) {
    const int size = x->size;
    assert(size == y->size);
    vector z = new_vector(size);
    for (int i = 1; i <= size; i++)
        vget(z, i) = vgetp(x, i) + vgetp(y, i);
    return z;
}

/**
 * Returns x - y (vector subtraction)
 */
vector vector_sub(const vector* x, const vector* y) {
    const int size = x->size;
    assert(size == y->size);
    vector z = new_vector(size);
    for (int i = 1; i <= size; i++)
        vget(z, i) = vgetp(x, i) - vgetp(y, i);
    return z;
}

/**
 * Returns the dot product of x and y
 */
double vector_dot_mult(const vector* x, const vector* y) {
    const int size = x->size;
    assert(size == y->size);
    double z = 0.0;
    for (int i = 1; i <= size; i++)
        z += vgetp(x, i) * vgetp(y, i);
    return z;
}

/**
 * Frees a vector's memory
 */
void delete_vector(vector* vec) {
    free(vec->data);
    vec->data = NULL;
    vec->size = 0;
}

/**
 * Prints a scalar value with variable name (not used in this code, but provided for completeness)
 */
void print_scalar_full(const double* z, char* varname) {
    printf("%s = %10.3e\n", varname, *z);
}

/**
 * Returns A*x (matrix-vector multiplication)
 */
vector matrix_vector_mult(const matrix* A, const vector* x) {
    const int rows = A->rows, cols = A->cols, size = x->size;
    assert(cols == size);
    vector Ax = new_vector(rows);
    for (int i = 1; i <= rows; i++) {
        double tmp = 0.0;
        for (int j = 1; j <= size; j++)
            tmp += mgetp(A, i, j) * vgetp(x, j);
        vget(Ax, i) = tmp;
    }
    return Ax;
}

/**
 * Solves Ax = b for x using Gaussian elimination with partial pivoting.
 * (A and b are overwritten during solution!)
 */
vector solve(const matrix* A, const vector* b) {
    const int rows = A->rows, cols = A->cols, size = b->size;
    assert(rows == cols && rows == size);

    // Copy A and b to local storage (since overwritten)
    matrix Acopy = new_matrix(rows, cols);
    for (int i = 1; i <= rows; i++)
        for (int j = 1; j <= cols; j++)
            mget(Acopy, i, j) = mgetp(A, i, j);
    vector bcopy = new_vector(size);
    for (int i = 1; i <= size; i++)
        vget(bcopy, i) = vgetp(b, i);

    vector x = new_vector(rows);

    // Forward elimination with partial pivoting
    for (int i = 1; i <= size - 1; i++) {
        // Find largest pivot in current column
        int p = i;
        double maxA = -100.0;
        for (int j = i; j <= size; j++) {
            double tmp = fabs(mget(Acopy, j, i));
            if (tmp > maxA) { p = j; maxA = tmp; }
        }
        // Check for near singular matrix
        if (maxA <= 1.0e-14) {
            printf("Cannot invert system\n");
            exit(1);
        }
        // Pivot rows
        if (p != i) {
            for (int j = 1; j <= size; j++) {
                double tmp1 = mget(Acopy, i, j);
                mget(Acopy, i, j) = mget(Acopy, p, j);
                mget(Acopy, p, j) = tmp1;
            }
            double tmp2 = vget(bcopy, i);
            vget(bcopy, i) = vget(bcopy, p);
            vget(bcopy, p) = tmp2;
        }
        // Elimination
        for (int j = i + 1; j <= size; j++) {
            double dm = mget(Acopy, j, i) / mget(Acopy, i, i);
            for (int k = 1; k <= size; k++)
                mget(Acopy, j, k) -= dm * mget(Acopy, i, k);
            vget(bcopy, j) -= dm * vget(bcopy, i);
        }
    }

    // Backward substitution
    vget(x, size) = vget(bcopy, size) / mget(Acopy, size, size);
    for (int j = 1; j <= size - 1; j++) {
        double sum = 0.0;
        for (int k = size - j + 1; k <= size; k++)
            sum += mget(Acopy, size - j, k) * vget(x, k);
        vget(x, size - j) = (vget(bcopy, size - j) - sum) / mget(Acopy, size - j, size - j);
    }

    delete_matrix(&Acopy);
    delete_vector(&bcopy);
    return x;
}

/**
 * Computes the dominant eigenvalue using power iteration.
 * Returns the eigenvalue, and b converges to the eigenvector.
 */
double power_iteration(const matrix* A, vector* b, const double tol, const int max_iter) {
    double lambda_old = 0.0, lambda_new = 0.0, diff;
    int iter = 0;
    // Normalize initial vector
    double norm_b = sqrt(vector_dot_mult(b, b));
    for (int i = 1; i <= b->size; i++) vgetp(b, i) /= norm_b;

    do {
        // Ab = A * b
        vector Ab = matrix_vector_mult(A, b);
        lambda_new = vector_dot_mult(b, &Ab); // Rayleigh quotient
        // Normalize
        double norm_Ab = sqrt(vector_dot_mult(&Ab, &Ab));
        for (int i = 1; i <= b->size; i++) vgetp(b, i) = vgetp(b, i) / norm_Ab;
        // Check convergence
        diff = fabs(lambda_new - lambda_old);
        lambda_old = lambda_new;
        delete_vector(&Ab);
        iter++;
    } while (diff > tol && iter < max_iter);

    return lambda_new;
}

/**
 * Computes the eigenvalue closest to mu using inverse iteration.
 * Returns the eigenvalue, and b converges to the eigenvector.
 */
double inverse_iteration(const matrix* A, vector* b, const double mu, const double tol, const int max_iter) {
    int n = b->size;
    double lambda_old = mu, lambda_new = mu, diff;
    int iter = 0;
    // Build (A - mu I)
    matrix A_shift = new_matrix(A->rows, A->cols);
    for (int i = 1; i <= A->rows; i++)
        for (int j = 1; j <= A->cols; j++) {
            double val = mgetp(A, i, j);
            if (i == j) val -= mu;
            mget(A_shift, i, j) = val;
        }
    // Normalize initial vector
    double norm_b = sqrt(vector_dot_mult(b, b));
    for (int i = 1; i <= n; i++) vgetp(b, i) /= norm_b;

    do {
        // y = solve (A - mu I) y = b
        vector y = solve(&A_shift, b);
        // Normalize y
        double norm_y = sqrt(vector_dot_mult(&y, &y));
        for (int i = 1; i <= n; i++) vgetp(b, i) = vget(y, i) / norm_y;
        // Rayleigh quotient
        vector c = matrix_vector_mult(A, b);
        lambda_new = vector_dot_mult(b, &c);
        delete_vector(&c);
        diff = fabs(lambda_new - lambda_old);
        lambda_old = lambda_new;
        delete_vector(&y);
        iter++;
    } while (diff > tol && iter < max_iter);

    delete_matrix(&A_shift);
    return lambda_new;
}

/**
 * Performs QR decomposition of A (A = QR), using Modified Gram-Schmidt process.
 * Q and R must be allocated before calling (as n x n matrices).
 */
void qr_decompose(const matrix* A, matrix* Q, matrix* R) {
    int n = A->rows;
    assert(A->rows == A->cols);
    // Temporary local matrices
    matrix q = new_matrix(n, n);
    matrix r = new_matrix(n, n);
    matrix a = new_matrix(n, n);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            mget(a, i, j) = mgetp(A, i, j);
    // Modified Gram-Schmidt
    for (int k = 1; k <= n; k++) {
        // Copy column k
        for (int i = 1; i <= n; i++)
            mget(q, i, k) = mget(a, i, k);
        for (int j = 1; j < k; j++) {
            double dot = 0.0;
            for (int i = 1; i <= n; i++)
                dot += mget(q, i, j) * mget(a, i, k);
            mget(r, j, k) = dot;
            for (int i = 1; i <= n; i++)
                mget(q, i, k) -= dot * mget(q, i, j);
        }
        double norm = 0.0;
        for (int i = 1; i <= n; i++)
            norm += mget(q, i, k) * mget(q, i, k);
        norm = sqrt(norm);
        mget(r, k, k) = norm;
        for (int i = 1; i <= n; i++)
            mget(q, i, k) /= norm;
    }
    // Copy to output
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++) {
            mgetp(Q, i, j) = mget(q, i, j);
            mgetp(R, i, j) = mget(r, i, j);
        }
    delete_matrix(&q);
    delete_matrix(&r);
    delete_matrix(&a);
}

/**
 * Computes all eigenvalues of a real symmetric matrix using QR iteration.
 * A is overwritten and converges to nearly diagonal; eigenvalues are on the diagonal.
 * No output is printed; you can print A after to see eigenvalues.
 */
void qr_eigenvalues(matrix* A, double tol, int max_iter) {
    int n = A->rows;
    assert(A->rows == A->cols);

    matrix Ak = new_matrix(n, n);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            mget(Ak, i, j) = mgetp(A, i, j);

    matrix Q = new_matrix(n, n);
    matrix R = new_matrix(n, n);
    matrix Ak_prev = new_matrix(n, n);

    for (int iter = 1; iter <= max_iter; iter++) {
        // Save previous Ak
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                mget(Ak_prev, i, j) = mget(Ak, i, j);

        qr_decompose(&Ak, &Q, &R);

        // Ak+1 = R * Q
        matrix RQ = matrix_mult(&R, &Q);
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                mget(Ak, i, j) = mget(RQ, i, j);

        // Check convergence: Frobenius norm of difference
        double delta = 0.0;
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++) {
                double d = mget(Ak, i, j) - mget(Ak_prev, i, j);
                delta += d * d;
            }
        delta = sqrt(delta);

        delete_matrix(&RQ);
        if (delta < tol) break;
    }

    // Copy result back to input matrix (A becomes nearly diagonal)
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            mgetp(A, i, j) = mget(Ak, i, j);

    delete_matrix(&Ak);
    delete_matrix(&Ak_prev);
    delete_matrix(&Q);
    delete_matrix(&R);
}