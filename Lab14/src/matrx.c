#include "matrix.h"

matrix new_matrix ( const int rows , const int cols){
    matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    assert (rows >0);
    assert (cols >0);
    mat.data = ( double *) malloc ( sizeof(double) * rows * cols);

    for (int i=0; i <(rows * cols); i++){
        mat.data[i] = 0.0;
    }
    return mat;
}

void print_matrix_full ( const matrix * mat , char* varname ){
    assert (mat ->rows >0);
    assert (mat ->cols >0);

    printf ("\n %.100s =\n", &varname[0] );
    for(int i=1; i <=mat ->rows; i++ ) {
        printf (" | ");
        for(int j=1; j <=mat ->cols; j++)
        {
            printf (" %10.3e",mgetp (mat ,i,j));
            if (j<mat ->cols) { printf (", ");}
            else { printf (" ");}
        }
        printf ("|\n");
    }
    printf ("\n");
}

matrix matrix_add ( const matrix * A, const matrix * B){
    const int rows = A->rows;
    const int cols = A->cols;

    assert (rows ==B->rows);
    assert (cols ==B->cols);
    matrix C = new_matrix (rows ,cols);

    for (int i=1; i <= rows; i++)
        for (int j=1; j <= cols; j++)
        {
            mget(C,i,j) = mgetp (A,i,j) + mgetp (B,i,j);
        }

    return C;
}

matrix matrix_sub ( const matrix * A, const matrix * B){
    const int rows = A->rows;
    const int cols = A->cols;

    assert (rows ==B->rows);
    assert (cols ==B->cols);
    matrix C = new_matrix (rows ,cols);

    for (int i=1; i <= rows; i++)
        for (int j=1; j <= cols; j++)
        {
            mget(C,i,j) = mgetp (A,i,j) - mgetp (B,i,j);
        }

    return C;
}

matrix matrix_mult ( const matrix * A, const matrix * B){
    const int rowsA = A->rows;
    const int colsA = A->cols;
    const int rowsB = B->rows;
    const int colsB = B->cols;

    assert ( colsA == rowsB );
    matrix C = new_matrix (rowsA , colsB );

    for (int i=1; i <= rowsA ; i++) {
        for (int j=1; j <= colsB ; j++) {
            for (int k=1; k <= colsA ; k++) {
                mget(C,i,j) += mgetp (A,i,k)* mgetp (B,k,j);
            }
        }
    }

    return C;
}

matrix matrix_dot_mult ( const matrix * A, const matrix * B){
    const int rows = A->rows;
    const int cols = A->cols;

    assert (rows ==B->rows);
    assert (cols ==B->cols);
    matrix C = new_matrix (rows ,cols);

    for (int i=1; i <= rows; i++){
        for (int j=1; j <= cols; j++){
            mget(C,i,j) = mgetp (A,i,j)* mgetp (B,i,j);
        }
    }

    return C;
}

matrix matrix_transpose ( const matrix * A){
    const int rows = A->rows;
    const int cols = A->cols;

    matrix At = new_matrix (cols ,rows);

    for (int i=1; i <= rows; i++){
        for (int j=1; j <= cols; j++){
            mget(At,j,i) = mgetp (A,i,j);
        }
    }

    return At;
}

void delete_matrix ( matrix * mat){
    free (mat ->data);
    mat ->data = NULL;
    mat ->rows = 0;
    mat ->cols = 0;
}

vector new_vector ( const int size){
    vector vec;
    vec.size = size;
    assert (size >0);

    vec.data = ( double *) malloc ( sizeof(double) * size);
    for (int i=0; i <size; i++){
        vec.data[i] = 0.0;
    }
    return vec;
}

void print_vector_full ( const vector * vec , char* varname ){
    assert (vec ->size >0);

    printf ("\n");
    printf (" %.100s =\n", & varname[0] );
    printf (" | ");

    for(int i=1; i <=vec ->size; i++ ) {
        printf (" %10.3e",vgetp (vec ,i));
        if (i<vec ->size) { printf (", ");}
    }

    printf (" |^T\n\n");
}

vector vector_add ( const vector * x, const vector * y){
    const int size = x->size;
    assert (size ==y->size);
    vector z = new_vector (size);

    for (int i=1; i <= size; i++) {
        vget(z,i) = vgetp (x,i) + vgetp (y,i);
    }

    return z;
}

vector vector_sub ( const vector * x, const vector * y){
    const int size = x->size;
    assert (size ==y->size);
    vector z = new_vector (size);

    for (int i=1; i <= size; i++) {
        vget(z,i) = vgetp (x,i) - vgetp (y,i);
    }

    return z;
}

double vector_dot_mult ( const vector * x, const vector * y){
    const int size = x->size; 
    assert (size ==y->size);

    double z = 0.0;
    for (int i=1; i <= size; i++){
        z += vgetp (x,i) * vgetp (y,i); 
    }

    return z;
}

void delete_vector ( vector * vec){
    free (vec ->data);
    vec ->data = NULL;
    vec ->size = 0;
}

void print_scalar_full ( const double * z, char* varname ){

}

vector matrix_vector_mult ( const matrix * A, const vector * x){
    const int rows = A->rows; const int cols = A->cols;
    const int size = x->size;
    assert (cols == size);
    vector Ax = new_vector (rows);

    for (int i=1; i <= rows; i++){
        double tmp = 0.0;
        for (int j=1; j <= size; j++){
            tmp += mgetp(A,i,j)* vgetp (x,j); 
        }
        vget(Ax ,i) = tmp;
    }

    return Ax;
}

vector solve ( const matrix * A, const vector * b){
    const int rows = A->rows;
    const int cols = A->cols;
    const int size = b->size;

    assert (rows == cols); assert (rows == size);
    vector x = new_vector (rows);

    for (int i=1; i <=( size -1); i++) {
        // Select largest pivot in current column
        int p=i; double maxA = -100.0;
        for (int j=i; j <= size; j++) {
            double tmp = fabs( mgetp (A,j,i));
            if ( tmp > maxA ){ p = j; maxA = tmp; }
        }

        // See if matrix is singular
        if (maxA <= 1.0e-14){
            printf (" Cannot invert system \n"); 
            exit (1); 
        }

        // Pivot (aka interchange rows)
        if (p!=i) {
            for (int j=1; j <= size; j++) {
                double tmp1 = mgetp (A,i,j);
                mgetp (A,i,j) = mgetp (A,p,j); mgetp (A,p,j) = tmp1;
            }

            double tmp2 = vgetp (b,i);
            vgetp (b,i) = vgetp (b,p); vgetp (b,p) = tmp2;
        }

        // Eliminate below diagonal
        for (int j=(i+1); j <= size; j++) {
            double dm = mgetp (A,j,i)/ mgetp (A,i,i);
            for (int k=1; k <= size; k++) {
                mgetp (A,j,k) = mgetp (A,j,k) - dm* mgetp (A,i,k); 
            }

            vgetp (b,j) = vgetp (b,j) - dm* vgetp (b,i);
        }
    }

    // Backward substitution

    vget(x,size) = vgetp (b,size)/ mgetp (A,size ,size);
    for (int j=1; j <=( size -1); j++){
        double sum = 0.0;

        for (int k=( size -j+1); k <= size; k++){
            sum = sum + mgetp (A,size -j,k)*vget(x,k); 
        }

        vget(x,size -j) = ( vgetp (b,size -j) - sum) / mgetp (A,size -j,size -j);
    }

    return x;
}

double power_iteration(const matrix * A, vector * b, const double tol, const int max_iter){
    double lambda_old = 0.0, lambda_new = 0.0;
    double diff;
    int iter = 0;

    // Normalize initial vector b
    double norm_b = sqrt(vector_dot_mult(b, b));
    for (int i = 1; i <= b->size; i++) vgetp(b, i) /= norm_b;

    do {
        // Multiply A * b
        vector Ab = matrix_vector_mult(A, b);

        // Compute new eigenvalue estimate λ = (bᵀ A b)
        lambda_new = vector_dot_mult(b, &Ab);

        // Normalize new vector
        double norm_Ab = sqrt(vector_dot_mult(&Ab, &Ab));
        for (int i = 1; i <= b->size; i++) vgetp(b, i) = vgetp(b, i) / norm_Ab;

        // Check for convergence
        diff = fabs(lambda_new - lambda_old);
        lambda_old = lambda_new;

        delete_vector(&Ab);
        iter++;

    } while (diff > tol && iter < max_iter);

    return lambda_new;       
}

double inverse_iteration(const matrix * A, vector * b, const double mu, const double tol, const int max_iter){
    int n = b->size;
    double lambda_old = mu, lambda_new = mu;
    double diff;
    int iter = 0;

    // Build (A - μI)
    matrix A_shift = new_matrix(A->rows, A->cols);
    for (int i = 1; i <= A->rows; i++) {
        for (int j = 1; j <= A->cols; j++) {
            double val = mgetp(A, i, j);
            if (i == j) val -= mu;
            mgetp((&A_shift), i, j) = val;
        }
    }

    // Normalize initial vector
    double norm_b = sqrt(vector_dot_mult(b, b));
    for (int i = 1; i <= n; i++) vgetp(b, i) /= norm_b;

    do {
        // Solve (A - μI) y = b
        vector y = solve(&A_shift, b);

        // Normalize y
        double norm_y = sqrt(vector_dot_mult(&y, &y));
        for (int i = 1; i <= n; i++) vgetp(b, i) = vgetp((&y), i) / norm_y;

        // Rayleigh quotient estimate: λ ≈ bᵀ A b
        vector c = matrix_vector_mult(A, b);
        lambda_new = vector_dot_mult(b, (&c));
        delete_vector(&c);

        diff = fabs(lambda_new - lambda_old);
        lambda_old = lambda_new;

        delete_vector(&y);
        iter++;

    } while (diff > tol && iter < max_iter);

    delete_matrix(&A_shift);

    return lambda_new;
}
