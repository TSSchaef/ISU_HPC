/*
 *
 * Compute w(i) = sum_{j=1..K} A(i,j) * x(j) without storing A or x,
 * where A(i,j) = 1/(i + j - 1)  (Hilbert-like) and x(j) = 1/j.
 *
 * So: w(i) = sum_{j=1..K} (1.0 / (i + j - 1)) * (1.0 / j)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s N K num_threads\n", argv[0]);
        return 1;
    }

    long N = atol(argv[1]);
    long K = atol(argv[2]);
    int num_threads = atoi(argv[3]);

    if (N <= 0 || K <= 0 || num_threads <= 0) {
        fprintf(stderr, "All arguments must be positive: N K num_threads\n");
        return 1;
    }

    double *w = (double*) malloc(sizeof(double) * (size_t)N);
    if (!w) {
        fprintf(stderr, "Allocation failed for output vector w (size %ld)\n", N);
        return 1;
    }

    omp_set_num_threads(num_threads);

    double t0 = omp_get_wtime();

    /* Parallelize over rows i. Each w[i] is independent. */
    #pragma omp parallel for schedule(static)

    for (long ii = 0; ii < N; ++ii) {
        double sum = 0.0;
        long i = ii + 1; /* convert to 1-based index for formula */
        for (long j = 1; j <= K; ++j) {
            sum += (1.0 / (double)(i + j - 1)) * (1.0 / (double)j);
        }
        w[ii] = sum;
    }

    double t1 = omp_get_wtime();
    double elapsed = t1 - t0;

    /* Print first and last elements and timing */
    if (N >= 1) {
        printf("N=%ld K=%ld threads=%d time(s)=%.6f\n", N, K, num_threads, elapsed);
        printf("w[0]   = %.15e\n", w[0]);
        printf("w[%ld] = %.15e\n", N-1, w[N-1]);
    }

    free(w);
    return 0;
}
