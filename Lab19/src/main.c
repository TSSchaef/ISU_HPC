#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <fine|coarse> <vector_size>\n", argv[0]);
        return 1;
    }

    char *mode = argv[1];
    int N = atoi(argv[2]);
    double *v = malloc(N * sizeof(double));
    double *u = malloc(N * sizeof(double));
    if (!v || !u) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // Initialize vector with some values
    for (int i = 0; i < N; i++) v[i] = i + 1.0;

    double norm_sq = 0.0;

    double t_start = omp_get_wtime();

    if (strcmp(mode, "fine") == 0) {
        // ---------- FINE-GRAIN PARALLELISM ----------
        #pragma omp parallel for reduction(+:norm_sq)
        for (int i = 0; i < N; i++) {
            norm_sq += v[i] * v[i];
        }
        double norm = sqrt(norm_sq);

        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            u[i] = v[i] / norm;
        }

    } else if (strcmp(mode, "coarse") == 0) {
        // ---------- COARSE-GRAIN PARALLELISM ----------
        int num_threads;
        double norm_sq_local = 0.0;

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
            #pragma omp single
            num_threads = nthreads;

            int chunk = (N + nthreads - 1) / nthreads;
            int start = tid * chunk;
            int end = (start + chunk > N) ? N : start + chunk;

            double local_sum = 0.0;
            for (int i = start; i < end; i++) {
                local_sum += v[i] * v[i];
            }

            #pragma omp atomic
            norm_sq += local_sum;
        }

        double norm = sqrt(norm_sq);

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
            int chunk = (N + nthreads - 1) / nthreads;
            int start = tid * chunk;
            int end = (start + chunk > N) ? N : start + chunk;

            for (int i = start; i < end; i++) {
                u[i] = v[i] / norm;
            }
        }

    } else {
        printf("Invalid mode: %s. Use 'fine' or 'coarse'.\n", mode);
        free(v);
        free(u);
        return 1;
    }

    double t_end = omp_get_wtime();

    printf("Mode: %s | N: %d | Time: %.6f s\n", mode, N, t_end - t_start);
    printf("First 5 normalized values: ");
    for (int i = 0; i < (N < 5 ? N : 5); i++) {
        printf("%.4f ", u[i]);
    }
    printf("\n");

    free(v);
    free(u);
    return 0;
}

