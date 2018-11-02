#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include<sys/resource.h>

#define N 2000
#define RND

void lup_doolittle(double *a, int *p, int n) {
    int max_pivot_row;
    int i, j, k;
    double max_pivot, tmp;
    for(i = 0; i < n; i++)
        p[i] = i;
    for(i = 0; i < n; i++) {
        max_pivot_row = i;
        max_pivot = *(a + (i * n + i));
        for(j = i + 1; j < n; j++) {
            tmp = fabs(*(a + (j * n + i)));
            if(tmp > max_pivot) {
                max_pivot_row = j;
                max_pivot = tmp;
            }
        }
        if(i != max_pivot_row) {
            tmp = p[i];
            p[i] = p[max_pivot_row];
            p[max_pivot_row] = tmp;
            for(k = 0; k < n; k++) {
                double *tmp_pointer1 = a + (i * n + k);
                double *tmp_pointer2 = a + (max_pivot_row * n + k);
                tmp = *tmp_pointer1;
                *tmp_pointer1 = *tmp_pointer2;
                *tmp_pointer2 = tmp;
            }
        }
        for(k = i + 1; k < n; k++) {
            *(a + (k * n + i)) /= *(a + (i * n + i));
            for(j = i + 1; j < n; j++) {
                *(a + (k * n + j)) -=  *(a + (k * n + i)) * *(a + (i * n + j));
            }
        }
    }
}

void lup_solve(double *a, int *p, double *b, double *x, int n) {
    for (int i = 0; i < n; i++) {
        x[i] = b[p[i]];
        for (int j = 0; j < i; j++)
            x[i] -= *(a + (i * n + j)) * x[j];
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++)
            x[i] -= *(a + (i * n + j)) * x[j];
        x[i] = x[i] / *(a + (i * n + i));
    }
}

void gaussian_elimination(double *a, double *b, double *x, int n) {

    int column, row, diagonal, max_pivot_row, j;
    double max_pivot, tmp;
    for(diagonal = 0; diagonal < n; diagonal++) {
        max_pivot_row = diagonal;
        max_pivot = *(a + (diagonal * n + diagonal));    // i,ith element of the matrix
        for(row = diagonal + 1; row < n; row++) {
            tmp = fabs(*(a + (row * n + diagonal)));
            if(tmp > max_pivot) {
                max_pivot_row = row;
                max_pivot = tmp;
            }
        }

        if(diagonal != max_pivot_row) {
            for(int k = 0; k < n; k++) {
                double *tmp_pointer1 = a + (diagonal * n + k);
                double *tmp_pointer2 = a + (max_pivot_row * n + k);
                tmp = *tmp_pointer1;
                *tmp_pointer1 = *tmp_pointer2;
                *tmp_pointer2 = tmp;
            }
            tmp = b[diagonal];
            b[diagonal] = b[max_pivot_row];
            b[max_pivot_row] = tmp;
        }

        for(row = diagonal + 1; row < n; row++) {
            tmp = *(a + (row * n + diagonal)) / *(a + (diagonal * n + diagonal));
            for(column = diagonal + 1; column < n; column++) {
                *(a + (row * n + column)) -= tmp * *(a + (diagonal * n + column));
            }
            *(a + (row * n + diagonal)) = 0;
            b[row] -= tmp * b[diagonal];
        }
    }

    for(row = n - 1; row >= 0; row--) {
        tmp = b[row];
        for(j = n - 1; j > row; j--) {
            tmp -= x[j] * *(a + (row * n + j));
        }
        x[row] = tmp / *(a + (row * n + row));
    }

}

long clctd(struct timeval start, struct timeval end)
{
    return (long)((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
}

int main(void) {

    static double a[N * N], a_prim[N * N];
    static double b[N], b_prim[N];
    static double x[N], x_prim[N];
    static int P[N];
    int s;
    struct timeval start_r, end_r;

    FILE *a_file = fopen("A_matrix_1.txt", "r");
    FILE *b_file = fopen("B_matrix_1.txt", "r");
    if(a_file == NULL || b_file == NULL) {
        fprintf(stderr, "Error while opening files.\n");
        exit(1);
    }

    srand(time(NULL));

    // printf("Matrix A:\n");
    for(int i = 0; i < N * N; i++) {
#ifndef RND
        fscanf(a_file ,"%lf", &a[i]);
#else
        a[i] = (double)random()/RAND_MAX * 5;
#endif
        a_prim[i] = a[i];
        // printf("%g ", a[i]);
        // if((i+1) % N == 0 && i) {
        //     printf("\n");
        // }
    }
    // printf("\n");

    // printf("Matrix B:\n");
    for(int i = 0; i < N; i++) {
#ifndef RND
        fscanf(b_file ,"%lf", &b[i]);
#else
        b[i] = (double)random()/RAND_MAX * 5;
#endif
        b_prim[i] = b[i];
        x[i] = 0.0;
        x_prim[i] = 0.0;
        P[i] = 0;
        // printf("%g\n", b[i]);
    }
    // printf("\n");

    gettimeofday(&start_r, NULL);
    gaussian_elimination(a, b, x, N);
    gettimeofday(&end_r, NULL);

    printf("Gaussian elimination with partial pivoting: %ld microseconds\n", clctd(start_r, end_r));

    gettimeofday(&start_r, NULL);
    lup_doolittle(a_prim, P, N);
    lup_solve(a_prim, P, b_prim, x_prim, N);
    gettimeofday(&end_r, NULL);

    printf("LU decomposition and solve (Doolittle): %ld microseconds\n\n", clctd(start_r, end_r));

    // printf("Gaussian elimination with partial pivoting result:\n");
    // for(int i = 0; i < N; i++) {
    //     printf("%.15g\n", x[i]);
    // }

    // printf("\nLU result:\n");
    // for(int i = 0; i < N; i++) {
    //     printf("%.15g\n", x_prim[i]);
    // }
    fclose(a_file);
    fclose(b_file);

    return 0;

}