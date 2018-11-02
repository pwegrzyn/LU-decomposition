#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void lup_cholesky(double *a, int n) {
    int i, j, k;
    double sum = 0.0;
    double *l = (double*)calloc(n * n, sizeof(double));
    for(i = 0; i < n; i++) {
        for(j = 0; j < (i + 1); j++) {
            sum = 0.0;
            for(k = 0; k < j; k++) {
                sum += *(l + (i * n + k)) * *(l + (j * n + k));
            }
            if (i == j) {
                *(l + (i * n + j)) = sqrt(*(a + (i * n + i)) - sum);
            } else {
                *(l + (i * n + j)) = 1.0 / *(l + (j * n + j)) * (*(a + (i * n + j)) - sum);
            }
        }
    }
    for(i = 0; i < n * n; i++) {
            a[i] = l[i];
    }
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            if(i <= j) {
                *(a + (i * n + j)) = *(l + (j * n + i));
            }
        }
    }
    free(l);
}

void lup_solve(double *a, double *b, double *x, int n) {
    for (int i = 0; i < n; i++) {
        x[i] = b[i];
        for (int j = 0; j < i; j++)
            x[i] -= *(a + (i * n + j)) * x[j];
        x[i] = x[i] / *(a + (i * n + i));
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++)
            x[i] -= *(a + (i * n + j)) * x[j];
        x[i] = x[i] / *(a + (i * n + i));
    }
}

void print_matrix(double *a, int n) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%f ", *(a + (i * n + j)));
        }
        printf("\n");
    }
}

int main(void) {

    static double A[9] = {4, 12, -16, 12, 37, -43, -16, -43, 98};
    static double B[3] = {3.3, 4.4, 2.5};
    static double X[3] = {0};

    printf("Matrix A: \n");
    print_matrix(A, 3);
    printf("\n");

    lup_cholesky(A, 3);

    printf("LUP decomposition (Cholesky) of A: \n");
    print_matrix(A, 3);
    printf("\n");

    lup_solve(A, B, X, 3);

    printf("B vector: \n");
    for(int i = 0; i < 3; i++) {
        printf("%f ", B[i]);
    }
    printf("\n\n");
    
    printf("Solution for Ax=B: \n");
    for(int i = 0; i < 3; i++) {
        printf("%f ", X[i]);
    }
    printf("\n\n");

    return 0;

}