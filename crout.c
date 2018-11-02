#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void lup_crout(double *a, int n) {
    int i, j, k;
    for(i = 0; i < n; i++) {
        for(j = i + 1; j < n; j++) {
            *(a + (i * n + j)) /= *(a + (i * n + i));
            for(k = i + 1; k < n; k++) {
                *(a + (k * n + j)) -= *(a + (i * n + j)) * *(a + (k * n + i));
            }
        }
    }
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

    lup_crout(A, 3);

    printf("LUP decomposition (Crout) of A: \n");
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