#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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
    static int P[3] = {0};
    static double X[3] = {0};

    printf("Matrix A: \n");
    print_matrix(A, 3);
    printf("\n");

    lup_doolittle(A, P, 3);

    printf("LUP decomposition (Doolittle) of A: \n");
    print_matrix(A, 3);
    printf("\n");

    printf("Columns in which the permutation matrix P has 1s for A: \n");
    for(int i = 0; i < 3; i++) {
        printf("%d ", P[i]);
    }
    printf("\n\n");

    lup_solve(A, P, B, X, 3);

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