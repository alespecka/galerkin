#include <cmath>

void lsolve(double *x, double **A, double *b, int N) {
// Gaussian elimination with partial pivoting
    
    double b0[N];
    double *A0[N];
    for (int i = 0; i < N; ++i) {
        b0[i] = b[i];
        A0[i] = A[i];
    }
    
    for (int p = 0; p < N; p++) {
        // find pivot row and swap
        int max = p;
        for (int i = p + 1; i < N; i++) {
            if (std::abs(A0[i][p]) > std::abs(A0[max][p])) {
                max = i;
            }
        }

        double *temp = A0[p];
        A0[p] = A0[max];
        A0[max] = temp;
        
        double t = b0[p];
        b0[p] = b0[max];
        b0[max] = t;

        // pivot within A and b
        for (int i = p + 1; i < N; i++) {
            double alpha = A0[i][p] / A0[p][p];
            b0[i] -= alpha * b0[p];
            for (int j = p; j < N; j++) {
                A0[i][j] -= alpha * A0[p][j];
            }
        }
    }

    // back substitution
    for (int i = N - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < N; j++) {
            sum += A0[i][j] * x[j];
        }
        x[i] = (b0[i] - sum) / A0[i][i];
    }
}

