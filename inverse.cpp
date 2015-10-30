#include <cmath>
#include <iostream>

// Method to carry out the partial-pivoting Gaussian
// elimination. Here index[] stores pivoting order.

void gaussian(double *A, int *index, int n){
    //int n = index.length;
    double c[n];

    // Initialize the index
    for (int i=0; i<n; ++i) index[i] = i;

    // Find the rescaling factors, one from each row
    for (int i=0; i<n; ++i) {
        double c1 = 0;
            for (int j=0; j<n; ++j) {
//                double c0 = std::abs(A[i][j]);
                double c0 = std::abs(A[i*n+j]);
                if (c0 > c1) c1 = c0;
            }
        c[i] = c1;
    }

    // Search the pivoting element from each column
    int k = 0;
    for (int j=0; j<n-1; ++j) {
        double pi1 = 0;
        for (int i=j; i<n; ++i) {
//            double pi0 = std::abs(A[index[i]][j]);
            double pi0 = std::abs(A[index[i]*n+j]);
            pi0 /= c[index[i]];
            if (pi0 > pi1) {
                pi1 = pi0;
                k = i;
            }
        }

        // Interchange rows according to the pivoting order
        int itmp = index[j];
        index[j] = index[k];
        index[k] = itmp;
        for (int i=j+1; i<n; ++i) {
//            double pj = A[index[i]][j]/A[index[j]][j];
            double pj = A[index[i]*n+j]/A[index[j]*n+j];

            // Record pivoting ratios below the diagonal
//        A[index[i]][j] = pj;
            A[index[i]*n+j] = pj;

            // Modify other elements accordingly
            for (int l=j+1; l<n; ++l)
//            A[index[i]][l] -= pj*A[index[j]][l];
            A[index[i]*n+l] -= pj*A[index[j]*n+l];
        }
    }
}

void invert(double **invA, double **A, int n) { // inverze matice
//    double **A0 = new double*[n];
//    for(int i = 0; i < n; i++){
//        A0[i] = new double[n];
//        for(int j = 0; j < n; j++){
//            A0[i][j] = A[i][j];
//        }
//    }
    
    double A0[n*n];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A0[i*n+j] = A[i][j];
        }
    }
//    std::fill(invA, invA + n * n, 0);  // perhaps not necessary
//    double invA[n][n];
    double b[n][n];
    int index[n];
    for (int i=0; i<n; ++i) {
        std::fill(&b[i][0], &b[i][0] + n, 0);
        b[i][i] = 1;
    }

    // Transform the matrix into an upper triangle
    gaussian(A0, index, n);

    // Update the matrix b[i][j] with the ratios stored
    for (int i=0; i<n-1; ++i) {
        for (int j=i+1; j<n; ++j) {
            for (int k=0; k<n; ++k) {
//                b[index[j]][k] -= A0[index[j]][i] * b[index[i]][k];
                b[index[j]][k] -= A0[index[j]*n+i] * b[index[i]][k];
            }
        }
    }

    // Perform backward substitutions
    for (int i=0; i<n; ++i) {
//        invA[n-1][i] = b[index[n-1]][i]/A0[index[n-1]][n-1];
        invA[n-1][i] = b[index[n-1]][i]/A0[index[n-1]*n+(n-1)];
        for (int j=n-2; j>=0; --j) {
//            invA[j][i] = b[index[j]][i];
            invA[j][i] = b[index[j]][i];
            for (int k=j+1; k<n; ++k) {
//                invA[j][i] -= A0[index[j]][k]*invA[k][i];
                invA[j][i] -= A0[index[j]*n+k]*invA[k][i];
            }
//            invA[j][i] /= A0[index[j]][j];
            invA[j][i] /= A0[index[j]*n+j];
        }
    }
    
//    return invA;
}

