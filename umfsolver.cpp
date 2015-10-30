//#include <cstdlib>
//#include <stdio.h>
#include "umfpack.h"

void quickSort(int* a, int* b, double* c, int left, int right) {
      int i = left, j = right;
      int tmp;
      double tmp_d;
      int pa = a[(left + right)/2];
      int pb = b[(left + right)/2];

      /* partition */
      while (i <= j) {
            while (a[i] < pa or (a[i] == pa and b[i] < pb))
                  i++;
            while (a[j] > pa or (a[j] == pa and b[j] > pb))
                  j--;
            if (i <= j) {
                  tmp = a[i];
                  a[i] = a[j];
                  a[j] = tmp;
                  tmp = b[i];
                  b[i] = b[j];
                  b[j] = tmp;
                  tmp_d = c[i];
                  c[i] = c[j];
                  c[j] = tmp_d;
                  i++;
                  j--;
            }
      };

      /* recursion */
      if (left < j)
            quickSort(a, b, c, left, j);
      if (i < right)
            quickSort(a, b, c, i, right);
}


void UMFPACK(int n, int *Ap, int *Ai, double *Ax, double *b, double *x)
{
	void *Symbolic, *Numeric;
	int i;

	/* symbolic analysis */
	umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);

	/* LU factorization */
	umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
	umfpack_di_free_symbolic(&Symbolic);

	/* solve system */
	umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);
	umfpack_di_free_numeric(&Numeric);
}


void umfSolve(int nnz, int n, int* J, int* I, double* H, double* b, double* x) {
    int i;

    // serazeni
    quickSort(I,J,H,0,nnz-1);

    // ccs format ridke matice
    int *Ic = new int[n+1];
//    int Ic[n+1];
    Ic[0] = 0;
    int s = 1;
    i = 0;
    while(i <= nnz){
    	while(I[i] == I[i+1]){
            i++;
    	}
        i++;
        Ic[s] = i;
        s++;
    }
    Ic[n] = nnz;

    UMFPACK(n,Ic,J,H,b,x);
    
    delete Ic;
}

//int main(){
//
//    // example
//    // |1 0 2|   |x0|    |1|
//    // |0 4 3| * |x1| =  |2|
//    // |1 0 1|   |x2|    |3|
//    
//    int n = 3;
//    int nnz = 6;
//    int I[] = {0,0,1,1,2,2} ;
//    int J[] = {0,2,1,2,0,2} ;
//    double H[] = {1,2,4,3,1,1} ;
//    double b[] = {1,2,3} ;
//    double x[n];
//
//    solve(nnz,n,I,J,H,b,x);
//    
//    for(int i = 0; i < n; i++)
//        printf("x%i = %6.3f\n",i, x[i]);
//    
//    return 0;
//}
