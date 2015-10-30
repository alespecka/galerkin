/* 
 * File:   sparseMatrix.h
 * Author: ales
 *
 * Created on April 30, 2015, 5:30 PM
 */

#ifndef SPARSEMATRIX_H
#define	SPARSEMATRIX_H

#include <vector>
typedef std::vector<int> intVector;
typedef std::vector<double> dVector;

class SparseMatrix {
public:
    int position;
    int size;
    intVector rowIdx;
    intVector colIdx;
    dVector values;
    
    SparseMatrix(int size);
//    SparseMatrix(int *rowIdx, int *colIdx, double *values, int size, int position, int dim);
//    SparseMatrix(const SparseMatrix& orig);
    virtual ~SparseMatrix();
    void insert(int rowIdx, int colIdx, double value);
    void print(int start, int end);
private:

};

#endif	/* SPARSEMATRIX_H */

