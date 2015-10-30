/* 
 * File:   sparseMatrix.cpp
 * Author: ales
 * 
 * Created on April 30, 2015, 5:30 PM
 */

#include "sparseMatrix.h"

#include <assert.h>
#include <vector>
#include <iostream>
#include <stdio.h>

SparseMatrix::SparseMatrix(int size) {
    this->size = size;
    position = 0;
    
    rowIdx.reserve(size);
    colIdx.reserve(size);
    values.reserve(size);
}

//sparseMatrix::sparseMatrix(int *rowIdx, int *colIdx, double *values, int size, int position, int dim) {
//    this->size = size;
//    this->position = position;
//    this->dim = dim;
//    this->rowIdx = rowIdx;
//    this->colIdx = colIdx;
//    this->values = values;
//}

void SparseMatrix::insert(int rowIndex, int colIndex, double value) {
//    assert(position < size);
    rowIdx[position] = rowIndex;
    colIdx[position] = colIndex;
    values[position] = value;
    ++position;
}

void SparseMatrix::print(int start, int end) {
    assert(start >= 0);
    for (int i = start; i < end && i <= position && i < size; ++i) {
        //std::cout << "(" << rowIdx[i] << " " << colIdx[i] << ") " << values[i] << "\n";
        printf("(%d, %d) %e\n", rowIdx[i], colIdx[i], values[i]);
    }
    std::cout << std::endl;
}

//SparseMatrix::SparseMatrix(const SparseMatrix& orig) {
//    delete rowIdx;
//}

SparseMatrix::~SparseMatrix() {
}

