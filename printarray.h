/*
 * printarray.h
 *
 *  Created on: Nov 16, 2014
 *      Author: ales
 */

#pragma once
#include <iostream>
#include <fstream>

template<typename Type>
void printColumn(Type* array, int length) {
    for (int i = 0; i < length; ++i) {
//        std::cout << array[i] << " ";
        printf("%.4e\n", array[i]);
    }
    printf("\n");
}

template<typename Type>
void printArray(Type* array, int length) {
    for (int i = 0; i < length; ++i) {
        std::cout << array[i] << ' ';
//        printf("%.4e ", array[i]);
    }
    std::cout << std::endl << std::endl;
}

template<typename Type>
void printArray(Type* array, int nRows, int nColumns) {
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nColumns; ++j) {
//            std::cout << array[i * nColumns + j] << " ";
            printf("%.4e   ", array[i*nColumns + j]);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template<typename Type>
void printColMajorMatrix(Type* array, int nRows, int nColumns) {    
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nColumns; ++j) {
//            std::cout << array[i + j*nRows] << " ";
            printf("%.4e   ", array[i + j*nRows]);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template<typename Type>
void printJaggedArray(Type **array, int nRows, int nColumns) {
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nColumns; ++j) {
            std::cout << array[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template<typename Type>
void printColMajorJaggedArray(Type **array, int nRows, int nColumns) {
    for (int j = 0; j < nColumns; ++j) {
        for (int i = 0; i < nRows; ++i) {
            std::cout << array[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template<typename Type>
int readArray(Type *&array, const std::string &fileName) {
    std::ifstream file(fileName);
    
    if (!file) {
        std::stringstream what;
        what << "cannot open file " << fileName;
        throw std::runtime_error(what.str());
    }
    
    Type temp;
    int nValues = 0;
    while (file >> temp) {
        ++nValues;
    }
    
    file.close();
    file.clear();
    file.open(fileName.c_str());
    
    array = new Type[nValues]; 
    
    for (int i = 0; i < nValues; ++i) {
        file >> array[i];
    }
    file.close();
    
    return nValues;
}

template<typename Type>
int readArray(Type* array, int len, const std::string &fileName) {
    std::ifstream file(fileName);
    
    if (!file) {
        std::cout << "Cannot open file " << fileName << "/n";
        return -1;
    }
    
    Type temp;
    int nValues = 0;
    while (file >> temp) {
        ++nValues;
    }
    
    file.close();
    file.clear();
    file.open(fileName.c_str());
    
    for (int i = 0; i < nValues && i < len; ++i) {
        file >> array[i];
    }
    file.close();
    
    if (nValues <= len) {
        return nValues;
    }
    else {
        std::cout << "Warning: reading file " << fileName << " was not finished.\n";
        return len;
    }
} 
