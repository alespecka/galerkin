/* 
 * File:   GaussQuadTri.cpp
 * Author: ales
 * 
 * Created on April 29, 2015, 1:51 PM
 */

#include "GaussQuadTri.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>

GaussQuadTri::GaussQuadTri(const std::string &fileName) {
    std::ifstream file(fileName);
    
    if (!file) {
        std::stringstream what;
        what << "cannot open file " << fileName;
        throw std::runtime_error(what.str()); 
    }
    
    file >> nWeights;

//    double coords[nWeights][nFaces];
//    double weights[nWeights];
    weights = new double[nWeights];
    coords = new double*[nWeights];
    for (int i = 0; i < nWeights; ++i) {
        coords[i] = new double[nFaces];
    }
//    *pWeights = weights;
//    *pCoords = coords;
    
    for (int i = 0; i < nWeights; ++i) {
        for (int j = 0; j < nFaces; ++j) {
            if (!(file >> coords[i][j])) {
                file.close();
                delete this;
                
                std::stringstream what;
                what << "file " << fileName << " has a wrong format";
                throw std::runtime_error(what.str());
            }
        }
    }
    for (int i = 0; i < nWeights; ++i) {
        if (!(file >> weights[i])) {
            file.close();
            delete this;
            
            std::stringstream what;
            what << "file " << fileName << " has a wrong format";
            throw std::runtime_error(what.str());
        }
    }
    
    file.close();
}

//GaussQuadTri::GaussQuadTri(const GaussQuadTri& orig) {
//}

GaussQuadTri::~GaussQuadTri() {
    delete weights;
    // CHYBA !!!
//    for (int i = 0; i < nWeights; ++i) {
//        delete coords[i];
//    }
//    delete coords;
}

void GaussQuadTri::getCoords(double &x, double &y, int indx, double (&nodeX)[nFaces], double (&nodeY)[nFaces]) {
    x = 0;
    y = 0;
    for (int l = 0; l < nFaces; ++l) {
        x += coords[indx][l] * nodeX[l];
        y += coords[indx][l] * nodeY[l];
    }
}

