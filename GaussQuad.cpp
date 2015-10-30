/* 
 * File:   GaussQuad.cpp
 * Author: ales
 * 
 * Created on April 28, 2015, 5:40 PM
 */

#include "GaussQuad.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>

GaussQuad::GaussQuad(const std::string &fileName) {
    std::ifstream file(fileName);
    if (!file) {
        std::stringstream what;
        what << "cannot open file " << fileName;
        throw std::runtime_error(what.str());      
    }
    
    file >> nPoints;
    weights = new double[nPoints];
    points = new double[nPoints];
    for (int i = 0; i < nPoints; ++i) {
        if (!(file >> points[i])) {
            file.close();
            delete this;
            
            std::stringstream what;
            what << "file " << fileName << " has a wrong format";
            throw std::runtime_error(what.str());
        }
    }
    for (int i = 0; i < nPoints; ++i) {
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

//GaussQuad::GaussQuad(const GaussQuad& orig) {
//}

GaussQuad::~GaussQuad() {
    delete weights;
    delete points;
}

double GaussQuad::getCoord(int indx, double pointA, double pointB) {
    return (pointB - pointA) / 2 * points[indx] + (pointA + pointB) / 2;
}

