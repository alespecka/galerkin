/* 
 * File:   GaussQuadTri.h
 * Author: ales
 *
 * Created on April 29, 2015, 1:51 PM
 */

#ifndef GAUSSQUADTRI_H
#define	GAUSSQUADTRI_H

#include <string>

class GaussQuadTri {
public:
    static const int nFaces = 3;
    int nWeights;
    double *weights;
    double **coords;
    
    GaussQuadTri(const std::string &fileName);
//    GaussQuadTri(const GaussQuadTri& orig);
    virtual ~GaussQuadTri();
    void getCoords(double &x, double &y, int indx, double (&nodeX)[nFaces], double (&nodeY)[nFaces]);
};

#endif	/* GAUSSQUADTRI_H */

