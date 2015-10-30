/* 
 * File:   GaussQuad.h
 * Author: ales
 *
 * Created on April 28, 2015, 5:40 PM
 */

#ifndef GAUSSQUAD_H
#define	GAUSSQUAD_H

#include <string>

class GaussQuad {
public:
    GaussQuad(const std::string &fileName);
//    GaussQuad(const GaussQuad& orig);
    virtual ~GaussQuad();
double getCoord(int indx, double pointA, double pointB);
//private:
    int nPoints;
    double *weights;
    double *points;
};

#endif	/* GAUSSQUAD_H */

