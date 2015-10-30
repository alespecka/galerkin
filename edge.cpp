/*
 * edge.cpp
 *
 *  Created on: Nov 15, 2014
 *      Author: ales
 */

#include "edge.h"
#include <cmath>
#include <iostream>

//Edge::Edge(double ax, double ay, double bx, double by) {
void Edge::set(double ax, double ay, double bx, double by, double centreX,
               double centreY, double neighIndex) {
    normalX = by - ay;
    normalY = ax - bx;
    length = std::sqrt(normalX * normalX + normalY * normalY);

    normalX = normalX / length;
    normalY = normalY / length;
    
//    vectX = (ax + bx)/2 - centreX;
//    vectY = (ay + by)/2 - centreY;
    
//    neighEdge = NULL;
    boundaryType = neighIndex;
}
