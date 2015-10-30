/*
 * edge.h
 *
 *  Created on: Nov 15, 2014
 *      Author: ales
 */

#pragma once

struct Edge {
//    static const int dim = 4;
    
//    Edge* neighEdge;  // je to potreba ???
    double length;
    double normalX, normalY;
//    double vectX, vectY;
//    double w[dim];
    
    /* boundaryType is less than 0 if the edge is part of the boundary
     * in which case it has a value of either -1 (WALL), -2 (INLET) or -3 (OUTLET) */
    double boundaryType;

//	Edge(double ax, double ay, double bx, double by);
    void set(double ax, double ay, double bx, double by, double centreX, double centreY, double neighIndex);
};
