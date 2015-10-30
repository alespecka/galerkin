/*
 * mesh.cpp
 *
 *  Created on: Nov 16, 2014
 *      Author: ales
 */

#include "mesh.h"
#include <iostream>
#include <list>
#include <cmath>

Element** Mesh::createMesh(double* points, int nTriangles, int* triangles, int* neighbours) {
    Element** mesh = new Element*[nTriangles];
    for (int i = 0; i < nTriangles; ++i) {
        mesh[i] = new Element(i, points, triangles, neighbours);
    }

    return mesh;
}


void Mesh::destructMesh(Element** mesh, int nTriangles) {
	for (int i = 0; i < nTriangles; ++i) {
		delete mesh[i];
	}
	delete[] mesh;
}

void Mesh::solvePDE(Element** mesh, int nTriangles) {
    
}


//Element* Mesh::createMesh(double *points, int nTriangle, int *triangles,
//                                int *neighbours) {
//    Element **mesh = new Element*[nTriangle];
//
//    mesh[0] = new Element(0, points, triangles);
//    for (int i = 1; i < nTriangle; ++i) {
//        mesh[i] = new Element(i, points, triangles, neighbours);
//        mesh[i-1]->next = mesh[i];
//    }
//
//    Element *triangle = mesh[0];
//    delete mesh;
//
//    return triangle;
//}
//
//void Mesh::destructMesh(Element *triangle, int nTriangles) {
//    Element *tri, *pom;
//    for (tri = triangle, pom = tri; tri != NULL; tri = pom) {
//        pom = pom->next;
//        delete tri;
//    }
//}
