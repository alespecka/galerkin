/* 
 * File:   element.h
 * Author: ales
 *
 * Created on 13. listopadu 2014, 13:43
 */

#pragma once

#define LAGRANGE_BASIS
// #define TAYLOR_BASIS

#include "edge.h"
#include "GaussQuad.h"
#include "GaussQuadTri.h"
#include "sparseMatrix.h"
#include "NumericalParameters.h"
#include "PhysicalParameters.h"

#include <omp.h>
#include <vector>
#include <armadillo>

class Element {
public:
    static const int nEqns = PhysicalParameters::nEqns;  // number of equations
    static const int nFaces = 3;    
    const int nBasis;
    
    typedef std::vector<Element *> ElemVector;
    typedef std::vector<double> DVector;    
    typedef arma::vec::fixed<nEqns> Vector4;
    typedef arma::vec::fixed<2> Vector2;
    typedef arma::mat::fixed<nEqns,2> Matrix42;    
    
    static bool isInitialised;
    static PhysicalParameters* param;
    static NumericalParameters* numParam;
    static void initialiseClass(PhysicalParameters* _physParam, NumericalParameters* _numParam);
    
    #ifdef _OPENMP
    static omp_lock_t* lock;
    static void initialiseClass(PhysicalParameters* _physParam, NumericalParameters* _numParam,
                                omp_lock_t* _lock);
    #endif /* _OPENMP */
    
    Element(int index, double* points, int* triangles, int* neighbours, int nBasis);
//    void setNeighbours(elemVector &mesh);
//    void findNeighFaces(ElemVector &mesh);
//    void print();
    
    double *basisConst;
    double **lagrangeCoeffs;
  
    double area;
    
    double centreX;
    double centreY;
    
    double nodeX[nFaces];
    double nodeY[nFaces];
    
    int neighbour[nFaces];
private:
    int index;
    double indiameter;
    int neighFace[nFaces];
    Edge edge[nFaces];

private:
    static double calculateArea(const double (&nodeX)[nFaces], const double (&nodeY)[nFaces]);
    static void calculateCentre(double &centreX, double &centreY, const double (&nodeX)[nFaces],
                         const double (&nodeY)[nFaces]);
    static double calculateIndiameter(const Edge (&edge)[nFaces], const double area);
    
//    double pressure;
//    double speedOfSound;
//    double derivX[dim];
//    double derivY[dim];
    /*
     * methods to be excluded from this class !!!
     */
    static void matrixTimesVector(double *x, double **A, double *b, int m, int n);
    static int devideElement(double ***pCoords, int n);
public:
//    double w[dim*nBasis];
//    double neww[dim*nBasis];

    double CFL(double *w);
    void combine(double *combinedW, double x, double y, double *w);
    void getW(arma::vec &w, double x, double y, const arma::mat &coeffs);
    void getGradW(double (&gradW)[2][nEqns], double x, double y, double *coeffs);
    void getGradW(arma::mat::fixed<nEqns,2> &gradW, double x, double y,
                      const arma::mat &coeffs);
    double basis(int indx, double x, double y);
    double derBasisX(int indx, double x, double y);
    double derBasisY(int indx, double x, double y);
    void lagrangianCoefficients();
    
    void artificialDamping(DVector &W, GaussQuadTri &quad, DVector &centreX, DVector &centreY);
    void massMatrix(double **M, GaussQuadTri &quad);
    arma::mat computeMassMatrix(GaussQuadTri &quad);
    void volumeIntegral(double *integral, DVector &W, GaussQuadTri &quad);
    arma::mat& volumeIntegral(arma::mat &integral, arma::mat &W, GaussQuadTri &quad);
    void lineIntegral(double *integral, int faceNumber, DVector &W, ElemVector &mesh, GaussQuad &quad);
    arma::mat& internalLineIntegral(arma::mat &integral, int faceNumber, arma::mat &W1,
                                    arma::mat &W2, ElemVector &mesh, GaussQuad &quad);
    arma::mat& boundaryLineIntegral(arma::mat &integral, int faceNumber, arma::mat &W1,
                                    ElemVector &mesh, GaussQuad &quad);
    void residualVector(DVector &RW, DVector &W, ElemVector &mesh, GaussQuad &quad, GaussQuadTri &quadTri);
    void residualVector(arma::mat &RW, arma::cube &W, ElemVector &mesh, GaussQuad &quad, GaussQuadTri &quadTri);
    void jacobi(SparseMatrix &jacobi, DVector &RW, DVector &W, ElemVector &mesh, double timeStep,
                GaussQuad &quad, GaussQuadTri &quadTri);
    void jacobi(SparseMatrix &jacobi, arma::mat &RW, arma::cube &W, ElemVector &mesh,
                   double timeStep, GaussQuad &quad, GaussQuadTri &quadTri);
    void insertJacobiBlock(SparseMatrix &jacobi, arma::mat &jacobiBlock, int rowNum, int colNum);
    
    void viscousNumericalFlux(double (&flux)[nEqns], double (&w1)[nEqns], double (&gradW1)[2][nEqns],
                                   double (&w2)[nEqns], double (&gradW2)[2][nEqns], const Edge& edge);
    void laxFriedrichs(double (&flux)[nEqns], const double (&w1)[nEqns], const double (&w2)[nEqns],
                       const Element& elem1, const Element& elem2, const Edge& edge);
    void vanLeer(double (&flux)[nEqns], double *w1, double *w2, const Edge& edge);
    arma::vec& vanLeer(arma::vec &flux, const arma::vec &w1, const arma::vec &w2, const Edge& edge);
    void AUSM(double (&flux)[nEqns], double *w1, double *w2, const Edge& edge);
    void boundary(double (&flux)[nEqns], double *w, const Edge& edge);
    void boundary(double (&flux)[nEqns], double (&w)[nEqns], double (&gradW)[2][nEqns], const Edge& edge);    
    arma::vec& boundary(Vector4 &flux, const Vector4 &w, const Edge& edge);
    arma::vec& viscousWall(Vector4 &flux, const Edge &edge, const Vector4 &w,
                           const Matrix42 &gradW);
    void physicalFluxes(double (&xflux)[nEqns], double (&yflux)[nEqns], double (&w)[nEqns]);    
    void physicalFluxes(double (&xflux)[nEqns], double (&yflux)[nEqns], double (&w)[nEqns], double p);
    void physicalFluxes(arma::vec &fluxX, arma::vec &fluxY, arma::vec &w);
    void physicalFluxes(arma::vec &fluxX, arma::vec &fluxY, arma::vec &w, double p);
    
    void stressTensor(arma::mat::fixed<2,2> &stress, const Vector2 &gradU, const Vector2 &gradV);
    void derivatives(const Vector4 &w, const Matrix42 &gradW, double p,
                          Vector2 &gradU, Vector2 &gradV, Vector2 &gradPOverRho);
    void viscousFluxes(Vector4 &fluxX, Vector4 &fluxY, const Vector4 &w,
                            const Matrix42 &gradW);
    arma::vec& viscousNumericalFlux(Vector4 &flux, const Edge &edge, const Vector4 &w1,
                               const Vector4 &w2, const Matrix42 &gradW1, const Matrix42 &gradW2);   
    
    void grad(double (&gradRho)[2], double (&gradU)[2], double (&gradV)[2], double (&gradPOverRho)[2],
              double (&w)[nEqns], double (&gradW)[2][nEqns]);
    void stressTensor(double (&tensor)[2][2], double (&gradU)[2], double (&gradV)[2]);
    void viscousFluxes(double (&fluxX)[nEqns], double (&fluxY)[nEqns],
                       double (&w)[nEqns], double (&gradW)[2][nEqns]);
    void computeWallW(Vector4 &wallW, const Vector4 &w);
//    void negativeDensityCheck(dVector &W);
};
