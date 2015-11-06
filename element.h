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
    static std::shared_ptr<PhysicalParameters> param;
    static std::shared_ptr<NumericalParameters> numParam;
    static std::shared_ptr<GaussQuad> quad;
    static std::shared_ptr<GaussQuadTri> quadTri;
    static void initialiseClass(std::shared_ptr<PhysicalParameters> _physParam,
                                std::shared_ptr<NumericalParameters> _numParam,
                                std::shared_ptr<GaussQuad> _quad,
                                std::shared_ptr<GaussQuadTri> _quadTri);
    
    #ifdef _OPENMP
        static omp_lock_t* lock;
        static void initialiseClass(std::shared_ptr<PhysicalParameters> _physParam,
                                    std::shared_ptr<NumericalParameters> _numParam,
                                    std::shared_ptr<GaussQuad> _quad,
                                    std::shared_ptr<GaussQuadTri> _quadTri,
                                    omp_lock_t* _lock);
    #endif /* _OPENMP */
    
    Element(int index, double* points, int* triangles, int* neighbours, int nBasis);
//    void setNeighbours(elemVector &mesh);
//    void findNeighFaces(ElemVector &mesh);
//    void print();  
    
    double amountOfViscosity;
    
//    double *basisConst;   
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
    static int devideElement(double ***pCoords, int n);
public:
//    double w[dim*nBasis];
//    double neww[dim*nBasis];

    double CFL(double *w);
    void getW(arma::vec &w, double x, double y, const arma::mat &coeffs);
    void getGradW(double (&gradW)[2][nEqns], double x, double y, double *coeffs);
    void getGradW(arma::mat::fixed<nEqns,2> &gradW, double x, double y,
                      const arma::mat &coeffs);
    double basis(int indx, double x, double y);
    double derBasisX(int indx, double x, double y);
    double derBasisY(int indx, double x, double y);
    void lagrangianCoefficients();
    
    void shockIndicator(arma::mat &W);    
    arma::mat computeMassMatrix();    
    arma::mat& volumeIntegral(arma::mat &integral, arma::mat &W);
    arma::mat& internalLineIntegral(arma::mat &integral, int faceNumber, arma::mat &W1,
                                    arma::mat &W2, ElemVector &mesh);
    arma::mat& boundaryLineIntegral(arma::mat &integral, int faceNumber, arma::mat &W1,
                                    ElemVector &mesh);
    void residualVector(arma::mat &RW, arma::cube &W, ElemVector &mesh);
    void jacobi(SparseMatrix &jacobi, arma::mat &RW, arma::cube &W, ElemVector &mesh,
                   double timeStep);
    void insertJacobiBlock(SparseMatrix &jacobi, arma::mat &jacobiBlock, int rowNum, int colNum);
    arma::vec& laxFriedrichs(Vector4& flux, const Vector4& w1, const Vector4& w2, const Edge& edge);
    arma::vec& vanLeer(arma::vec &flux, const arma::vec &w1, const arma::vec &w2, const Edge& edge);    
    arma::vec& boundary(Vector4 &flux, const Vector4 &w, const Edge& edge);
    arma::vec& viscousWall(Vector4 &flux, const Edge &edge, const Vector4 &w,
                           const Matrix42 &gradW);
    void cartesianPhysicalFluxes(Matrix42& flux, const Vector4& w);
    void cartesianPhysicalFluxes(Matrix42& flux, const Vector4& w, double p);
    void normalPhysicalFlux(Vector4& flux, const Vector4& w, const Edge& edge, double p);
    
    void computeWallW(Vector4 &wallW, const Vector4 &w);
    void stressTensor(arma::mat::fixed<2,2> &stress, const Vector2 &gradU, const Vector2 &gradV);
    void derivatives(const Vector4 &w, const Matrix42 &gradW, double p,
                          Vector2 &gradU, Vector2 &gradV, Vector2 &gradPOverRho);
    void viscousFluxes(Matrix42& cartesianFlux, const Vector4& w, const Matrix42& gradW);    
    arma::vec& viscousNumericalFlux(Vector4 &flux, const Edge &edge, const Vector4 &w1,
        const Vector4 &w2, const Matrix42 &gradW1, const Matrix42 &gradW2, double amountOfViscosity);       
//    void negativeDensityCheck(dVector &W);
};
