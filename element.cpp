#include "element.h"
#include "inverse.h"
#include "printarray.h"

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <ctime>
//#include <iostream>

#define BOUNDARY -1
#define WALL -1
#define INLET -2
#define OUTLET -3

//// row-major acces to matrices
//#define AT2(i, j)     (i)*n + (j)
//#define AT3(i, j, k)  (i)*n + (j) + (k)*m*n

// column-major access to matrices
#define AT2(i, j)     (i) + (j)*nBasis
#define AT3(i, j, k)  (i) + (j)*nBasis + (k)*nBasis*nEqns

// overloading AT macro
#define GET_MACRO(_1,_2,_3,NAME,...) NAME
#define AT(...) GET_MACRO(__VA_ARGS__, AT3, AT2)(__VA_ARGS__)


bool Element::isInitialised = false;

PhysicalParameters* Element::param = nullptr;
NumericalParameters* Element::numParam = nullptr;

void Element::initialiseClass(PhysicalParameters* _physParam, NumericalParameters* _numParam)
{
    param = _physParam;
    numParam = _numParam;
    isInitialised = true;
}

#ifdef _OPENMP
omp_lock_t* Element::lock = nullptr;

void Element::initialiseClass(PhysicalParameters* _physParam, NumericalParameters* _numParam,
                                omp_lock_t* _lock)
{
    param = _physParam;
    numParam = _numParam;
    lock = _lock;
    isInitialised = true;
}
#endif /* _OPENMP */

Element::Element(int index, double* points, int* triangles, int* neighbours, int nBasis) :
    nBasis(nBasis)
{
    if (!isInitialised) {
        throw std::runtime_error("class Element has not been initialised, use static method initialiseClass(...) first");
    }
    
    // double (&w0)[dim],
    this->index = index;

    for (int j = 0; j < nFaces; ++j) {
        neighbour[j] = neighbours[nFaces * index + j]; // set indexes of neighbours

        int pointIndex = triangles[nFaces * index + j];
        nodeX[j] = points[2 * pointIndex];
        nodeY[j] = points[2 * pointIndex + 1];
    }

    // std::fill(neighFace, neighFace + Element::nFaces, -1);

    calculateCentre(centreX, centreY, nodeX, nodeY);

    for (int j = 0; j < nFaces; ++j) {
        edge[j].set(nodeX[j], nodeY[j], nodeX[(j + 1) % 3], nodeY[(j + 1) % 3], centreX,
                centreY, neighbours[nFaces * index + j]);
    }

    area = calculateArea(nodeX, nodeY);
    indiameter = calculateIndiameter(edge, area);
    numParam->domainArea += area;

    for (int j = 0; j < nFaces; ++j) {
        if (neighbour[j] <= BOUNDARY) {
            ++numParam->nBoundaryEdges;
        }
    }

    // memory leaks   
    lagrangianCoefficients();

    //    pressure = (kapa-1) * (w[3] - 0.5*(w[1]*w[1]+w[2]*w[2])/w[0]);
    //    speedOfSound = std::sqrt(kapa * pressure / w[0]);
}

//void Element::findNeighFaces(ElemVector &mesh) {
//    for (int j = 0; j < nFaces; ++j) {
//        //        double faceIndex = -1;
//        for (int k = 0; k < nFaces; ++k) {
//            if (neighbour[j] >= 0 && index == mesh[neighbour[j]]->neighbour[k]) {
//                //                faceIndex = k;
//                edge[j].neighEdge = &mesh[neighbour[j]]->edge[k];
//                break;
//            }
//        }
//        //        neighFace[j] = faceIndex;
//    }
//}

double Element::calculateArea(const double (&nodeX)[nFaces], const double (&nodeY)[nFaces]) {
    double area = 0;
    for (int j = 1; j < nFaces; ++j) {
        area += nodeX[j - 1] * nodeY[j] - nodeX[j] * nodeY[j - 1];
    }
    area += nodeX[nFaces - 1] * nodeY[0] - nodeX[0] * nodeY[nFaces - 1];

    return 0.5 * std::abs(area);
    //    return 0.5 * std::abs((nodeX[1]-nodeX[0])*(nodeY[2]-nodeY[0]) - (nodeY[1]-nodeY[0])*(nodeX[2]-nodeX[0]));
}

void Element::calculateCentre(double &centreX, double &centreY, const double (&nodeX)[nFaces],
        const double (&nodeY)[nFaces]) {
    centreX = 0;
    centreY = 0;
    for (int j = 0; j < nFaces; ++j) {
        centreX += nodeX[j];
        centreY += nodeY[j];
    }
    centreX = centreX / nFaces;
    centreY = centreY / nFaces;
}

double Element::calculateIndiameter(const Edge(&edge)[nFaces], const double area) {
    double circumference = 0;
    for (int j = 0; j < nFaces; ++j) {
        circumference += edge[j].length;
    }
    return 4 * area / circumference;
}

double Element::CFL(double *w) {
    //    double w[dim];
    //    combine(w, centreX, centreY, longW, nBasis);
    double pressure = (param->kapa - 1) * (w[3] - 0.5 * (w[1] * w[1] + w[2] * w[2]) / w[0]);

//    if (pressure < numParam->tol) {
//        std::cerr << "pressure is too low" << std::endl;
//        pressure = numParam->tol;
//    }

    double speedOfSound = std::sqrt(param->kapa * pressure / w[0]);

    double lambdaX = std::abs(w[1] / w[0]) + speedOfSound;
    double lambdaY = std::abs(w[2] / w[0]) + speedOfSound;
    
    if (param->isFlowViscous) {
        return 1 / ((lambdaX + lambdaY) / indiameter + 2 / param->reynolds / (indiameter * indiameter));
    }
    else {
        return indiameter / (lambdaX + lambdaY);
    }
}

void Element::getW(arma::vec &w, double x, double y, const arma::mat &coeffs) {
    w.zeros();

    for (int i = 0; i < nBasis; ++i) {
        double basisValue = basis(i, x, y);
        // can be vectorised
        for (int d = 0; d < nEqns; ++d) {
            w(d) += basisValue * coeffs(i, d);
        }
    }
}

void Element::combine(double *w, double x, double y, double *coeffs) {
    std::fill(w, w + nEqns, 0);

    for (int i = 0; i < nBasis; ++i) {
        double bas = basis(i, x, y);
        for (int d = 0; d < nEqns; ++d) {
            w[d] += bas * coeffs[d * nBasis + i];
        }
    }
}

void Element::getGradW(Matrix42 &gradW, double x, double y,
                           const arma::mat &coeffs)
{
    gradW.zeros();
    
    for (int i = 0; i < nBasis; ++i) {
        double derBasX = derBasisX(i, x, y);
        double derBasY = derBasisY(i, x, y);
        for (int d = 0; d < nEqns; ++d) {
            gradW(d,0) += derBasX * coeffs(i,d);
            gradW(d,1) += derBasY * coeffs(i,d);
        }
    }    
}

void Element::getGradW(double (&gradW)[2][nEqns], double x, double y, double *coeffs) {
    std::fill(gradW[0], gradW[0] + nEqns, 0);
    std::fill(gradW[1], gradW[1] + nEqns, 0);
    
    for (int i = 0; i < nBasis; ++i) {
        double derBasX = derBasisX(i, x, y);
        double derBasY = derBasisY(i, x, y);
        for (int d = 0; d < nEqns; ++d) {
            gradW[0][d] += derBasX * coeffs[d * nBasis + i];
            gradW[1][d] += derBasY * coeffs[d * nBasis + i];
        }
    }    
}

#ifdef TAYLOR_BASIS
double Element::basis(int indx, double x, double y) {
    int m = numParam->powers[0][indx];
    int n = numParam->powers[1][indx];
    
    return std::pow(x - centreX, m) * std::pow(y - centreY, n); // - basisConst[indx];
}

double Element::derBasisX(int indx, double x, double y) {
    int m = numParam->powers[0][indx];
    
    if (m >= 1) {
        int n = numParam->powers[1][indx];
        return m * std::pow(x - centreX, m - 1) * std::pow(y - centreY, n);
    }
    else {
        return 0;
    }
}

double Element::derBasisY(int indx, double x, double y) {
    int n = numParam->powers[1][indx];

    if (n >= 1) {
        int m = numParam->powers[0][indx];
        return n * std::pow(x - centreX, m) * std::pow(y - centreY, n - 1);
    }
    else {
        return 0;
    }
}
#endif /* TAYLOR_BASIS */

#ifdef LAGRANGE_BASIS
double Element::basis(int indx, double x, double y) {
    double value = 0;
    int k = 0;
    for (int order = 0; order < numParam->orderOfOccuracy; ++order) {
        for (int i = 0; i < order + 1; ++i) {
            value += lagrangeCoeffs[k][indx] * std::pow(x, order - i) * std::pow(y, i);
            ++k;
        }
    }

    return value;
}

double Element::derBasisX(int indx, double x, double y) {
    double value = 0;
    int k = 0;
    for (int order = 0; order < numParam->orderOfOccuracy; ++order) {
        for (int i = 0; i < order; ++i) {
            value += (order - i) * lagrangeCoeffs[k][indx] * std::pow(x, order - i - 1) * std::pow(y, i);
            ++k;
        }
        ++k;
    }

    return value;
}

double Element::derBasisY(int indx, double x, double y) {
    double value = 0;
    int k = 0;
    for (int order = 0; order < numParam->orderOfOccuracy; ++order) {
        ++k;
        for (int i = 1; i < order + 1; ++i) {
            value += i * lagrangeCoeffs[k][indx] * std::pow(x, order - i) * std::pow(y, i - 1);
            ++k;
        }
    }

    return value;
}
#endif /* LAGRANGE_BASIS */

/*
 * calculate coefficients for the Lagrangian basis
 */
void Element::lagrangianCoefficients() {
    double **coords;
    int len = devideElement(&coords, numParam->orderOfOccuracy); // memory leak
    assert(len == nBasis);

    double x[len];
    double y[len];
    std::fill(x, x + len, 0);
    std::fill(y, y + len, 0);
    for (int j = 0; j < len; ++j) {
        for (int i = 0; i < nFaces; ++i) {
            x[j] += nodeX[i] * coords[i][j];
            y[j] += nodeY[i] * coords[i][j];
        }
    }

    double **Vand = new double*[len];
    lagrangeCoeffs = new double*[len]; // memory leak
    for (int i = 0; i < len; ++i) {
        Vand[i] = new double[len];
        lagrangeCoeffs[i] = new double[len];  // memory leak
    }

    int k = 0;
    for (int order = 0; order < numParam->orderOfOccuracy; ++order) {
        for (int i = 0; i < order + 1; ++i) {
            for (int j = 0; j < len; ++j) {
                Vand[j][k] = std::pow(x[j], order - i) * std::pow(y[j], i);
            }

            ++k;
        }
    }

    invert(lagrangeCoeffs, Vand, len);

    for (int i = 0; i < len; ++i) {
        delete[] Vand[i];
    }
    delete[] Vand;
}

arma::mat Element::computeMassMatrix(GaussQuadTri &quad) {
    arma::mat M(nBasis, nBasis);
    M.zeros();

    for (int k = 0; k < quad.nWeights; ++k) {
        double x = 0;
        double y = 0;
        for (int l = 0; l < nFaces; ++l) {
            x += quad.coords[k][l] * nodeX[l];
            y += quad.coords[k][l] * nodeY[l];
        }

        for (int i = 0; i < nBasis; ++i) {
            for (int j = i; j < nBasis; ++j) {
//                M[i*n+j] += quad.weights[k] * area * basis(i, x, y) * basis(j, x, y);
//                M[i][j] += quad.weights[k] * area * basis(i, x, y) * basis(j, x, y);
                M(i,j) += quad.weights[k] * area * basis(i, x, y) * basis(j, x, y);
            }
        }
    }

    for (int i = 0; i < nBasis; ++i) {
        for (int j = i + 1; j < nBasis; ++j) {
//            M[j*n+i] = M[i*n+j];
//            M[j][i] = M[i][j];
            M(j,i) = M(i,j);
        }
    }
    
    return M;
}

void Element::massMatrix(double **M, GaussQuadTri &quad) {
    //    std::fill(M, M + n*n, 0);
    for (int i = 0; i < nBasis; ++i) {
        std::fill(M[i], M[i] + nBasis, 0);
    }

    for (int k = 0; k < quad.nWeights; ++k) {
        double x = 0;
        double y = 0;
        for (int l = 0; l < nFaces; ++l) {
            x += quad.coords[k][l] * nodeX[l];
            y += quad.coords[k][l] * nodeY[l];
        }

        for (int i = 0; i < nBasis; ++i) {
            for (int j = i; j < nBasis; ++j) {
                //                M[i*n+j] += quad.weights[k] * area * basis(i, x, y) * basis(j, x, y);
                M[i][j] += quad.weights[k] * area * basis(i, x, y) * basis(j, x, y);
            }
        }
    }

    for (int i = 0; i < nBasis; ++i) {
        for (int j = i + 1; j < nBasis; ++j) {
            //            M[j*n+i] = M[i*n+j];
            M[j][i] = M[i][j];
        }
    }
}

void Element::artificialDamping(DVector &W, GaussQuadTri &quad,
                                DVector &centreX, DVector &centreY) {
    int pos = index * nEqns * nBasis;

    double averageRho = W[pos];

    double x, y;
    double integral1 = 0;
    double integral2 = 0;

    for (int k = 0; k < quad.nWeights; ++k) {
        quad.getCoords(x, y, k, nodeX, nodeY);
        double rho = 0;
        for (int i = 0; i < nBasis; ++i) {
            rho += W[pos + i] * basis(i, x, y);
        }
        integral1 += quad.weights[k] * rho * rho;
        integral2 += quad.weights[k] * (rho - averageRho) * (rho - averageRho);
    }
    integral1 *= area;
    integral2 *= area;

    double rez = integral2 / integral1;

    double f = 1 - (std::atan(numParam->dampConst * (rez - numParam->dampTol)) + M_PI / 2)
                 / M_PI * (1 - numParam->dampMin);

    if (rez > numParam->dampTol) {
        centreX.push_back(this->centreX);
        centreY.push_back(this->centreY);
    }

    for (int d = 0; d < nEqns; ++d) {
        for (int i = 1; i < nBasis; ++i) {
            W[pos + d * nBasis + i] *= f;
        }
    }
}

void Element::derivatives(const Vector4 &w, const Matrix42 &gradW, double p,
                          Vector2 &gradU, Vector2 &gradV, Vector2 &gradPOverRho) 
{
    double rho = w[0];
    double u = w[1] / rho;
    double v = w[2] / rho;
    
    /* i=0 is the derivative with respect to x and i=1 is the derivative with respect to y*/
    for (int i = 0; i < 2; ++i) {
        double dRho = gradW(0,i);
        gradU(i) = (gradW(1,i) - dRho*u) / rho;
        gradV(i) = (gradW(2,i) - dRho*v) / rho;

        double dp = (param->kapa-1) * (gradW(3,i) - 0.5*dRho*(u*u+v*v)
                                      - rho*(u*gradU(i)+v*gradV(i)));

        gradPOverRho(i) = (dp*rho - p*dRho) / (rho*rho);
    }
}

void Element::viscousFluxes(Vector4 &fluxX, Vector4 &fluxY, const Vector4 &w,
                            const Matrix42 &gradW)
{    
    double p = (param->kapa-1) * (w[3] - 0.5 * (w[1]*w[1] + w[2]*w[2]) / w[0]);
    
    Vector2 gradU, gradV, gradPOverRho;    
    derivatives(w, gradW, p, gradU, gradV, gradPOverRho);
        
    arma::mat::fixed<2,2> stress;    
    stressTensor(stress, gradU, gradV);            
    
    double constant = param->kapa/(param->kapa-1)/param->prandtl;  // = 4.8611
    double u = w[1] / w[0];
    double v = w[2] / w[0];
    
    fluxX[0] = 0;
    fluxX[1] = stress(0,0);
    fluxX[2] = stress(0,1);
    fluxX[3] = u * stress(0,0) + v * stress(0,1) + constant * gradPOverRho(0);
    
    fluxY[0] = 0;
    fluxY[1] = stress(1,0);
    fluxY[2] = stress(1,1);
    fluxY[3] = u * stress(1,0) + v * stress(1,1) + constant * gradPOverRho(1);
}

void Element::computeWallW(Vector4 &wallW, const Vector4 &w)
{    
    wallW[0] = w[0];
    wallW[1] = 0.0;
    wallW[2] = 0.0;
    wallW[3] = w[3] - 0.5 * (w[1]*w[1] + w[2]*w[2]) / w[0];
}

void Element::stressTensor(arma::mat::fixed<2,2> &stress, const Vector2 &gradU, const Vector2 &gradV)
{
    stress(0,0) = 4.0/3.0 * gradU(0) - 2.0/3.0 * gradV(1);
    stress(0,1) = gradU(1) + gradV(0);
    stress(1,0) = stress(0,1);
    stress(1,1) = 4.0/3.0 * gradV(1) - 2.0/3.0 * gradU(0);
}

arma::mat& Element::volumeIntegral(arma::mat &integral, arma::mat &W, GaussQuadTri &quad) {
    integral.zeros();
    
    Vector4 w;
    Vector4 fluxX;
    Vector4 fluxY;
    Vector4 tempFluxX;
    Vector4 tempFluxY;
    double x, y;    
    
    for (int k = 0; k < quad.nWeights; ++k) {
        quad.getCoords(x, y, k, nodeX, nodeY);

        getW(w, x, y, W);  // carry out linear combination
        physicalFluxes(fluxX, fluxY, w);        
        if (param->isFlowViscous) {
            Matrix42 gradW;
            getGradW(gradW, x, y, W);
            viscousFluxes(tempFluxX, tempFluxY, w, gradW);
            fluxX -= 1 / param->reynolds * tempFluxX;
            fluxY -= 1 / param->reynolds * tempFluxY;
        }
        
        for (int d = 0; d < nEqns; ++d) {
            for (int i = 0; i < nBasis; ++i) {
                integral(i,d) += quad.weights[k] * (fluxX(d) * derBasisX(i, x, y)
                               + fluxY(d) * derBasisY(i, x, y)) * area;
            }
        }
    }
    
    return integral;
}

void Element::volumeIntegral(double *integral, DVector &W, GaussQuadTri &quad) {
    double w[nEqns];
    double gradW[2][nEqns];
    double fluxX[nEqns], fluxY[nEqns];
    double x, y;

    for (int k = 0; k < quad.nWeights; ++k) {
        quad.getCoords(x, y, k, nodeX, nodeY);

        combine(w, x, y, &W[index * nEqns * nBasis]);
        physicalFluxes(fluxX, fluxY, w);
        
        if (param->isFlowViscous) {
            getGradW(gradW, x, y, &W[index * nEqns * nBasis]);
            viscousFluxes(fluxX, fluxY, w, gradW);            
        }
        
        for (int d = 0; d < nEqns; ++d) {
            for (int i = 0; i < nBasis; ++i) {
                integral[d * nBasis + i] += quad.weights[k] * (fluxX[d] * derBasisX(i, x, y)
                        + fluxY[d] * derBasisY(i, x, y)) * area;
            }
        }
    }
}

arma::mat& Element::boundaryLineIntegral(arma::mat &integral, int faceNumber, arma::mat &W1,
                                         ElemVector &mesh, GaussQuad &quad)
{
    integral.zeros();
    
    Vector4 w1;
    Vector4 w2;
    Vector4 flux;
    Vector4 tempFlux;    
    for (int k = 0; k < quad.nPoints; ++k) {
        double x = quad.getCoord(k, nodeX[faceNumber], nodeX[(faceNumber + 1) % nFaces]);
        double y = quad.getCoord(k, nodeY[faceNumber], nodeY[(faceNumber + 1) % nFaces]);

        getW(w1, x, y, W1);
        flux = boundary(tempFlux, w1, edge[faceNumber]);            

        if (param->isFlowViscous && neighbour[faceNumber] == WALL) {
            Matrix42 gradW1;
            getGradW(gradW1, x, y, W1);

            viscousWall(tempFlux, edge[faceNumber], w1, gradW1);                
            flux -= tempFlux / param->reynolds;
            
            computeWallW(w2, w1);
            flux -= numParam->penalty * (w2 - w1);  // penalty
        }                               

        for (int d = 0; d < nEqns; ++d) {
            for (int i = 0; i < nBasis; ++i) {
                integral(i, d) -= quad.weights[k] * flux(d) * basis(i, x, y) * edge[faceNumber].length / 2;
            }
        }
    }
    
    return integral;
}

arma::mat& Element::internalLineIntegral(arma::mat &integral, int faceNumber, arma::mat &W1,
                                 arma::mat &W2, ElemVector &mesh, GaussQuad &quad)
{
    integral.zeros();
    
    Vector4 w1;
    Vector4 w2;
    Vector4 flux;
    Vector4 tempFlux;    
    for (int k = 0; k < quad.nPoints; ++k) {
        double x = quad.getCoord(k, nodeX[faceNumber], nodeX[(faceNumber + 1) % nFaces]);
        double y = quad.getCoord(k, nodeY[faceNumber], nodeY[(faceNumber + 1) % nFaces]);

        getW(w1, x, y, W1);
        mesh[neighbour[faceNumber]]->getW(w2, x, y, W2);
        flux = vanLeer(tempFlux, w1, w2, edge[faceNumber]);

        if (param->isFlowViscous) {
            Matrix42 gradW1, gradW2;
            getGradW(gradW1, x, y, W1);
            mesh[neighbour[faceNumber]]->getGradW(gradW2, x, y, W2);

            flux -= 1 / param->reynolds * viscousNumericalFlux(tempFlux, edge[faceNumber],
                                          w1, w2, gradW1, gradW2);

            flux -= numParam->penalty * (w2 - w1);  // penalty
        }                                           

        for (int d = 0; d < nEqns; ++d) {
            for (int i = 0; i < nBasis; ++i) {
                integral(i, d) -= quad.weights[k] * flux(d) * basis(i, x, y) * edge[faceNumber].length / 2;
            }
        }
    }
    
    return integral;
}


void Element::lineIntegral(double *integral, int faceNumber, DVector &W, ElemVector &mesh, GaussQuad &quad) {
    double flux[nEqns];
    double w1[nEqns];
    double w2[nEqns];

    for (int k = 0; k < quad.nPoints; ++k) {
        double x = quad.getCoord(k, nodeX[faceNumber], nodeX[(faceNumber + 1) % nFaces]);
        double y = quad.getCoord(k, nodeY[faceNumber], nodeY[(faceNumber + 1) % nFaces]);

        combine(w1, x, y, &W[index * nEqns * nBasis]);
        if (neighbour[faceNumber] > BOUNDARY) {
            mesh[neighbour[faceNumber]]->combine(w2, x, y, &W[neighbour[faceNumber] * nEqns * nBasis]);
            vanLeer(flux, w1, w2, edge[faceNumber]);
            
            if (param->isFlowViscous) {
                double gradW1[2][nEqns], gradW2[2][nEqns];
                getGradW(gradW1, x, y, &W[index * nEqns * nBasis]);
                mesh[neighbour[faceNumber]]->getGradW(gradW2, x, y, &W[neighbour[faceNumber] * nEqns * nBasis]);
                viscousNumericalFlux(flux, w1, gradW1, w2, gradW2, edge[faceNumber]);   
            }
            
        } else {
            if (!param->isFlowViscous) {
            boundary(flux, w1, edge[faceNumber]);
            } else {
                double gradW1[2][nEqns];
                getGradW(gradW1, x, y, &W[index * nEqns * nBasis]);
                boundary(flux, w1, gradW1, edge[faceNumber]);
            }
        }

        for (int d = 0; d < nEqns; ++d) {
            for (int i = 0; i < nBasis; ++i) {
                integral[d * nBasis + i] -= quad.weights[k] * flux[d] * basis(i, x, y) * edge[faceNumber].length / 2;
            }
        }
    }
}

/**
 * Compute residual vector R(W) = M^-1 * ( volumeIntegral - lineIntegral )
 * 
 * @param RW residual vector which is the output parameter
 * @param W 
 * @param mesh
 * @param quad
 * @param quadTri
 */
void Element::residualVector(arma::mat &RW, arma::cube &W, ElemVector &mesh,
                              GaussQuad &quad, GaussQuadTri &quadTri)
{    
    arma::mat integral(nBasis, nEqns);    
    arma::mat tempIntegral(nBasis, nEqns);
    arma::mat &W1 = W.slice(index);
    
    integral = volumeIntegral(tempIntegral, W1, quadTri);    
    for (int j = 0; j < nFaces; ++j) {
        if (neighbour[j] > BOUNDARY) {
            integral += internalLineIntegral(tempIntegral, j, W1, W.slice(neighbour[j]), mesh, quad);
        } else {
            integral += boundaryLineIntegral(tempIntegral, j, W1, mesh, quad);
        }
    }
       
    arma::mat invMassMatrix = arma::inv(computeMassMatrix(quadTri));   
    
    for (int d = 0; d < nEqns; ++d) {
        RW.unsafe_col(d) = invMassMatrix * integral.unsafe_col(d);
    }
}

void Element::residualVector(DVector &RW, DVector &W, ElemVector &mesh, GaussQuad &quad, GaussQuadTri &quadTri) {
    double integral[nEqns * nBasis];
    std::fill(integral, integral + nEqns * nBasis, 0);
    
    volumeIntegral(integral, W, quadTri);
    for (int j = 0; j < nFaces; ++j) {
        lineIntegral(integral, j, W, mesh, quad);
    }    

    double** M = new double*[nBasis];
    for (int i = 0; i < nBasis; ++i) {
        M[i] = new double[nBasis];
    }
    massMatrix(M, quadTri);
    invert(M, M, nBasis);
    
    for (int d = 0; d < nEqns; ++d) {
        matrixTimesVector(&RW[index * nEqns * nBasis + d * nBasis], M, &integral[d * nBasis], nBasis, nBasis);
    }

    for (int i = 0; i < nBasis; ++i) {
        delete[] M[i];
    }
    delete[] M;
}

/* TO DO: place in sparse matrix class */
void Element::insertJacobiBlock(SparseMatrix &jacobi, arma::mat &jacobiBlock,
                                int rowNum, int colNum)
{
    for (int j = 0; j < jacobiBlock.n_cols; ++j) {
        for (int i = 0; i < jacobiBlock.n_rows; ++i) {
            jacobi.insert(rowNum + i, colNum + j, jacobiBlock(i,j));
        }
    }
}
void Element::jacobi(SparseMatrix &jacobi, arma::mat &RW, arma::cube &W, ElemVector &mesh,
                      double timeStep, GaussQuad &quad, GaussQuadTri &quadTri)
{    
    arma::cube lineInt(nBasis, nEqns, nFaces);    
    arma::mat W1 = W.slice(index);  // copy basis coefficients of the current element into W1
    
    /* compute residual vector R(W) */    
    volumeIntegral(RW, W1, quadTri);    
    for (int j = 0; j < nFaces; ++j) {
        if (neighbour[j] > BOUNDARY) {
            internalLineIntegral(lineInt.slice(j), j, W1, W.slice(neighbour[j]), mesh, quad);            
        } else {
            boundaryLineIntegral(lineInt.slice(j), j, W1, mesh, quad);
        }
        
        RW += lineInt.slice(j);
    }
    
    double h = numParam->diffTol;    
    arma::mat jacobiBlock(nBasis * nEqns, nBasis * nEqns);    
    arma::mat RW_h(nBasis, nEqns); /* residual vector R(W + h*e_r) */            
    
    arma::mat tempIntegral(nBasis, nEqns);    
    
    /* compute residual vector R(W + h*e_r) */
    for (int r = 0; r < W1.n_elem; ++r) {
        W1(r) += h;
        
        volumeIntegral(RW_h, W1, quadTri);    
        for (int j = 0; j < nFaces; ++j) {
            if (neighbour[j] > BOUNDARY) {
                internalLineIntegral(tempIntegral, j, W1, W.slice(neighbour[j]), mesh, quad);
            } else {
                boundaryLineIntegral(tempIntegral, j, W1, mesh, quad);                
            }
            RW_h += tempIntegral;
        }                                    
        
        /* insert "-(R(W) - R(W + h*e_r)) / h" into local Jacobi matrix */
        jacobiBlock.col(r) = arma::vectorise( -(RW_h-RW)/h );
        
        W1(r) -= h;
    }
    
    /* add mass matrix to the diagonal blocks of the local Jacobi matrix */
    arma::mat massMatrix = computeMassMatrix(quadTri);
    for (int d = 0; d < nEqns; ++d) {
        jacobiBlock(arma::span(d*nBasis, (d+1)*nBasis-1), arma::span(d*nBasis,(d+1)*nBasis-1))
                   += massMatrix / timeStep;
    }    
    
    omp_set_lock(lock);
    insertJacobiBlock(jacobi, jacobiBlock, index*nBasis*nEqns, index*nBasis*nEqns);
    omp_unset_lock(lock);
        
    /*** *** ***/    
    
    for (int l = 0; l < nFaces; ++l) {
        if (neighbour[l] <= BOUNDARY) {
            continue;
        }
        
        arma::mat W2 = W.slice(neighbour[l]);

        for (int r = 0; r < W2.n_elem; ++r) {                        
            W2(r) += h;

            internalLineIntegral(tempIntegral, l, W1, W2, mesh, quad);   
            jacobiBlock.col(r) = arma::vectorise( -(tempIntegral - lineInt.slice(l)) / h );

            W2(r) -= h;
        }
        
//        std::cout << jacobiBlock << std::endl;
        omp_set_lock(lock);
        insertJacobiBlock(jacobi, jacobiBlock, index*nBasis*nEqns, neighbour[l]*nBasis*nEqns);
        omp_unset_lock(lock);
    }
}

void Element::jacobi(SparseMatrix &jacobi, DVector &RW, DVector &W, ElemVector &mesh,
        double timeStep, GaussQuad &quad, GaussQuadTri &quadTri) {
    /*
     * mass matrix
     */
    //    double currentTime;
    //    std::clock_t startTime = std::clock();

    double** M = new double*[nBasis];
    for (int i = 0; i < nBasis; ++i) {
        M[i] = new double[nBasis];
    }

    //    currentTime = (std::clock() - startTime) / (double) CLOCKS_PER_SEC;
    //    std::cout << currentTime << "\n";
    massMatrix(M, quadTri);

    /*
     * residual vector
     */
    double integral[nEqns * nBasis];
    std::fill(integral, integral + nEqns * nBasis, 0);
    double lineInt[nFaces][nEqns * nBasis];
    for (int j = 0; j < nFaces; ++j) {
        std::fill(lineInt[j], lineInt[j] + nEqns * nBasis, 0);
    }

    volumeIntegral(integral, W, quadTri);
    for (int j = 0; j < nFaces; ++j) {
        lineIntegral(lineInt[j], j, W, mesh, quad);
    }
    for (int i = 0; i < nEqns * nBasis; ++i) {
        for (int j = 0; j < nFaces; ++j) {
            integral[i] += lineInt[j][i];
        }
    }

    for (int i = 0; i < nEqns * nBasis; ++i) {
        RW[index * nEqns * nBasis + i] = integral[i];
    }
    
    /*
     * 
     */
    double h = numParam->diffTol;
    double integral0[nEqns * nBasis];
    double localJacobi[nEqns * nBasis][nEqns * nBasis];
    for (int r = 0; r < nEqns * nBasis; ++r) {
        W[index * nEqns * nBasis + r] += h;

        std::fill(integral0, integral0 + nEqns * nBasis, 0);

        volumeIntegral(integral0, W, quadTri);
        for (int j = 0; j < nFaces; ++j) {
            lineIntegral(integral0, j, W, mesh, quad);
        }

        for (int l = 0; l < nEqns * nBasis; ++l) {
            localJacobi[l][r] = -(integral0[l] - integral[l]) / h;
        }

        W[index * nEqns * nBasis + r] -= h;
    }
    
    for (int d = 0; d < nEqns; ++d) {
        for (int i = 0; i < nBasis; ++i) {
            for (int j = 0; j < nBasis; ++j) {
                localJacobi[d * nBasis + i][d * nBasis + j] += M[i][j] / timeStep;
            }
        }
    }        

    for (int r = 0; r < nEqns * nBasis; ++r) {
        for (int l = 0; l < nEqns * nBasis; ++l) {
            jacobi.insert(index * nEqns * nBasis + l, index * nEqns * nBasis + r, localJacobi[l][r]);
        }
    }        

    /*
     * 
     */
    for (int j = 0; j < nFaces; ++j) {
        if (neighbour[j] > BOUNDARY) {
            for (int r = 0; r < nEqns * nBasis; ++r) {
                W[neighbour[j] * nEqns * nBasis + r] += h;

                std::fill(integral0, integral0 + nEqns * nBasis, 0);
                lineIntegral(integral0, j, W, mesh, quad);
                
                for (int l = 0; l < nEqns * nBasis; ++l) {
                    double value = -(integral0[l] - lineInt[j][l]) / h;
                    jacobi.insert(index * nEqns * nBasis + l, neighbour[j] * nEqns * nBasis + r, value);
                }

                W[neighbour[j] * nEqns * nBasis + r] -= h;
            }
        }
    }    

    //////////////////////////////////////////////////////////

    //    currentTime = (std::clock() - startTime) / (double) CLOCKS_PER_SEC;
    //    printf("%.1e\n", currentTime);

    for (int i = 0; i < nBasis; ++i) {
        delete M[i];
    }
    delete M;
}

arma::vec& Element::viscousNumericalFlux(Vector4 &flux, const Edge &edge, const Vector4 &w1,
                           const Vector4 &w2, const Matrix42 &gradW1, const Matrix42 &gradW2)
                                                                                  
{
    Vector4 flux1X, flux1Y, flux2X, flux2Y;
    viscousFluxes(flux1X, flux1Y, w1, gradW1);
    viscousFluxes(flux2X, flux2Y, w2, gradW2);
    
    flux = 0.5  *  ((flux1X * edge.normalX + flux1Y * edge.normalY)
                  + (flux2X * edge.normalX + flux2Y * edge.normalY));
    return flux;
}

void Element::viscousNumericalFlux(double (&flux)[nEqns], double (&w1)[nEqns], double (&gradW1)[2][nEqns],
                                   double (&w2)[nEqns], double (&gradW2)[2][nEqns], const Edge& edge) {
    double nx = edge.normalX;
    double ny = edge.normalY;
    
    double fluxX1[nEqns], fluxY1[nEqns], fluxX2[nEqns], fluxY2[nEqns];
    
    std::fill(fluxX1, fluxX1 + nEqns, 0);
    std::fill(fluxY1, fluxY1 + nEqns, 0);
    std::fill(fluxX2, fluxX2 + nEqns, 0);
    std::fill(fluxY2, fluxY2 + nEqns, 0);
    
    viscousFluxes(fluxX1, fluxY1, w1, gradW1);
    viscousFluxes(fluxX2, fluxY2, w2, gradW2);
    for (int d = 0; d < nEqns; ++d) {
        flux[d] += 0.5 * ((fluxX1[d] * nx + fluxY1[d] * ny) + (fluxX2[d] * nx + fluxY2[d] * ny))
                 - numParam->penalty * (w2[d] - w1[d]);
    }
    
//    std::fill(fluxX, fluxX + dim, 0);
//    std::fill(fluxY, fluxY + dim, 0);
//    viscousFluxes(fluxX, fluxY, w2, gradW2);
//    for (int d = 0; d < dim; ++d) {
//        flux[d] += 0.5 * (fluxX[d] * nx + fluxY[d] * ny);
//    }
}

arma::vec& Element::vanLeer(arma::vec &flux, const arma::vec &w1, const arma::vec &w2, const Edge& edge)
{        
    double nx = edge.normalX;
    double ny = edge.normalY;
    double kapa = param->kapa;
    
    // konstanty
    double rho, u, v, Vn, q2, M, a, p;
    rho = w1[0];
    u = w1[1] / rho;
    v = w1[2] / rho;
    Vn = u * nx + v * ny; // normalova rychlost
    q2 = u * u + v*v;
    p = (kapa - 1) * (w1[3] - 0.5 * rho * (u * u + v * v));
//    if (p < numParam->tol) {
//        std::cerr << "pressure is too low" << std::endl;
//        p = numParam->tol;
//    }
    a = std::sqrt(kapa * p / rho); // rychlost zvuku
    M = Vn / a;

    double fPlus[nEqns];
    if (std::abs(M) < 1) {
        double fm = rho * a * ((M + 1)*(M + 1)) / 4; // mass flux
        fPlus[0] = fm;
        fPlus[1] = fm * (u + (-Vn + 2 * a) / kapa * nx);
        fPlus[2] = fm * (v + (-Vn + 2 * a) / kapa * ny);
        fPlus[3] = fm * ((q2 - Vn * Vn) / 2 + (((kapa - 1) * Vn + 2 * a))*(((kapa - 1) * Vn
                + 2 * a)) / (2 * (kapa * kapa - 1)));
    } else if (M >= 1) {
        fPlus[0] = w1[0] * Vn;
        fPlus[1] = w1[1] * Vn + p * nx;
        fPlus[2] = w1[2] * Vn + p * ny;
        fPlus[3] = (w1[3] + p) * Vn;
    } else {
        fPlus[0] = 0;
        fPlus[1] = 0;
        fPlus[2] = 0;
        fPlus[3] = 0;
    }

    //    Edge &edge2 = *edge.neighEdge;
    rho = w2[0];
    u = w2[1] / rho;
    v = w2[2] / rho;
    Vn = u * nx + v * ny; // normalova rychlost
    q2 = u * u + v*v;
    p = (kapa - 1) * (w2[3] - 0.5 * rho * (u * u + v * v));
//    if (p < numParam->tol) {
//        std::cerr << "pressure is too low" << std::endl;
//        p = numParam->tol;
//    }
    a = std::sqrt(kapa * p / rho); // rychlost zvuku
    M = Vn / a;

    double fMinus[nEqns];
    if (std::abs(M) < 1) {
        double fm = -rho * a * ((M - 1)*(M - 1)) / 4; // mass flux
        fMinus[0] = fm;
        fMinus[1] = fm * (u + (-Vn - 2 * a) / kapa * nx);
        fMinus[2] = fm * (v + (-Vn - 2 * a) / kapa * ny);
        fMinus[3] = fm * ((q2 - Vn * Vn) / 2 + (((kapa - 1) * Vn - 2 * a)*((kapa - 1) * Vn
                - 2 * a)) / (2 * (kapa * kapa - 1)));
    } else if (M <= -1) {
        fMinus[0] = w2[0] * Vn;
        fMinus[1] = w2[1] * Vn + p * nx;
        fMinus[2] = w2[2] * Vn + p * ny;
        fMinus[3] = (w2[3] + p) * Vn;
    } else {
        fMinus[0] = 0;
        fMinus[1] = 0;
        fMinus[2] = 0;
        fMinus[3] = 0;
    }

    for (int d = 0; d < nEqns; ++d) {
        flux[d] = fPlus[d] + fMinus[d];
    }
    
    return flux;
}

void Element::vanLeer(double (&flux)[nEqns], double *w1, double *w2, const Edge& edge) {
    // konstanty
    double nx = edge.normalX;
    double ny = edge.normalY;

    double kapa = param->kapa;
    double rho, u, v, Vn, q2, M, a, p;
    rho = w1[0];
    u = w1[1] / rho;
    v = w1[2] / rho;
    Vn = u * nx + v * ny; // normalova rychlost
    q2 = u * u + v*v;
    p = (kapa - 1) * (w1[3] - 0.5 * rho * (u * u + v * v));
//    if (p < numParam->tol) {
//        std::cerr << "pressure is too low" << std::endl;
//        p = numParam->tol;
//    }
    a = std::sqrt(kapa * p / rho); // rychlost zvuku
    M = Vn / a;

    double fPlus[nEqns];
    if (std::abs(M) < 1) {
        double fm = rho * a * ((M + 1)*(M + 1)) / 4; // mass flux
        fPlus[0] = fm;
        fPlus[1] = fm * (u + (-Vn + 2 * a) / kapa * nx);
        fPlus[2] = fm * (v + (-Vn + 2 * a) / kapa * ny);
        fPlus[3] = fm * ((q2 - Vn * Vn) / 2 + (((kapa - 1) * Vn + 2 * a))*(((kapa - 1) * Vn
                + 2 * a)) / (2 * (kapa * kapa - 1)));
    } else if (M >= 1) {
        fPlus[0] = w1[0] * Vn;
        fPlus[1] = w1[1] * Vn + p * nx;
        fPlus[2] = w1[2] * Vn + p * ny;
        fPlus[3] = (w1[3] + p) * Vn;
    } else {
        fPlus[0] = 0;
        fPlus[1] = 0;
        fPlus[2] = 0;
        fPlus[3] = 0;
    }

    //    Edge &edge2 = *edge.neighEdge;
    rho = w2[0];
    u = w2[1] / rho;
    v = w2[2] / rho;
    Vn = u * nx + v * ny; // normalova rychlost
    q2 = u * u + v*v;
    p = (kapa - 1) * (w2[3] - 0.5 * rho * (u * u + v * v));
//    if (p < numParam->tol) {
//        std::cerr << "pressure is too low" << std::endl;
//        p = numParam->tol;
//    }
    a = std::sqrt(kapa * p / rho); // rychlost zvuku
    M = Vn / a;

    double fMinus[nEqns];
    if (std::abs(M) < 1) {
        double fm = -rho * a * ((M - 1)*(M - 1)) / 4; // mass flux
        fMinus[0] = fm;
        fMinus[1] = fm * (u + (-Vn - 2 * a) / kapa * nx);
        fMinus[2] = fm * (v + (-Vn - 2 * a) / kapa * ny);
        fMinus[3] = fm * ((q2 - Vn * Vn) / 2 + (((kapa - 1) * Vn - 2 * a)*((kapa - 1) * Vn
                - 2 * a)) / (2 * (kapa * kapa - 1)));
    } else if (M <= -1) {
        fMinus[0] = w2[0] * Vn;
        fMinus[1] = w2[1] * Vn + p * nx;
        fMinus[2] = w2[2] * Vn + p * ny;
        fMinus[3] = (w2[3] + p) * Vn;
    } else {
        fMinus[0] = 0;
        fMinus[1] = 0;
        fMinus[2] = 0;
        fMinus[3] = 0;
    }

    for (int d = 0; d < nEqns; ++d) {
        flux[d] = fPlus[d] + fMinus[d];
    }
}

void Element::AUSM(double (&flux)[nEqns], double *w1, double *w2, const Edge& edge) {
    //    double p1 = elem1.pressure;
    //    double p2 = elem2.pressure;
    //    double a1 = elem1.speedOfSound;
    //    double a2 = elem2.speedOfSound;
    double kapa = param->kapa;
    double p1 = (kapa - 1) * (w1[3] - 0.5 * (w1[1] * w1[1] + w1[2] * w1[2]) / w1[0]);
    double a1 = std::sqrt(kapa * p1 / w1[0]);
    double p2 = (kapa - 1) * (w2[3] - 0.5 * (w2[1] * w2[1] + w2[2] * w2[2]) / w2[0]);
    double a2 = std::sqrt(kapa * p2 / w2[0]);

    double nx = edge.normalX;
    double ny = edge.normalY;

    //    Edge edge2 = *edge.neighEdge;

    double velocity1, velocity2, mach1, mach2;
    velocity1 = (w1[1] * nx + w1[2] * ny) / w1[0];
    velocity2 = (w2[1] * nx + w2[2] * ny) / w2[0];
    mach1 = velocity1 / a1;
    mach2 = velocity2 / a2;

    double splitMach1, P1;
    if (mach1 > 1 || mach1 < -1) {
        splitMach1 = (mach1 + std::abs(mach1)) / 2;
        P1 = splitMach1 / mach1;
    } else {
        splitMach1 = (mach1 + 1)*(mach1 + 1) / 4 + (mach1 * mach1 - 1)*(mach1 * mach1 - 1) / 4;
        //        P1 = (mach1+1)*(mach1+1)/4 * (2-mach1);
        P1 = (1 + mach1) / 2;
    }

    double splitMach2, P2;
    if (mach2 > 1 || mach2 < -1) {
        splitMach2 = (mach2 - std::abs(mach2)) / 2;
        P2 = splitMach2 / mach2;
    } else {
        splitMach2 = -(mach2 - 1)*(mach2 - 1) / 4 - (mach2 * mach2 - 1)*(mach2 * mach2 - 1) / 4;
        //        P2 = -(mach2-1)*(mach2-1)/4 * (2+mach2);
        P2 = (1 - mach2) / 2;
    }

    double p = P1 * p1 + P2 * p2;
    double mach = splitMach1 + splitMach2;
    double machPlus = std::max(mach, 0.);
    double machMinus = std::min(mach, 0.);

    double flux1[nEqns];
    flux1[0] = a1 * machPlus * w1[0];
    flux1[1] = a1 * machPlus * w1[1];
    flux1[2] = a1 * machPlus * w1[2];
    flux1[3] = a1 * machPlus * (w1[3] + p1);

    double flux2[nEqns];
    flux2[0] = a2 * machMinus * w2[0];
    flux2[1] = a2 * machMinus * w2[1];
    flux2[2] = a2 * machMinus * w2[2];
    flux2[3] = a2 * machMinus * (w2[3] + p2);

    flux[0] = flux1[0] + flux2[0];
    flux[1] = flux1[1] + flux2[1] + p * nx;
    flux[2] = flux1[2] + flux2[2] + p * ny;
    flux[3] = flux1[3] + flux2[3];
}

void Element::boundary(double (&flux)[nEqns], double (&w1)[nEqns], double (&gradW)[2][nEqns], const Edge& edge) {
    boundary(flux, w1, edge);
    
    if (edge.boundaryType == WALL) {
        double gradRho[2], gradU[2], gradV[2], gradPOverRho[2];
        grad(gradRho, gradU, gradV, gradPOverRho, w1, gradW);
        double tensor[2][2];
        stressTensor(tensor, gradU, gradV);
        flux[1] -= 1. / param->reynolds * (tensor[0][0] * edge.normalX + tensor[0][1] * edge.normalY);
        flux[2] -= 1. / param->reynolds * (tensor[1][0] * edge.normalX + tensor[1][1] * edge.normalY);
        
        double p = (param->kapa - 1) * (w1[3] - 0.5 * (w1[1] * w1[1] + w1[2] * w1[2]) / w1[0]);
//        if (p < numParam->tol) {
//            std::cerr << "pressure is too low" << std::endl;
//            p = numParam->tol;
//        }
        
        double w2[nEqns];
        w2[0] = w1[0];
        w2[1] = 0.;
        w2[2] = 0.;
        w2[3] = p / (param->kapa - 1);
//        w2[3] = w1[3] - 0.5 * (w1[1] * w1[1] + w1[2] * w1[2]) / w1[0];
        
        for (int d = 0; d < nEqns; ++d) {
            flux[d] -= numParam->penalty * (w2[d] - w1[d]);
        }
    }
}

arma::vec& Element::viscousWall(Vector4 &flux, const Edge &edge, const Vector4 &w,
                                const Matrix42 &gradW)
{
    double rho = w[0];
    double u = w[1] / rho;
    double v = w[2] / rho;
    
    Vector2 gradU, gradV;
    for (int i = 0; i < 2; ++i) {
        double dRho = gradW(0,i);
        gradU(i) = (gradW(1,i) - dRho*u) / rho;
        gradV(i) = (gradW(2,i) - dRho*v) / rho;
    }
    
    arma::mat::fixed<2,2> stress;
    stressTensor(stress, gradU, gradV);
    
    flux[0] = 0;
    flux[1] = stress(0,0) * edge.normalX + stress(0,1) * edge.normalY;
    flux[2] = stress(1,0) * edge.normalX + stress(1,1) * edge.normalY;
    flux[3] = 0;
}

arma::vec& Element::boundary(Vector4 &flux, const Vector4 &w, const Edge& edge) {
    double kapa = param->kapa;
    if (edge.boundaryType == WALL) {
        double p = (kapa - 1) * (w[3] - 0.5 * (w[1] * w[1] + w[2] * w[2]) / w[0]);
//        if (p < numParam->tol) {
//            std::cerr << "pressure is too low" << std::endl;
//            p = numParam->tol;
//        }
        flux[0] = 0;
        flux[1] = p * edge.normalX;
        flux[2] = p * edge.normalY;
        flux[3] = 0;
    } else if (edge.boundaryType == INLET) {
        //double p = (kapa-1) * (w[3] - 0.5*(w[1]*w[1]+w[2]*w[2])/w[0]);
        double pIn = (kapa - 1) * (w[3] - 0.5 * (w[1] * w[1] + w[2] * w[2]) / w[0]);
//        if (pIn < numParam->tol) {
//            std::cerr << "pressure is too low" << std::endl;
//            pIn = numParam->tol;
//        }
        double machIn = sqrt(2 / (kapa - 1) * (pow(param->pIn0 / pIn, ((kapa - 1) / kapa)) - 1));
        double rhoIn = param->rhoIn0 * pow(1 + (kapa - 1) / 2 * machIn * machIn, 1 / (1 - kapa));
        double velocity = std::abs(machIn * sqrt(kapa * pIn / rhoIn));

        double u = velocity * cos(param->angleIn);
        double v = velocity * sin(param->angleIn);

        double E = pIn / (kapa - 1) + 0.5 * rhoIn * velocity * velocity;

        double w0[nEqns];
        w0[0] = rhoIn;
        w0[1] = rhoIn * u;
        w0[2] = rhoIn * v;
        w0[3] = E;

        double fluxX[nEqns], fluxY[nEqns];
        physicalFluxes(fluxX, fluxY, w0, pIn);

        for (int d = 0; d < nEqns; ++d) {
            flux[d] = fluxX[d] * edge.normalX + fluxY[d] * edge.normalY;
        }
    } else if (edge.boundaryType == OUTLET) {
        double Eout = param->pOut0 / (kapa - 1) + 0.5 * (w[1] * w[1] + w[2] * w[2]) / w[0];

        double w0[nEqns];
        w0[0] = w[0];
        w0[1] = w[1];
        w0[2] = w[2];
        w0[3] = Eout;

        double fluxX[nEqns], fluxY[nEqns];
        physicalFluxes(fluxX, fluxY, w0, param->pOut0);

        for (int d = 0; d < nEqns; ++d) {
            flux[d] = fluxX[d] * edge.normalX + fluxY[d] * edge.normalY;
        }
    }

    return flux;
}

void Element::boundary(double (&flux)[nEqns], double *w, const Edge& edge) {
    double kapa = param->kapa;
    if (edge.boundaryType == WALL) {
        double p = (kapa - 1) * (w[3] - 0.5 * (w[1] * w[1] + w[2] * w[2]) / w[0]);
//        if (p < numParam->tol) {
//            std::cerr << "pressure is too low" << std::endl;
//            p = numParam->tol;
//        }
        flux[0] = 0;
        flux[1] = p * edge.normalX;
        flux[2] = p * edge.normalY;
        flux[3] = 0;
    } else if (edge.boundaryType == INLET) {
        //double p = (kapa-1) * (w[3] - 0.5*(w[1]*w[1]+w[2]*w[2])/w[0]);
        double pIn = (kapa - 1) * (w[3] - 0.5 * (w[1] * w[1] + w[2] * w[2]) / w[0]);
//        if (pIn < numParam->tol) {
//            std::cerr << "pressure is too low" << std::endl;
//            pIn = numParam->tol;
//        }
        double machIn = sqrt(2 / (kapa - 1) * (pow(param->pIn0 / pIn, ((kapa - 1) / kapa)) - 1));
        double rhoIn = param->rhoIn0 * pow(1 + (kapa - 1) / 2 * machIn * machIn, 1 / (1 - kapa));
        double velocity = std::abs(machIn * sqrt(kapa * pIn / rhoIn));

        double u = velocity * cos(param->angleIn);
        double v = velocity * sin(param->angleIn);

        double E = pIn / (kapa - 1) + 0.5 * rhoIn * velocity * velocity;

        double w0[nEqns];
        w0[0] = rhoIn;
        w0[1] = rhoIn * u;
        w0[2] = rhoIn * v;
        w0[3] = E;

        double fluxX[nEqns], fluxY[nEqns];
        physicalFluxes(fluxX, fluxY, w0, pIn);

        for (int d = 0; d < nEqns; ++d) {
            flux[d] = fluxX[d] * edge.normalX + fluxY[d] * edge.normalY;
        }
    } else if (edge.boundaryType == OUTLET) {
        double Eout = param->pOut0 / (kapa - 1) + 0.5 * (w[1] * w[1] + w[2] * w[2]) / w[0];

        double w0[nEqns];
        w0[0] = w[0];
        w0[1] = w[1];
        w0[2] = w[2];
        w0[3] = Eout;

        double fluxX[nEqns], fluxY[nEqns];
        physicalFluxes(fluxX, fluxY, w0, param->pOut0);

        for (int d = 0; d < nEqns; ++d) {
            flux[d] = fluxX[d] * edge.normalX + fluxY[d] * edge.normalY;
        }
    }

}

//void Element::fFun(double (&flux)[dim], const double (&w)[dim], const double p) {
//    flux[0] = w[1];
//    flux[1] = w[1] * w[1] / w[0] + p;
//    flux[2] = w[1] * w[2] / w[0];
//    flux[3] = w[1]        / w[0] * (w[3]+p);
//}
//
//void Element::gFun(double (&flux)[dim], const double (&w)[dim], const double p) {
//    flux[0] = w[2];
//    flux[1] = w[2] * w[1] / w[0];
//    flux[2] = w[2] * w[2] / w[0] + p;
//    flux[3] = w[2]        / w[0] * (w[3]+p);
//}

void Element::physicalFluxes(arma::vec &fluxX, arma::vec &fluxY, arma::vec &w) {
    double p = (param->kapa - 1) * (w[3] - 0.5 * (w[1] * w[1] + w[2] * w[2]) / w[0]);
    physicalFluxes(fluxX, fluxY, w, p);
}

void Element::physicalFluxes(arma::vec &fluxX, arma::vec &fluxY, arma::vec &w, double p) {
    fluxX[0] = w[1];
    fluxX[1] = w[1] * w[1] / w[0] + p;
    fluxX[2] = w[1] * w[2] / w[0];
    fluxX[3] = w[1] / w[0] * (w[3] + p);

    fluxY[0] = w[2];
    fluxY[1] = w[2] * w[1] / w[0];
    fluxY[2] = w[2] * w[2] / w[0] + p;
    fluxY[3] = w[2] / w[0] * (w[3] + p);
}

void Element::physicalFluxes(double (&fluxX)[nEqns], double (&fluxY)[nEqns], double (&w)[nEqns]) {
    double p = (param->kapa - 1) * (w[3] - 0.5 * (w[1] * w[1] + w[2] * w[2]) / w[0]);
//    if (p < numParam->tol) {
//        std::cerr << "pressure is too low" << std::endl;
//        p = numParam->tol;
//    }
    physicalFluxes(fluxX, fluxY, w, p);
}

void Element::physicalFluxes(double (&fluxX)[nEqns], double (&fluxY)[nEqns], double (&w)[nEqns], double p) {
    fluxX[0] = w[1];
    fluxX[1] = w[1] * w[1] / w[0] + p;
    fluxX[2] = w[1] * w[2] / w[0];
    fluxX[3] = w[1] / w[0] * (w[3] + p);

    fluxY[0] = w[2];
    fluxY[1] = w[2] * w[1] / w[0];
    fluxY[2] = w[2] * w[2] / w[0] + p;
    fluxY[3] = w[2] / w[0] * (w[3] + p);
}

void Element::grad(double (&gradRho)[2], double (&gradU)[2], double (&gradV)[2], double (&gradPOverRho)[2],
        double (&w)[nEqns], double (&gradW)[2][nEqns]) {
    double rho = w[0];
    double u = w[1] / rho;
    double v = w[2] / rho;
    double E = w[3];

    double p = (param->kapa - 1) * (E - 0.5 * rho * (u * u + v * v));
    for (int k = 0; k < 2; ++k) {
        gradRho[k] = gradW[k][0];
        gradU[k] = 1. / rho * (gradW[k][1] - gradRho[k] * u);
        gradV[k] = 1. / rho * (gradW[k][2] - gradRho[k] * v);

        double derURhoU = u * (rho * gradU[k] + gradW[k][1]);
//        double derURhoU = gradRho[k]*u*u + rho*2*u*gradU[k];
        double derVRhoV = v * (rho * gradV[k] + gradW[k][2]);
//        double derVRhoV = gradRho[k]*v*v + rho*2*v*gradV[k];

        double derP = (param->kapa - 1) * (gradW[k][3] - 0.5 * (derURhoU + derVRhoV));

        gradPOverRho[k] = (derP * rho - p * gradRho[k]) / (rho * rho);
    }
}

void Element::stressTensor(double (&tensor)[2][2], double (&gradU)[2], double (&gradV)[2]) {
    tensor[0][0] = 4./3. * gradU[0] - 2./3. * gradV[1];
    tensor[0][1] = gradU[1] + gradV[0];
    tensor[1][0] = tensor[0][1];
    tensor[1][1] = 4./3. * gradV[1] - 2./3. * gradU[0];
}

void Element::viscousFluxes(double (&fluxX)[nEqns], double (&fluxY)[nEqns],
                            double (&w)[nEqns], double (&gradW)[2][nEqns]) {
    double gradRho[2], gradU[2], gradV[2], gradPOverRho[2];
    grad(gradRho, gradU, gradV, gradPOverRho, w, gradW);
    
    double tensor[2][2];
    stressTensor(tensor, gradU, gradV);
    
    double u = w[1] / w[0];
    double v = w[2] / w[0];
    double constant = param->kapa / ((param->kapa - 1.) * param->prandtl);
    double reynolds = param->reynolds;
    
    fluxX[1] -= 1/reynolds * tensor[0][0];
    fluxX[2] -= 1/reynolds * tensor[0][1];
    fluxX[3] -= 1/reynolds * (u * tensor[0][0] + v * tensor[0][1] + constant * gradPOverRho[0]);
    
    fluxY[1] -= 1/reynolds * tensor[1][0];
    fluxY[2] -= 1/reynolds * tensor[1][1];
    fluxY[3] -= 1/reynolds * (u * tensor[1][0] + v * tensor[1][1] + constant * gradPOverRho[1]);
}

/*
 * 
 * methods to be excluded from this class !!!
 * 
 */

/**
 * Multiplies a m x n matrix A by a vector b of length n and returns resulting vector x of length m.
 * 
 * @param x  resulting n x 1 vector
 * @param A  m x n matrix to be multiplied
 * @param b  m x 1 vector to be multiplied
 * @param m  number of rows of the matrix A
 * @param n  number of columns of the matrix A
 */
void Element::matrixTimesVector(double *x, double **A, double *b, int m, int n) {
    for (int i = 0; i < m; ++i) {
        double val = 0;
        for (int j = 0; j < n; ++j) {
            //            val += A[i*n+j] * b[j];
            val += A[i][j] * b[j];
        }
        x[i] = val;
    }
}

int Element::devideElement(double ***pCoords, int n) {
    int len = n * (n + 1) / 2;
    int nFaces = 3;
    double **coords = new double*[nFaces];
    for (int i = 0; i < nFaces; ++i) {
        coords[i] = new double[len];
    }

    double y[n];
    double x[n];
    x[0] = 0;
    y[0] = 0;
    for (int i = 1; i < n; ++i) {
        x[i] = x[i - 1] + 1. / (n - 1);
        y[i] = x[i];
    }

    int k = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            coords[0][k] = 1 - x[i];
            coords[1][k] = y[j];
            coords[2][k] = x[i] - y[j];
            ++k;
        }
    }

    *pCoords = coords;

    return len;
}


//void Element::print() {
//    printf("Triangle %d:\n", index);
//    printf("\tcoordinates: ");
//    for (int j = 0; j < nFaces; ++j) {
//        printf("[%.1f, %.1f], ", nodeX[j], nodeY[j]);
//    }
//
//    printf("\n\tindexes of neighbours: ");
//    for (int j = 0; j < nFaces; ++j) {
//        printf("%d ", neighbour[j]);
//    }
////    printf("\n\tindexes of neighbouring face: ");
////    for (int j = 0; j < nFaces; ++j) {
////        printf("%d ", neighFace[j]);
////    }
//    
//    printf("\n\tarea: %.2f", area);
//    
//    printf("\n\tcentre: [%.2f %.2f]", centreX, centreY);
//
//    printf("\n\tindiameter: %.2f", indiameter);
//
//    printf("\n\tedge: \n");
//    for (int j = 0; j < nFaces; ++j) {
//    	printf("\t\t(%d) length = %.2f,\t", j, edge[j].length);
//        printf("centre vectors: [%.2f %.2f]\t", edge[j].vectX, edge[j].vectY);
//    	printf("normal = [%.2f, %.2f]\n", edge[j].normalX, edge[j].normalY);
//    }
//
//    printf("\n\tneighbour's edge: \n");
//    for (int j = 0; j < nFaces; ++j) {
//        if (edge[j].neighEdge != NULL) {
//            printf("\t\t(%d) length = %.2f,\t", j, edge[j].neighEdge->length);
//            printf("centre vectors: [%.2f %.2f]\t", edge[j].neighEdge->vectX, edge[j].neighEdge->vectY);
//            printf("normal = [%.2f, %.2f]\n", edge[j].neighEdge->normalX, edge[j].neighEdge->normalY);
//        }
//    }
//    
//    printf("\tw: [");
//    for (int d = 0; d < dim; ++d) {
//        printf("%.2f ", w[d]);
//    }
//    printf("]\n");
//}
