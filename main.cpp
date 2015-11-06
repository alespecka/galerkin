/* 
 * File:   main.cpp
 * Author: ales
 *
 * Created on 12. listopadu 2014, 13:10
 */

//#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <ctime>
#include <omp.h>
#include <sys/time.h>
#include <assert.h>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <chrono> // high resolution clock

#include "element.h"
#include "printarray.h"
#include "list"

#include "GaussQuad.h"
#include "GaussQuadTri.h"
#include "sparseMatrix.h"
#include "NumericalParameters.h"
#include "PhysicalParameters.h"

#include "umfsolver.h"

using namespace std;
using namespace std::chrono;

typedef std::vector<double> DVector;
typedef std::vector<Element *> ElemVector;

const std::string quadPath = "gauss/";
const std::string quadTriPath = "gauss_tri/";

const std::string pointFileName = "points.txt";
const std::string triangleFileName = "triangles.txt";
const std::string neighbourFileName = "neighbours.txt";

/**
 * Gets rid of spaces before and after the string str.
 * 
 * @param str string with potential spaces at the beginning or at the end
 * @param whitespace
 * @return trimmed string
 */
std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

void loadParameters(const std::string &parameterFileName, std::shared_ptr<PhysicalParameters>& physParam,
    std::shared_ptr<NumericalParameters>& numParam, std::string &meshPath, std::string &outFile)
{
    boost::property_tree::ptree tree;
    boost::property_tree::ptree childTree;
    boost::property_tree::read_xml(parameterFileName, tree);
    
    childTree = tree.get_child("parameters.physicalParameters");     
    physParam.reset(new PhysicalParameters(childTree));
    
    childTree = tree.get_child("parameters.numericalParameters");    
    numParam.reset(new NumericalParameters(childTree));    
        
    childTree = tree.get_child("parameters"); 
    
    // load matrices containing the data of the mesh
    meshPath = trim(childTree.get<std::string>("meshPath"));
    
    outFile = meshPath + "/";
    outFile += trim(childTree.get<std::string>("outputFilePrefix"));
    
    // copy the file with parameters
    std::ifstream src(parameterFileName, std::ios::binary);
    std::ofstream dest(outFile + "-parameteres.xml", std::ios::binary);
    dest << src.rdbuf();
}

int loadMeshMatrices(const std::string &meshPath, double *&points, 
                     int *&triangles, int *&neighbours) {               
          
    
    int nTriangles, arrayLength; //, nPoints;
        
    std::string pointFilePath =     meshPath + "/" + pointFileName;
    std::string triangleFilePath =  meshPath + "/" + triangleFileName;
    std::string neighbourFilePath = meshPath + "/" + neighbourFileName;
    
    arrayLength = readArray(points, pointFilePath);      
    arrayLength = readArray(triangles, triangleFilePath);
    readArray(neighbours, neighbourFilePath);
//    nPoints = arrayLength / 2 ; 
    nTriangles = arrayLength / 3;    
       
    return nTriangles;
}

void loadQuadratures(int orderOfOccuracy, std::shared_ptr<GaussQuad>& pQuad, 
                     std::shared_ptr<GaussQuadTri>& pQuadTri)
{
    int polynomialOrder = 2 * (orderOfOccuracy - 1);
        
    std::string quadFileName = quadPath + std::to_string(orderOfOccuracy) + ".txt";
    pQuad.reset(new GaussQuad(quadFileName));
    
    std::string quadTriFileName = quadTriPath + std::to_string(std::max(polynomialOrder, 1)) + ".txt";
    pQuadTri.reset(new GaussQuadTri(quadTriFileName));
}

/*
 * print result to text file
 */ 
void saveResults(arma::cube &W, ElemVector &mesh, const std::string &outFile) {   
    std::string outFileW =              outFile + ".txt";
    std::string outFileCentreSolution = outFile + "-centres.txt";
    std::string outFileViscosity      = outFile + "-artificialViscosity.txt";
    
    std::ofstream fout;
    int dim = Element::nEqns;
    int nTriangles = mesh.size();
    int nBasis = W.size() / (dim * nTriangles);
    Element::Vector4 w;
        
    fout.open(outFileCentreSolution.c_str());
    for (int i = 0; i < mesh.size(); ++i) {
        mesh[i]->getW(w, mesh[i]->centreX, mesh[i]->centreY, W.slice(i));
        for (int d = 0; d < dim; ++d) {
            fout << w[d] << " ";
        }
        fout << std::endl;
    }
    fout.close();

//    fout.open(outFileNodeSolution);
//    for (int i = 0; i < mesh.size(); ++i) {
//        for (int j = 0; j < Element::nFaces; ++j) {
//            mesh[i]->combine(w, mesh[i]->nodeX[j], mesh[i]->nodeY[j], &W[i * dim * nBasis]);
//            for (int d = 0; d < dim; ++d) {
//                fout << w[d] << " ";
//            }
//        }
//        fout << endl;
//    }
//    fout.close();

    fout.open(outFileW.c_str());
    for (int i = 0; i < nTriangles; ++i) {
        for (int j = 0; j < dim * nBasis; ++j) {
            fout << W[i*dim*nBasis + j] << " ";
        }
        fout << "\n";
    }
    fout.close();
    
    fout.open(outFileViscosity);
    for (int i = 0; i < mesh.size(); ++i) {
        fout << mesh[i]->amountOfViscosity << std::endl;
    }
    fout.close();
}

void initialConditions(double *w, PhysicalParameters param) {
    double kapa = param.kapa;
    double machIn = sqrt(2 / (kapa - 1) * (pow(param.pIn0 / param.pOut0, ((kapa - 1) / kapa)) - 1));
    double rhoIn = param.rhoIn0 * pow(1 + (kapa - 1) / 2 * machIn * machIn, 1 / (1 - kapa));
    double velocity = abs(machIn * sqrt(kapa * param.pOut0 / rhoIn));

    double u = velocity * cos(param.angleIn);
    double v = velocity * sin(param.angleIn);

    double E = param.pOut0 / (kapa - 1) + 0.5 * rhoIn * velocity * velocity;

    w[0] = rhoIn;
    w[1] = rhoIn * u;
    w[2] = rhoIn * v;
    w[3] = E;
}

void printCurrentState(int stepIndex, double time, double timeStep,
                       double maxResidue, struct timeval start) {
    struct timeval end;
    
    gettimeofday(&end, NULL);
    double currentTime = ((end.tv_sec  - start.tv_sec) * 1000000u + 
                  end.tv_usec - start.tv_usec) / 1.e6;
    printf("%d) ", stepIndex);
    printf("residuum: %.2e   ", maxResidue);
    printf("physical time: %.2f   ", time + timeStep);
    printf("computational time: %dm %2ds   ", (int) currentTime / 60, (int) currentTime % 60);
//            printf("cfl: %.2f   ", cfl);
    printf("step: %.5f\n", timeStep);    
}

void printCurrentState(int stepIndex, double time, double timeStep, double maxResidue,
                      struct timeval start, int nElements, NumericalParameters &numParam) {
    
    if (stepIndex % numParam.printFrequancy != 0) {
        printCurrentState(stepIndex, time, timeStep, maxResidue, start);
    }
    
    std::cout << "order of occuracy: " << numParam.orderOfOccuracy << std::endl;
    std::cout << "number of elements: " << nElements << std::endl;
    #ifdef _OPENMP
        std::cout << "number of threads: " << numParam.nThreads << std::endl;
    #endif    
}

void solvePDE(arma::cube& W, DVector& timeMesh, ElemVector& mesh, NumericalParameters& numParam) {    
    
    const int nEqns = Element::nEqns;
    const int nBasis = numParam.nBasis;
    double time = 0;
    double timeStep;        

    struct timeval start;  // std::clock_t startTime = std::clock();
    gettimeofday(&start, NULL);

    SparseMatrix *pJacobi;
    if (numParam.implicit) {
        int size = std::pow(nBasis * nEqns, 2) * ((Element::nFaces + 1) * mesh.size() - numParam.nBoundaryEdges);
        pJacobi = new SparseMatrix(size);
    }
    SparseMatrix &jacobi = *pJacobi;

    arma::cube newW(nBasis, nEqns, mesh.size());
    arma::cube RW(nBasis, nEqns, mesh.size());
    arma::mat centreW(nEqns, mesh.size());    
    arma::vec tempVect(nEqns);
    for (int i = 0; i < mesh.size(); ++i) {
        mesh[i]->getW(tempVect, mesh[i]->centreX, mesh[i]->centreY, W.slice(i));
        centreW.col(i) = tempVect;
    }   

    double residue[nEqns];
    double maxResidue = numParam.terminationCondition;
    double cfl;
    int k = 0;
    int n = 10;    
    while (k < timeMesh.size() && (k < numParam.minSteps || maxResidue >= numParam.terminationCondition)) {        
        if (k < n && numParam.implicit && !numParam.continueCalculation) {
            cfl = (k + 1) * numParam.maxCfl / n;
        }
        else {
            cfl = numParam.maxCfl;
        }
        
        /* calculate time step according to CFL condition */
        double timeStep = 10000; // random sufficiently large number
        for (int i = 0; i < mesh.size(); ++i) {
            double step = cfl / (2 * numParam.orderOfOccuracy - 1) * mesh[i]->CFL(&centreW[i * nEqns]);
            if (step < timeStep) {
                timeStep = step;
            }
        }
        
        /* compute the amount of artificial viscosity that needs to be added */
        #pragma omp parallel for
        for (int i = 0; i < mesh.size(); ++i) {
            mesh[i]->shockIndicator(W.slice(i));
        }
                
        if (numParam.implicit) {  /*** Backward Euler ***/            
//            high_resolution_clock::time_point start = high_resolution_clock::now();
            
            jacobi.position = 0;
            #pragma omp parallel for
            for (int i = 0; i < mesh.size(); ++i) {
                mesh[i]->jacobi(jacobi, RW.slice(i), W, mesh, timeStep);
            }            
            
            #ifdef UMFSOLVER
                umfSolve(jacobi.size, W.size(), &jacobi.rowIdx[0], &jacobi.colIdx[0], &jacobi.values[0], &RW[0], &newW[0]);
            #else
                std::cerr << "UMFPACK has not been activated" << std::endl;
                exit(1);
            #endif /* UMFSOLVER */                        
            
            for (int i = 0; i < W.size(); ++i) {
                newW[i] = W[i] + newW[i];                                
            }
                
//            high_resolution_clock::time_point end = high_resolution_clock::now();
//            duration<double> duration = end - start;            
//            std::cout << "duration: " << duration.count() << std::endl;
//            for (int i = 0; i < mesh.size(); ++i) {
//                mesh[i]->artificialDamping(newW, quadTri, numParam, centreX, centreY);
//            }
        }
        else {  /*** second order Runge-Kutta ***/            
            #pragma omp parallel for
            for (int i = 0; i < mesh.size(); ++i) {
                mesh[i]->residualVector(RW.slice(i), W, mesh);
                
                for (int pos = i*nEqns*nBasis; pos < (i+1)*nEqns*nBasis; ++pos) {
                    newW[pos] = W[pos] + 0.5 * timeStep * RW[pos];                    
                }                
//                  mesh[i]->artificialDamping(newW, quadTri);
            }           
            
            #pragma omp parallel for
            for (int i = 0; i < mesh.size(); ++i) {
                mesh[i]->residualVector(RW.slice(i), newW, mesh);
                
                for (int pos = i*nEqns*nBasis; pos < (i+1)*nEqns*nBasis; ++pos) {
                    newW[pos] = W[pos] + timeStep * RW[pos];
                }
//                std::cout << "W: " << std::endl << W.slice(i) << std::endl;
//                std::cout << "RW: " << std::endl << RW.slice(i) << std::endl;
//                std::cout << "newW: " << std::endl << newW.slice(i) << std::endl;
//                std::cout << "newW2: " << std::endl << (W.slice(i) + timeStep * RW.slice(i)) << std::endl;
//                mesh[i]->artificialDamping(newW, quadTri);
            }
        }
        
        /*** check whether the solution contains any NANs ***/
        for (int i = 0; i < mesh.size(); ++i) {   
            for (int pos = i*nEqns*nBasis; pos < (i+1)*nEqns*nBasis; ++pos) {
                if (newW[pos] != newW[pos]) { 
                    std::cerr << "NAN - time step: " << k
                              << ", triangle: " << i << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        
//        for (int i = 0; i < mesh.size(); ++i) {            
//            if (W[i * dim * nBasis] < numParam.tol) {
//                W[i * dim * nBasis] = numParam.tol;
//                std::cerr << "pressure is too low" << std::endl;
//            }
//        }

        /*** calculate residuum ***/
        Element::Vector4 neww;
        std::fill(residue, residue + nEqns, 0);
        for (int i = 0; i < mesh.size(); ++i) {
            mesh[i]->getW(neww, mesh[i]->centreX, mesh[i]->centreY, newW.slice(i));
            for (int d = 0; d < nEqns; ++d) {
                residue[d] += mesh[i]->area * std::abs(neww[d] - centreW[i * nEqns + d]) / numParam.domainArea;
                centreW[i * nEqns + d] = neww[d];
            }
        }
        maxResidue = *std::max_element(residue, residue + nEqns);

        for (int i = 0; i < W.size(); ++i) {
            W[i] = newW[i];
        }
        
//        std::cout << "W: " << std::endl << W << std::endl;

        /*** print some information about the progress of the computation ***/
        if ((k+1) % numParam.printFrequancy == 0 || k == 0) {
            printCurrentState(k+1, time, timeStep, maxResidue, start);
        }

        /*** take a time step ***/
        time += timeStep;
        timeMesh[k] = time;
        ++k;
    }
        
    printCurrentState(k, time, timeStep, maxResidue, start, mesh.size(), numParam);    
    timeMesh.resize(k - 1);    
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "exactly one parameter is required" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string paramFileName = argv[1];
    
    try {
        std::shared_ptr<PhysicalParameters> physParam;
        std::shared_ptr<NumericalParameters> numParam;
        std::string meshPath, outFile;
        loadParameters(paramFileName, physParam, numParam, meshPath, outFile);
        
        std::shared_ptr<GaussQuad> quad;
        std::shared_ptr<GaussQuadTri> quadTri;
        loadQuadratures(numParam->orderOfOccuracy, quad, quadTri);
                
        #ifdef _OPENMP
            omp_set_num_threads(numParam->nThreads); // specify the number of threads
            omp_lock_t lock;                        
            Element::initialiseClass(physParam, numParam, quad, quadTri, &lock);
        #else  /* _OPENMP */          
            Element::initialiseClass(physParam, numParam, quad, quadTri);
        #endif /* _OPENMP */

        /*
         * read data from text files
         */
        double *points;
        int *triangles;
        int *neighbours;        
        int nTriangles = loadMeshMatrices(meshPath, points, triangles, neighbours);        

        int nBasis = numParam->nBasis;
        int nEqns = PhysicalParameters::nEqns;
        arma::cube W(nBasis, nEqns, nTriangles); // solution vector

        /*
         * set initial conditions
         */ 
        if (numParam->continueCalculation == true) {
            readArray(&W[0], nBasis * nEqns * nTriangles, outFile + ".txt");
        }
        else {    
            double w0[nEqns];
            initialConditions(w0, *physParam);
            for (int i = 0; i < nTriangles; ++i) {
                for (int d = 0; d < nEqns; ++d) {
                    for (int j = 0; j < nBasis; ++j) {
                        W(j, d, i) = w0[d];
                    }
                }
            }
        }

        /*
         * create mesh
         */
        ElemVector mesh;
        mesh.reserve(nTriangles);
        for (int i = 0; i < nTriangles; ++i) {
            Element *elem = new Element(i, points, triangles, neighbours, nBasis);
            mesh.push_back(elem);
        }

        /*
         * solve problem
         */
        DVector timeMesh(numParam->nTimeSteps); // declare time mesh
        solvePDE(W, timeMesh, mesh, *numParam);                

        /*
         * print result to a text file
         */
        saveResults(W, mesh, outFile);

        for (Element *elem : mesh)
            delete elem;        
        
        delete points;
        delete triangles;
        delete neighbours;

        return EXIT_SUCCESS;
    
    } catch (boost::property_tree::ptree_error &ex) {
        std::cerr << "file " << paramFileName << " has a wrong format: "<< ex.what() << std::endl;
        exit(EXIT_FAILURE);
    } catch (std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

