/* 
 * File:   NumericalParameters.h
 * Author: ales
 *
 * Created on May 28, 2015, 4:19 PM
 */

#pragma once

#include <iostream>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

struct NumericalParameters {
    const bool continueCalculation;
    const bool implicit;
    
    const int printFrequancy;
    
    const double terminationCondition;
    const int nTimeSteps;
    const int minSteps;
    const double maxCfl;
    const int orderOfOccuracy;
    int nBoundaryEdges;
    double domainArea;
    
    const double dampMax;
    const double dampBegin;
    const double dampEnd;
    
    const double tol;
    const double diffTol;  // differentiation tolerance
    const double penalty;
    const int nThreads;
    
    const int nBasis;   // number of basis functions
    int **powers;
    
    NumericalParameters(bool continueCalculation, bool implicit, int printFrequancy, double terminationCondition,
        int nTimeSteps, int minSteps, double maxCfl, double tol, double diffTol, int orderOfOccuracy, double dampMin,
        double dampTol, double dampConst, double penalty, int nThreads) :
        
        continueCalculation(continueCalculation), terminationCondition(terminationCondition), nTimeSteps(nTimeSteps),
        minSteps(minSteps), maxCfl(maxCfl), tol(tol), diffTol(diffTol), orderOfOccuracy(orderOfOccuracy),
        dampMax(dampMin), dampBegin(dampTol), dampEnd(dampConst), implicit(implicit), printFrequancy(printFrequancy),
        penalty(penalty), nThreads(nThreads), nBasis(0.5 * orderOfOccuracy * (orderOfOccuracy + 1)), nBoundaryEdges(0),
        domainArea(0.0)
    {
        initialise();
    }
    
    NumericalParameters(boost::property_tree::ptree params) :
        continueCalculation(params.get<bool>("continue")), 
    
        terminationCondition(params.get<double>("terminationCondition")),
        orderOfOccuracy(params.get<int>("orderOfOccuracy")),
        nThreads(params.get<int>("numThreads")),
        
        implicit(params.get<bool>("implicit")),
        maxCfl(params.get<double>("cfl")),
        nTimeSteps(params.get<int>("maxNumSteps")),
        minSteps(params.get<int>("minNumSteps")),
        printFrequancy(params.get<int>("printFrequancy")),

        tol(params.get<double>("tol")),
        diffTol(params.get<double>("diffTol")),
        penalty(params.get<double>("penalty")),

        dampMax(params.get<double>("damping.max")),
        dampBegin(params.get<double>("damping.begin")),    
        dampEnd(params.get<double>("damping.end")),
    
        nBasis(0.5 * orderOfOccuracy * (orderOfOccuracy + 1)),
        nBoundaryEdges(0),
        domainArea(0.0)
    {    
        initialise();
    }
    
    ~NumericalParameters() {
        delete[] powers[0];
        delete[] powers[1];
        delete[] powers;
    }
    
private:
    void initialise() {        
        powers = new int*[2];
        powers[0] = new int[nBasis];
        powers[1] = new int[nBasis];
        
        int n = 0;
        for (int degree = 0; degree < orderOfOccuracy; ++degree) {
            for (int j = 0; j <= degree; ++j) {
                powers[0][n] = degree - j;
                powers[1][n] = j;
                ++n;
            }
        }        
    }
};

