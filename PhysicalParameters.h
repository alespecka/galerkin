/* 
 * File:   PhysicalParameters.h
 * Author: ales
 *
 * Created on May 28, 2015, 4:16 PM
 */

#pragma once

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

struct PhysicalParameters {
    static const int nEqns = 4;  // number of components of Euler equations in 2D
    const double kapa;
    const double reynolds;
    const double prandtl;
    
    bool isFlowViscous;

    // inlet boundary condition:
    const double pIn0;
    const double angleIn;
    const double rhoIn0;
    // outlet boundary condition:
    double pOut0;
    
    PhysicalParameters(double kapa, double reynolds, double prandtl,
        double pIn0, double angleIn, double rhoIn0, double pOut0) :
        
        kapa(kapa), reynolds(reynolds), prandtl (prandtl), pIn0(pIn0),
        angleIn(angleIn), rhoIn0(rhoIn0), pOut0(pOut0)
    {
    
        if (reynolds > 1e-2) {
            isFlowViscous = true;
        }
        else {
            isFlowViscous = false;
        }
    }
    
    PhysicalParameters(boost::property_tree::ptree paramTree) :
        kapa(paramTree.get<double>("kapa")),
        reynolds(paramTree.get<double>("reynolds")),
        prandtl(paramTree.get<double>("prandtl")),
        
        pIn0(paramTree.get<double>("pIn")),
        rhoIn0(paramTree.get<double>("rhoIn")),
        angleIn(paramTree.get<double>("angleIn"))                 
    {        
        boost::optional<double> pOutOptional = paramTree.get_optional<double>("pOut");
        boost::optional<double> machOptional = paramTree.get_optional<double>("mach");
        
        if (pOutOptional && machOptional) {
            throw boost::property_tree::ptree_error("only one of the nodes pOut or mach must by defined");
        }
        
        if (pOutOptional) {
            pOut0 = pOutOptional.get();
        } else if (machOptional) {
            double mach0 = machOptional.get();
            pOut0 = pIn0 / std::pow((kapa-1.)/2.*mach0*mach0 + 1., kapa/(kapa-1.));
        } else {
            throw boost::property_tree::ptree_error("either pOut or mach must by defined");
        }
                
        if (reynolds > 1e-2) {
            isFlowViscous = true;
        } else {
            isFlowViscous = false;
        }
    }
};
