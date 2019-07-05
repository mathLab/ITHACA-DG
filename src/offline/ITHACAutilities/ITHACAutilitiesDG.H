/*------------------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗        ██████╗    ██████╗  
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗       ██╔═══██╗ ██╔════╝  
     ██║   ██║   ███████║███████║██║     ███████║█████╗ ██║   ██║ ██║  ██╗     
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝ ██║   ██║ ██║   ██╗ 
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║       ██████╔═╝ ╚██████╔╝ 
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝       ╚═════╝    ╚═════╝  
                                                                                 
 * In real Time Highly Advanced Computational Applications for Discontinuous Galerkin 
 * Copyright (C) 2019 by the ITHACA-DG authors         
--------------------------------------------------------------------------------------*\
                                                       
License                                                
    This file is part of ITHACA-DG

    ITHACA-DG is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ITHACA-DG is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-DG. If not, see <http://www.gnu.org/licenses/>.

Class
    ITHACAutilitiesDG

Description
    Utilities to manage files and check existence of folders

SourceFiles
    ITHACAutilitiesDG.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAutilitiesDG class.
/// \dir
/// Directory containing the header and source files for the ITHACAutilitiesDG class.

#ifndef ITHACAutilitiesDG_H
#define ITHACAutilitiesDG_H

 #include "dgCFD.H"
#include "IOmanip.H"
// #include "freestreamFvPatchField.H"
#include <sys/stat.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "../../thirdparty/Eigen/Eigen/Eigen"
#include "../../thirdparty/Eigen/Eigen/LU"
#pragma GCC diagnostic pop
#include <functional>
#include "./colormod.H"

using namespace std::placeholders;

// #include "Foam2Eigen.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class ITHACAutilitiesDG Declaration
\*---------------------------------------------------------------------------*/

/// Class to implement some useful functions inside ITHACA-FV
class ITHACAutilitiesDG
{
private:

public:

    //--------------------------------------------------------------------------
    /// Check if a folder exists
    ///
    /// @param[in]  folder  The folder you want to check.
    ///
    /// @return     1 if the folder exists 0 elsewhere.
    ///
    static bool check_folder(word folder);

    //--------------------------------------------------------------------------
    /// Check if the offline data folder exists
    ///
    /// @return     1 if the folder exists 0 elsewhere.
    ///
    static bool check_off();

    //--------------------------------------------------------------------------
    /// Check if the POD data folder exists
    ///
    /// @return     1 if the folder exists 0 elsewhere.
    ///
    static bool check_pod();

    //--------------------------------------------------------------------------
    /// Check if the supremizer folder exists
    ///
    /// @return     1 if the folder exists 0 elsewhere.
    ///
    static bool check_sup();


    //--------------------------------------------------------------------------
    /// Assign Boundary Condition to a dgScalarField
    ///
    /// @param[in]  s       field where you want to assign the BC.
    /// @param[in]  BC_ind  The BC index.
    /// @param[in]  value   The value you want to assign (it must be a double).
    ///
    static void assignBC(dgScalarField& s, label BC_ind, double& value);

    //--------------------------------------------------------------------------
    /// Assign Boundary Condition to a dgVectorField
    ///
    /// @param[in]  s       field where you want to assign the BC.
    /// @param[in]  BC_ind  The BC index.
    /// @param[in]  value   The value you want to assign (it must be an OpenFOAM vector).
    ///
    static void assignBC(dgVectorField& s, label BC_ind, Eigen::MatrixXd value);

    //--------------------------------------------------------------------------
    /// Assign Boundary Condition to a dgScalarField
    ///
    /// @param[in]  s       field where you want to assign the BC.
    /// @param[in]  BC_ind  The BC index.
    /// @param[in]  value   The value you want to assign (it must be a double).
    ///
    static void assignBC(dgScalarField& s, label BC_ind, Eigen::MatrixXd value);



#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcomment"
#pragma GCC diagnostic pop
};

#endif






