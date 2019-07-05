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
    elemReduction

SourceFiles
    elemReduction.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the elemReduction class.
/// \dir
/// Directory containing the header and source files for the elemReduction class.

#ifndef elemReduction_H
#define elemReduction_H
#include <iostream>
#include "IOmanip.H"
#include <stdio.h>
#include "OFstream.H"
#include "dgCFD.H"
#include <sys/stat.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
// #include "GenEigsSolver.h"
// #include <SymEigsSolver.h>
#pragma GCC diagnostic pop

#include "../../thirdparty/Eigen/Eigen/Eigen"
#include "../../thirdparty/Eigen/Eigen/LU"
// #include <functional>



/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class of the computation of the POD, it exploits the method of snapshots.
class elemReduction
{

public:

    // Public Members

    //--------------------------------------------------------------------------
    /// Compute the reduced field for a vector field
    ///
    /// @param[in]  snapshotsU  a PtrList of volVectorField snapshots.
    /// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
    ///
    static void reduce(PtrList<dgVectorField>& snapshotsU, bool sup = 0);

    //--------------------------------------------------------------------------
    /// Compute the reduced field  for a scalar field
    ///
    /// @param[in]  snapshotsP  a PtrList of volScalarField snapshots.
    /// @param[in]  sup         boolean variable 1 if you want to compute the supremizer modes 0 elsewhere.
    ///
    static void reduce(PtrList<dgScalarField>& snapshotsP, bool sup = 0);

    //--------------------------------------------------------------------------
    /// Export the basis for a vector field into the ITHACAOutput/POD or ITHACAOutput/supremizer
    ///
    /// @param[in] s        a PtrList of volVectorField where the snapshots associated with the bases are stored.
    /// @param[in] bases    a PtrList of volVectorField where the bases are stored.
    /// @param[in] sup      a boolean variable 1 if you want to export the supremizer bases (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
    ///
    static void exportBases(PtrList<dgVectorField>& s, PtrList<dgVectorField>& bases, bool sup = 0);

    //--------------------------------------------------------------------------
    /// Export the basis for a scalar field into the ITHACAOutput/ folder
    ///
    /// @param[in] s        a PtrList of volScalarField where the snapshots associated with the bases are stored.
    /// @param[in] bases    a PtrList of volScalarField where the bases are stored.
    /// @param[in] sup      a boolean variable 1 if you want to export the supremizer bases (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
    ///
    static void exportBases(PtrList<dgScalarField>& s, PtrList<dgScalarField>& bases, bool sup = 0);

    //--------------------------------------------------------------------------
    /// Export the reduced solution for a vector field into the ITHACAOutput/ folder
    ///
    /// @param[in] s        a PtrList of volVectorField where the snapshots  are stored.
    /// @param[in] name   name of the folder.
    /// @param[in] sup      a boolean variable 1 if you want to export the supremizer bases (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
    ///
    static void exportReducedSolution(dgVectorField& s, fileName name, bool sup = 0);

    //--------------------------------------------------------------------------
    /// Export the reduced solution for a scalar field into the ITHACAOutput/ folder
    ///
    /// @param[in] s        a PtrList of volScalarField where the snapshots  are stored.
    /// @param[in] name   name of the folder.
    /// @param[in] sup      a boolean variable 1 if you want to export the supremizer bases (in ITHACAOutput/supremizer) 0 elsewhere (Default is 0).
    ///
    static void exportReducedSolution(dgScalarField& s, fileName name, bool sup = 0);
    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif



