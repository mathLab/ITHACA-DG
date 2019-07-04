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
    reductionProblem

Description
    A general reduction problem class is implemented

SourceFiles
    reductionProblem.C
    reductionProblemTemplates.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reductionProblem class.
/// \dir
/// Directory containing the header, source and template files for the reductionProblem class.


#ifndef reductionProblem_H
#define reductionProblem_H

#include <random>
#include "dgCFD.H"
#include "IOmanip.H"
#include "OFstream.H"
#include <sys/stat.h>
#include "ITHACAutilitiesDG.H"
// #include "ITHACAstreamDG.H"
// #include "../../ITHACAstream/ITHACAparameters.H"
#include "../../thirdparty/Eigen/Eigen/Eigen"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// A general class for the implementation of a full order parametrized problem
class reductionProblem
{
private:

public:
    // Constructors

    /// Construct Null
    reductionProblem();
    virtual ~reductionProblem() {};

    // Members
    /// Number of parameters
    label Pnumber;
    /// Dimension of the training set (used only when gerating parameters without input)
    label Tnumber;
    /// Matrix of parameters
    Eigen::MatrixXd mu;
    /// Range of the parameter spaces
    Eigen::MatrixXd mu_range;
    /// Current value of the parameter
//     double mu_cur;
    /// Boolean variable, it is 1 if the POD has already been computed, else 0
    bool podex;
    /// Boolean variable, it is 1 if the Offline phase has already been computed, else 0
    bool offline;
    /// dictionary to store input output infos
   IOdictionary* ITHACAdict;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcomment"
    /// Matrix that contains informations about the inlet boundaries
    /// @details The dimension is: <br>
    /// Rows = Number of parametrized boundary conditions <br>
    /// Cols = 2 <br>
    /// Example: <br>
    /// example.inletIndex.resize(1, 2); <br>
    /// example.inletIndex(0, 0) = 0; <br>
    /// example.inletIndex(0, 1) = 0; <br>
    /// Means that there is one parametrized boundary conditions, the first col contains the
    /// index of the parametrized boundary, the second col contains the direction along which
    /// the BC is parametrized: <br>
    /// 0 for x <br>
    /// 1 for y <br>
    /// 2 for z <br>
    /// The Matrix should be implemented in the following way:
    /// \\[
    /// \\mbox{inletIndex}=\\begin{bmatrix}
    /// \\mbox{BC1_index} & \\mbox{Direction of BC1} \\\\ \
    /// \\mbox{BC2_index} & \\mbox{Direction of BC3} \\\\ \
    /// \\vdots & \\vdots \\\\ \
    /// \\mbox{BCN_index} & \\mbox{Direction of BCN} \\ \
    /// \\end{bmatrix}
    /// \\]
    Eigen::MatrixXi inletIndex;
#pragma GCC diagnostic pop

    // Functions
     
    //--------------------------------------------------------------------------
    /// Set Parameters Problems
    void setParameters();

    //--------------------------------------------------------------------------
    /// Generate Random Numbers
//     void genRandPar();

    //--------------------------------------------------------------------------
    /// Generate Random Numbers given the dimension of the training set
    ///
    /// @param[in]  tsize  Dimension of the training set.
    ///
//     void genRandPar(int tsize);

    //--------------------------------------------------------------------------
    /// Generate Equidistributed Numbers
    void genEquiPar();

    //--------------------------------------------------------------------------
    /// Perform a TruthSolve
    void truthSolve();

    
    Eigen::MatrixXi inletIndexT;

    //--------------------------------------------------------------------------
    /// Export a field (template function)
    ///
    /// @param[in]  s          field you want to export.
    /// @param[in]  subfolder  The subfolder used to save the snapshot.
    /// @param[in]  folder     The folder used to save the snapshot.
    ///
    /// @tparam     T          Type of parameter you want to export (dgVectorField or dgScalarField).
    ///
    template<typename T>
    void exportSolution(T& s, fileName subfolder, fileName folder = "./Offline");

    
    //--------------------------------------------------------------------------
    /// Assign Boundary Condition to a volVectorField
    ///
    /// @param[in]  s       field where you want to assign the BC.
    /// @param[in]  BC_ind  The BC index.
    /// @param[in]  value   The value you want to assign (it must be an OpenFOAM vector).
    ///
    void assignBC(dgVectorField& s, label BC_ind, Vector<double>& value);

    //--------------------------------------------------------------------------
    /// Assign Boundary Condition to a volScalarField
    ///
    /// @param[in]  s       field where you want to assign the BC.
    /// @param[in]  BC_ind  The BC index.
    /// @param[in]  value   The value you want to assign (it must be a double).
    ///
    void assignBC(dgScalarField& s, label BC_ind, double& value);

    //--------------------------------------------------------------------------
    /// Change the boundary condition type for a volVectorField
    ///
    /// @param[in]  field   The field
    /// @param[in]  BCtype  The BC type you want to set
    /// @param[in]  BC_ind  The BC index
    ///
    void changeBCtype(dgVectorField& field, word BCtype, label BC_ind);

    //--------------------------------------------------------------------------
    /// Change the boundary condition type for a volScalarField
    ///
    /// @param[in]  field   The field
    /// @param[in]  BCtype  The BC type you want to set
    /// @param[in]  BC_ind  The BC index
    void changeBCtype(dgScalarField& field, word BCtype, label BC_ind);
    //--------------------------------------------------------------------------
    /// Assign internal field condition
    ///
    /// @param[in,out]  s      field where you want to assign the internal field condition
    /// @param          value  The value you want to assign
    ///
    /// @tparam         T      type of field (volVectorField or volScalarField)
    /// @tparam         G      type of value you want to assign (OpenFOAM vector or scalar)
    ///
    template<typename T, typename G>
    void assignIF(T& s, G& value);

    //--------------------------------------------------------------------------
    /// Write out a list of scalar
    ///
    /// @param[in]  mu_now  The list of scalars.
    ///
    void writeMu(List<scalar> mu_now);
    
};


#ifdef NoRepository
#   include "reductionProblemTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// Additional Doxygen Documentation
/// \dir ../problems
/// Directory containing the implementation of the full order parametrized problems
/// and of the projection methods used to create the reduced matrices.

// Doxigen Documentation (Main Page)

/*! \mainpage In real Time Highly Advanced Computational Applications for Finite Volumes
 *
 * \section About About ITHACA-DG
 *
 * ITHACA-DG (In real Time Highly Advanced Computational Applications for
Finite Volumes) is C++ library based on the finite volume solver
[<b> OpenFOAM </b>](https://openfoam.org/). It consists of the implementation of several reduced order
modeling techniques for parametrized problems. ITHACA-DG is developed and maintained by
[<b> Dr. Andrea Lario </b>] under the supervision of [<b> Prof. Gianluigi Rozza </b>](http://people.sissa.it/~grozza/) at [<b> SISSA mathLab </b>](https://mathlab.sissa.it/) an applied math laboratory of the
[<b> International School for advanced studies </b>](https://www.sissa.it/).

 */






