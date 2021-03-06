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
    foam2eigen

Description
    Converts foam objects into Eigen and viceversa

SourceFiles
    foam2eigen.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the Foam2Eigen class.
/// \dir
/// Directory containing the header and source files for the Foam2Eigen class.

#ifndef Foam2Eigen_H
#define Foam2Eigen_H

#include "fvCFD.H"
#include "IOmanip.H"
#include <tuple>
#include <sys/stat.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "../../thirdparty/Eigen/Eigen/Eigen"
#include "../../thirdparty/Eigen/Eigen/LU"
#pragma GCC diagnostic pop
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
  Class foam2eigen Declaration
\*---------------------------------------------------------------------------*/

/// Class to converts OpenFOAM objects into Eigen and viceversa
class Foam2Eigen
{
private:

public:

    //--------------------------------------------------------------------------
    /// @brief      Convert a ldu OpenFOAM matrix into a Eigen Matrix A and a source vector b
    ///
    /// @param[in]      foam_matrix       The foam matrix can be fvScalarMatrix or fvVectorMatrix
    /// @param[out]     A                 The matrix or the list of matrices, can be dense or sparse
    /// @param[out]     b                 The source term vector or the list of source term vectors, always dense
    ///
    /// @tparam         type_f            The type of field, can be a volvectorfield or a volscalarfield
    /// @tparam         type_foam_matrix  The type of foam matrix can be fvScalarMatrix or fvVectorMatrix
    /// @tparam         type_A            The type of matrix or the list of matrices, can be dense (Eigen::MatrixXd) or sparse (Eigen::SparseMatrix<double>) or a List of them
    /// @tparam         type_B            The type source term vector or the list can be Eigen::VectorXd or List<Eigen::VectorXd>
    ///
    template <class type_foam_matrix, class type_A, class type_B>
    static void fvMatrix2Eigen(type_foam_matrix& foam_matrix, type_A& A, type_B& b);


    //--------------------------------------------------------------------------
    /// @brief      Convert a PtrList of snapshots to Eigen matrix (only internal field)
    ///
    /// @param[in]  fields The fields can be a PtrList<volScalarField> or PtrList<volVectorField>
    ///
    /// @tparam     type_f  Type of the fields can be PtrList<volScalarField> or PtrList<volVectorField>
    ///
    /// @return     An Eigen matrix containing as columns the snapshots
    ///
    ///
    ///
    template <class type_f>
    static Eigen::MatrixXd PtrList2Eigen(type_f& fields, int Nfields = 10000000);


    //--------------------------------------------------------------------------
    /// @brief      Convert an OpenFOAM field to an Eigen Vector
    ///
    /// @param[in]  field   The field
    ///
    /// @tparam     type_f  Type of the field can be volScalarField or volVectorField
    ///
    /// @return     return the dense Eigen vector
    ///
    template <class type_f>
    static Eigen::VectorXd field2Eigen(type_f& field);

    //--------------------------------------------------------------------------
    /// @brief      Convert an OpenFOAM field to a List of Eigen Vectors, one for each boundary
    ///
    /// @param[in]  field   The field
    ///
    /// @tparam     type_f  Type of the field can be volScalarField or volVectorField
    ///
    /// @return     return the dense Eigen vector
    ///
    template <class type_f>
    static List<Eigen::VectorXd> field2EigenBC(type_f& field);


    //--------------------------------------------------------------------------
    /// @brief      Converts an OpenFOAM PtrList to a list of Eigen Matrices, one for each boundary
    ///
    /// @param      fields   The fields
    /// @param[in]  Nfields  The Number of requested modes.
    ///
    /// @tparam     type_f   type of the field, can be volVectorField or volScalarField
    ///
    /// @return     List of Eigen Matrices
    ///
    template <class type_f>
    static List<Eigen::MatrixXd> PtrList2EigenBC(type_f& fields, int Nfields = 10000000);

    //--------------------------------------------------------------------------
    /// @brief      Convert an Eigen Vector to an OpenFOAM field (only the internal Field is copied inside the OpenFOAM field)
    ///
    /// @param[in,out]  field         The OpenFOAM field where to copy the Eigen Vector
    /// @param[in]      eigen_vector  The Eigen vector
    ///
    /// @tparam     type_f         Type of the field can be volScalarField or volVectorField
    ///
    /// @return     The OpenFoam Field
    ///
    template <class type_f>
    static type_f Eigen2field(type_f& field, Eigen::VectorXd& eigen_vector);

    //--------------------------------------------------------------------------
    /// @brief      Perform the projection of an OpenFOAM fvMatrix onto a set of modes using the Eigen matrix multiplication
    ///
    /// It is performed the standard projection used in reduced order models to obtain the reduced matrices:
    ///
    /// \f[  \mathbf{A} \in R^{NxN} , \mathbf{A_r} \in R^{nxn}  , \mathbf{b} \in R^{Nx1} , \mathbf{b_r} \in R^{nx1} , \mathbf{M} \in R^{Nxn} \f]
    ///
    /// \f[ \mathbf{A_r} = \mathbf{M^T A M} , \mathbf{b_r} = \mathbf{M^T b}  \f]
    ///
    /// @param[in]  matrix        The A OpenFOAM fvMatrix
    /// @param[in]  modes         The M OpenFOAM modes contained in a PtrList
    /// @param[in]  Nmodes        The n number of modes you want to use (The size of input modes can be bigger than the requested number of modes)
    ///
    /// @tparam     type_m        type of fvMatrix can be either scalar or vectorial
    /// @tparam     type_PtrList  the type of the PtrList for the modes can be volScalarField or volVectorField
    ///
    /// @return     A tuple containing the projected Ar matrix and vector br which are coming from the source term.
    ///
    template <class type_m, class type_PtrList>
    static std::tuple<Eigen::MatrixXd, Eigen::VectorXd> projectFvMatrix(type_m& matrix, type_PtrList& modes, int Nmodes);

    //--------------------------------------------------------------------------
    /// @brief      Perform the projection of an OpenFOAM field onto a set of modes using the Eigen matrix multiplication
    ///
    /// It is the standard Galerkin projection of a field onto a reduced basis:
    ///
    /// \f[  \mathbf{f} \in R^{Nx1} , \mathbf{f_r} \in R^{nx1} , \mathbf{M} \in R^{Nxn} \f]
    ///
    /// \f[ \mathbf{f_r} = \mathbf{M^T f}  \f]
    ///
    /// @param[in]  field         The f OpenFOAM field
    /// @param[in]  modes         The M OpenFOAM modes contained in a PtrList
    /// @param[in]  Nmodes        The n number of modes you want to use (The size of input modes can be bigger than the requested number of modes)
    ///
    /// @tparam     type_f        type of field can be either scalar or vectorial
    /// @tparam     type_PtrList  the type of the PtrList for the modes can be volScalarField or volVectorField
    ///
    /// @return     a vector fr with the projection
    ///
    template <class type_f, class type_PtrList>
    static Eigen::VectorXd projectField(type_f& field, type_PtrList& modes, int Nmodes);

    //--------------------------------------------------------------------------
    /// @brief      Obtain the Mass Matrix from a list of reduced basis
    ///
    /// @param      modes         The OpenFOAM modes contained in a PtrList
    /// @param[in]  Nmodes        The number of requested modes
    ///
    /// @tparam     type_PtrList  the type of the PtrList for the modes can be volScalarField or volVectorField
    ///
    /// @return     An Eigen matrix containing the mass matrix
    ///
    template <class type_PtrList>
    static Eigen::MatrixXd MassMatrix(type_PtrList& modes, int Nmodes);

    //--------------------------------------------------------------------------
    /// @brief      Convert a PtrList of OpenFOAM fvMatrix into a tuple of lists of Eigen Sparse Matrices and source vectors
    ///
    /// @param[in]      MatrixList   The matrix list
    ///
    /// @tparam[in]     type_matrix  type of fvMatrix can be either scalar or vectorial
    ///
    /// @return     a tuple containing a list of sparse matrices and a list of source vectors
    ///
    template<class type_matrix>
    static std::tuple<List<Eigen::SparseMatrix<double> >, List<Eigen::VectorXd> > LFvMatrix2LSM(PtrList<type_matrix>& MatrixList);


};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(volVectorField& field);

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(volScalarField& field);

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(fvMesh const& field);

template<>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(volVectorField& field);

template<>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(volScalarField& field);

template<class type_f>
List<Eigen::MatrixXd> PtrList2EigenBC(type_f& fields, int Nfields);

template<>
List<Eigen::MatrixXd> PtrList2EigenBC(PtrList<volScalarField>& fields, int Nfields);

template<>
List<Eigen::MatrixXd> PtrList2EigenBC(PtrList<volVectorField>& fields, int Nfields);


template<class type_f>
Eigen::VectorXd Foam2Eigen::field2Eigen(type_f& field)
{
    Eigen::VectorXd out;
    out.resize(int(field.size()));
    for (int l = 0; l < field.size(); l++)
    {
        out(l) = field[l];
    }
    return out;
}

template<class type_matrix>
std::tuple<List<Eigen::SparseMatrix<double> >, List<Eigen::VectorXd> > Foam2Eigen::LFvMatrix2LSM(PtrList<type_matrix>& MatrixList)
{
    List<Eigen::SparseMatrix<double> > SM_list;
    List<Eigen::VectorXd> V_list;
    int LSize =  MatrixList.size();
    SM_list.resize(LSize);
    V_list.resize(LSize);
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    for (int j = 0; j < LSize; j++)
    {
        fvMatrix2Eigen(MatrixList[j], A, b);
        SM_list[j] = A;
        V_list[j] = b;
    }
    std::tuple <List<Eigen::SparseMatrix<double> >, List<Eigen::VectorXd> > tupla;
    tupla = std::make_tuple(SM_list, V_list);
    return tupla;
}

template <class type_m, class type_PtrList>
std::tuple<Eigen::MatrixXd, Eigen::VectorXd> Foam2Eigen::projectFvMatrix(type_m& matrix, type_PtrList& modes, int Nmodes)
{
    Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd Ar;
    Eigen::VectorXd b;
    Eigen::VectorXd br;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Foam2Eigen::fvMatrix2Eigen(matrix, A, b);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Volumes.rows(), Nmodes);
    if (Volumes.rows() != Eig_Modes.rows())
    {
        VolumesN.resize(Eig_Modes.rows(), Nmodes);
    }
    if (Volumes.rows() == Eig_Modes.rows())
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i) = Volumes;
        }
    }
    else
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i).segment(0, Volumes.rows()) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() + 1, Volumes.rows() * 2) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() * 2 + 1, Volumes.rows() * 3) = Volumes;
        }
    }
    Ar = Eig_Modes.transpose() * A * Eig_Modes;
    br = Eig_Modes.transpose() * b;
    std::tuple <Eigen::MatrixXd, Eigen::VectorXd> tupla;
    tupla = std::make_tuple(Ar, br);
    return tupla;
}

template <class type_PtrList>
Eigen::MatrixXd Foam2Eigen::MassMatrix(type_PtrList& modes, int Nmodes)
{
    Eigen::MatrixXd Mr;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Volumes.rows(), Nmodes);
    if (Volumes.rows() != Eig_Modes.rows())
    {
        VolumesN.resize(Eig_Modes.rows(), Nmodes);
    }
    if (Volumes.rows() == Eig_Modes.rows())
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i) = Volumes;
        }
    }
    else
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i).segment(0, Volumes.rows()) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() + 1, Volumes.rows() * 2) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() * 2 + 1, Volumes.rows() * 3) = Volumes;

        }
    }
    Mr = Eig_Modes.transpose() * (Eig_Modes.cwiseProduct(VolumesN));
    return Mr;
}

template <class type_f, class type_PtrList>
Eigen::VectorXd Foam2Eigen::projectField(type_f& field, type_PtrList& modes, int Nmodes)
{
    Eigen::VectorXd fr;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Eigen::VectorXd f = Foam2Eigen::field2Eigen(field);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Volumes.rows(), 1);
    VolumesN = Volumes;
    if (Volumes.rows() != Eig_Modes.rows())
    {
        VolumesN.resize(Eig_Modes.rows(), 1);
        VolumesN.col(0).segment(0, Volumes.rows()) = Volumes;
        VolumesN.col(0).segment(Volumes.rows() + 1, Volumes.rows() * 2) = Volumes;
        VolumesN.col(0).segment(Volumes.rows() * 2 + 1, Volumes.rows() * 3) = Volumes;
    }
    fr = Eig_Modes.transpose() * (f.cwiseProduct(VolumesN));
    return fr;
}
#endif


