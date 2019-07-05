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
    ITHACAstreamDG

Description
    Methods for input and outputs from and to files

SourceFiles
    ITHACAstreamDG.C
    ITHACAstreamDGTemplates.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAstreamDG class, it contains the implementation of
/// several methods for input output operations.
/// \dir
/// Directory containing the header, source and template files for the ITHACAstreamDG class, ITHACAparameters and Assert Functions

#ifndef ITHACAstreamDG_H
#define ITHACAstreamDG_H

// #include "fvCFD.H"
#include "dgCFD.H"
#include "OFstream.H"
#include "IOmanip.H"
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <algorithm>
#include <fstream>
#include "ITHACAassert.H"
// #include "ITHACAutilitiesDG.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "../../thirdparty/Eigen/Eigen/Eigen"
#pragma GCC diagnostic pop

#define MAXBUFSIZE (static_cast<int> (1e6))


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class for input-output manipulation
class ITHACAstreamDG
{
private:

public:


    //--------------------------------------------------------------------------
    /// Funtion to read a list of dgVectorField from name of the field and casename
    ///
    /// @param[in]  Lfield      a PtrList of dgVectorField where you want to store the field.
    /// @param[in]  Name        The name of the field you want to read.
    /// @param[in]  caseName    The folder where the field is stored.
    /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
    /// @param[in]  n_snap      The number of snapshots you want to read.
    ///
    static void read_fields(PtrList<dgVectorField>& Lfield, word Name, fileName caseName, label first_snap = 0, label n_snap = 0);

    //--------------------------------------------------------------------------
    /// Funtion to read a list of dgScalarField from name of the field and folder name
    ///
    /// @param[in]  Lfield      a PtrList of dgScalarField where you want to store the field.
    /// @param[in]  Name        The name of the field you want to read.
    /// @param[in]  caseName    The folder where the field is stored.
    /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
    /// @param[in]  n_snap      The number of snapshots you want to read.
    ///
    static void read_fields(PtrList<dgScalarField>& Lfield, word Name, fileName caseName, label first_snap = 0, label n_snap = 0);

    //--------------------------------------------------------------------------
    /// Funtion to read a list of dgVectorField from name of the field a dgVectorField if it is already existing
    ///
    /// @param[in]  Lfield      a PtrList of dgVectorField where you want to store the field.
    /// @param[in]  field       The field used as template to read other fields.
    /// @param[in]  caseName    The folder where the field is stored.
    /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
    /// @param[in]  n_snap      The number of snapshots you want to read.
    ///
    static void read_fields(PtrList<dgVectorField>& Lfield, dgVectorField& field, fileName caseName, label first_snap = 0, label n_snap = 0);

    //--------------------------------------------------------------------------
    /// Funtion to read a list of dgScalarField from casename and a list of dgScalarField with the same name
    ///
    /// @param[in]  Lfield      a PtrList of dgScalarField where you want to store the field.
    /// @param[in]  field       The field used as template to read other fields.
    /// @param[in]  caseName    The folder where the field is stored.
    /// @param[in]  first_snap  The first snapshots from which you want to start reading the field.
    /// @param[in]  n_snap      The number of snapshots you want to read.
    static void read_fields(PtrList<dgScalarField>& Lfield, dgScalarField& field, fileName caseName, label first_snap = 0, label n_snap = 0);


    //--------------------------------------------------------------------------
    /// Export a generic field to file in a certain folder and subfolder
    /* This is a function to export a generic field into a certain subfolder */
    ///
    /// @param[in] s volVectorField or volScalarField.
    /// @param[in] subfolder string to indicated the subfolder where the field is stored.
    /// @param[in] folder string to indicated the folder where the field is stored.
    /// @param[in] fieldName string to indicate the name of the field/file.
    ///
    template<typename T>
    static void exportSolution(T& s, fileName subfolder, fileName folder, word fieldName);

};

template<typename T>
void ITHACAstreamDG::exportSolution(T& s, fileName subfolder, fileName folder, word fieldName)
{
    mkDir(folder + "/" + subfolder);
    T act(fieldName, s);
    fileName fieldname = folder + "/" + subfolder + "/" + fieldName;
    OFstream os(fieldname);
    act.writeHeader(os);
    os << act << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






