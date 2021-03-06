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
    ITHACAparameters

Description
    Some parameters for input output informations

SourceFiles
    IHACAparameters.H

\*---------------------------------------------------------------------------*/

#ifndef Parameters_H
#define Parameters_H

#include <iostream>
#include "dgCFD.H"

/// Class for the definition of some general parameters, the parameters must be defined from the file ITHACAdict inside the
/// system folder.
class ITHACAparameters
{
public:
    ITHACAparameters()
    {
        Foam::Time runTime(Foam::Time::controlDictName, "./", "");
        Foam::dgMesh mesh
        (
            Foam::IOobject
            (
                Foam::dgMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        );
        ITHACAdict = new IOdictionary
        (
            IOobject
            (
                "ITHACAdict",
                "./system",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
       precision = ITHACAdict->lookupOrDefault<int>("OutPrecision", 10);
       word typeout = ITHACAdict->lookupOrDefault<word>("OutType", "fixed");
       if (typeout == "fixed") {outytpe = std::ios_base::fixed;}
       if (typeout == "scientific") {outytpe = std::ios_base::scientific;}
//     eigensolver = ITHACAdict->lookupOrDefault<word>("EigenSolver", "spectra");
        eigensolver = ITHACAdict->lookupOrDefault<word>("EigenSolver", "spectra");

    }
    /// type of eigensolver used in the eigenvalue decomposition can be either be eigen or spectra
    word eigensolver;

    /// precision of the output Market Matrix objects (i.e. reduced matrices, eigenvalues, ...)
    int precision;

    /// type of output format can be fixed or scientific
    std::_Ios_Fmtflags outytpe;

    /// Dictionary for input objects from file
   IOdictionary* ITHACAdict;

};

#endif
