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

\*---------------------------------------------------------------------------*/

#include "ITHACAstreamDG.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ITHACAstreamDG::read_fields(PtrList<dgVectorField>& Lfield, word Name, fileName casename, label first_snap, label n_snap)
{
    Info << "######### Reading the Data for " << Name << " #########" << endl;
    fileName rootpath(".");
    label last_s;

    Foam::Time runTime2(Foam::Time::controlDictName, rootpath, casename);

    dgMesh mesh
    (
        Foam::IOobject
        (
            Foam::dgMesh::defaultRegion,
            casename + runTime2.timeName(),
            runTime2,
            Foam::IOobject::MUST_READ
        )
    );

    if (first_snap >= runTime2.times().size())
    {
        Info << "Error the index of the first snapshot must be smaller than the number of snapshots" << endl;
        exit(0);
    }

    if (n_snap == 0)
    {
        last_s = runTime2.times().size();
    }
    else
    {
        last_s = min(runTime2.times().size(), n_snap + 2);
    }

    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading " << Name << " number " << i - 1 << endl;
        dgVectorField tmp_field(
            IOobject
            (
                Name,
                runTime2.times()[i].name(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );
        Lfield.append(tmp_field);
    }
    std::cout << std::endl;
}

void ITHACAstreamDG::read_fields(PtrList<dgScalarField>& Lfield, word Name, fileName casename, label first_snap, label n_snap)
{
    Info << " ######### Reading the Data for " << Name << " #########" << endl;
    fileName rootpath(".");
    Foam::Time runTime2(Foam::Time::controlDictName, rootpath, casename);
    label last_s;

    dgMesh mesh
    (
        Foam::IOobject
        (
            Foam::dgMesh::defaultRegion,
            casename + runTime2.timeName(),
            runTime2,
            Foam::IOobject::MUST_READ
        )
    );

    if (first_snap >= runTime2.times().size())
    {
        Info << "Error the index of the first snapshot must be smaller than the number of snapshots" << endl;
        exit(0);
    }

    if (n_snap == 0)
    {
        last_s = runTime2.times().size();
    }
    else
    {
        last_s = min(runTime2.times().size(), n_snap + 2);
    }

    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading " << Name << " number " << i - 1 << endl;

        dgScalarField tmp_field(
            IOobject
            (
                Name,
                runTime2.times()[i].name(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );
        Lfield.append(tmp_field);
    }
    std::cout << std::endl;
}

void ITHACAstreamDG::read_fields(PtrList<dgScalarField>& Lfield, dgScalarField & field, fileName casename, label first_snap, label n_snap)
{
    Info << "######### Reading the Data for " << field.name() << " #########" << endl;
    fileName rootpath(".");
    Foam::Time runTime2(Foam::Time::controlDictName, rootpath, casename);
    label last_s;
    
    if (first_snap >= runTime2.times().size())
    {
        Info << "Error the index of the first snapshot must be smaller than the number of snapshots" << endl;
        exit(0);
    }

    Foam::dgMesh mesh2    (
        Foam::IOobject
        (
            Foam::dgMesh::defaultRegion,
            runTime2.timeName(),
            runTime2,
            Foam::IOobject::MUST_READ
        )
    );
   
    if (n_snap == 0)
    {
        last_s = runTime2.times().size();
    }
    else
    {
        last_s = min(runTime2.times().size(), n_snap + 1);
    }

    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading " << field.name() << " number " << i - 1 << endl;
        dgScalarField tmp_field(
            IOobject
            (
                field.name(),
//                 casename + runTime2.times()[i].name(),
                runTime2.times()[i].name(),
                mesh2,
                IOobject::MUST_READ
            ),
            mesh2
        );

        Lfield.append(tmp_field);
    }
    std::cout << std::endl;
}

void ITHACAstreamDG::read_fields(PtrList<dgVectorField>& Lfield, dgVectorField& field, fileName casename, label first_snap, label n_snap)
{
    Info << "######### Reading the Data for " << field.name() << " #########" << endl;
    fileName rootpath(".");
    Foam::Time runTime2(Foam::Time::controlDictName, rootpath, casename);

    label last_s;
    
    if (first_snap >= runTime2.times().size())
    {
        Info << "Error the index of the first snapshot must be smaller than the number of snapshots" << endl;
        exit(0);
    }
Info << runTime2.times().size() << endl;
    // Create Mesh
    Foam::dgMesh mesh2
    (
        Foam::IOobject
        (
            Foam::dgMesh::defaultRegion,
            runTime2.timeName(),
            runTime2,
            Foam::IOobject::MUST_READ
        )
    );

    if (n_snap == 0)
    {
        last_s = runTime2.times().size();
    }
    else
    {
        last_s = min(runTime2.times().size(), n_snap + 1);
    }

    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading " << field.name() << " number" << i - 1 << endl;
        dgVectorField tmp_field(
            IOobject
            (
                field.name(),
//                 casename + runTime2.times()[i].name(),
                 runTime2.times()[i].name(),
                mesh2,
                IOobject::MUST_READ
            ),
            mesh2
        );
        Lfield.append(tmp_field);
    }
    std::cout << std::endl;
}

// * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //

