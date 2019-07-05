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

#include "unsteadyNS_DG.H"
#include "pCorrectEquation.H"
#include "bCorrectEquation.H"

// Constructor
unsteadyNS::unsteadyNS() {}

// Initialization of the case
unsteadyNS::unsteadyNS(int argc, char *argv[])
{

	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
 	#include "createFieldInitialization.H"

    supex = ITHACAutilitiesDG::check_sup();

    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
   tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
    maxIter = ITHACAdict->lookupOrDefault<scalar>("maxIter", 1000);
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

void unsteadyNS::readSolution(fileName casename, label first_snap, label n_snap)
{
    fileName rootpath(".");
    Foam::Time runTime2(Foam::Time::controlDictName, rootpath, casename);
    label last_s;
    if (first_snap >= runTime2.times().size())
    {
        Info << "Error the index of the first snapshot must be smaller than the number of snapshots" << endl;
        exit(0);
    }
    
    _meshRead = autoPtr<dgMesh>
    (
        new dgMesh
        (
            Foam::IOobject
            (
                Foam::dgMesh::defaultRegion,
                runTime2.timeName(),
                runTime2,
                Foam::IOobject::MUST_READ
            )
        )
    );
    Foam::dgMesh& meshRead= _meshRead();
    if (n_snap == 0)
    {
        last_s = runTime2.times().size();
    }
    else
    {
        last_s = min(runTime2.times().size(), n_snap + 1);
    }
    // Read the snapshots for the variables U and p
    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading Snapshot nr " << i - 1 << endl;
        
        Info << "Reading U" << endl;
        dgVectorField U_tmp(
            IOobject
            (
                "U",
                 runTime2.times()[i].name(),
                meshRead,
                IOobject::MUST_READ
            ),
            meshRead
        );
        
        Info << "Reading p" << endl;
        dgScalarField p_tmp(
            IOobject
            (
                "p",
                runTime2.times()[i].name(),
                meshRead,
                IOobject::MUST_READ
            ),
        meshRead
        );
        // Create the U field and P field putting together all the snapshots
        Pfield.append(p_tmp);
        Ufield.append(U_tmp);
        
        Info << "    " << endl;
    }
}


/*------------------------------------------------------------------------------------
    Perform a full-order simulation
--------------------------------------------------------------------------------------- */

void unsteadyNS::truthSolve(List<scalar> mu)
{

    // Create the basic structures
    Time& runTime = _runTime();
    Foam::dgMesh& mesh = _mesh();

    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);

    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;

    //
    // Create  fields (reading from files contained in the directory ./0)
    //
    
    dgScalarField p = _p();
    dgVectorField U = _U();
    #include "createFields.H"

//     p.writeOpt()  = IOobject::AUTO_WRITE;
    U.writeOpt()  = IOobject::AUTO_WRITE;
    p.writeOpt()  = IOobject::AUTO_WRITE;
    
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
//     label outFlag;
//     label inFlag;
//     scalar pi = constant::mathematical::pi;

    scalar paraT1 = 1;
    scalar paraT2 = 1;
    scalar paraT  = 1;
    scalar paraT3 = 1;
    scalar paraT4 = 1;

//     #include "setparaT.H"

    while(runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;
        dimensionedScalar rDeltaT = 1.0/runTime.deltaT();

        //- step 1: velocity prediction
        convecUold = convecU;

        dg::solveEquation(dgm::Sp(convecU) == dgc::div(U, U) );
        
        UT = (a0*U.oldTime() + a1*U.oldTime().oldTime() - (b0*convecU + b1*convecUold) * runTime.deltaT()) / g0;

        UT.correctBoundaryConditions();

        //- step 2: pressure correction
        #include "pEqnCorrect.H"
        shared_ptr<dg::Equation<scalar>> result1 = make_shared<bCorrectEquation<scalar>>(dgm::laplacian(p), paraT1);

        (dg::solveEquation( result1 == (pEqnCorrect + dgc::div(UT) * (rDeltaT * g0 )) ));

        //- step 3: velocity correction
        dg::solveEquation(dgm::Sp(UT) + (dgc::grad(p) * (runTime.deltaT() / g0)) == UT );

        //- step 4: viscosity contribution
        dimensionedScalar tmp( g0/(runTime.deltaT() * nu) );

        shared_ptr<dg::Equation<vector>> result2 = make_shared<bCorrectEquation<vector>>(dgm::laplacian(U), paraT4);

        (dg::solveEquation( dgm::Sp(U) *tmp == result2 + UT*tmp));

        runTime.write();
//         #include "setBoundaryValues.H"

        if(a0<2.0){
            g0 = 1.5, a0 = 2.0, a1 = -0.5, b0 = 2.0, b1 = -1.0;
        }
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
    
            // Write the solution if the target time  is reached
        if (checkWrite(runTime))
        {
            std::ofstream of("./ITHACAoutput/OfflineDG/"+name(counter)+"/"+runTime.timeName());                
                //
                // Evaluation of the primitive variables U and p
//                 forAll(U , ind)
//                 {
//                                 U[ind].x() = rhoU[ind].x()/rho[ind];
//                                 U[ind].y() = rhoU[ind].y()/rho[ind];
//                 }
// 
//                 p = (gamma - dimensionedScalar("one",gamma.dimensions(),1.0))*(Ener - 0.5*(rho*magSqr(U)));
                
                //
                // Write all the output variables
                //
                // State variables
                // Primitive variables 
            exportSolution(p, name(counter), "./ITHACAoutput/OfflineDG");
            exportSolution(U, name(counter), "./ITHACAoutput/OfflineDG");

            Pfield.append(p);
            Ufield.append(U);

            counter++;
            nextWrite += writeEvery;
            writeMu(mu);
        }
    }
}


// Method to solve the supremizer problem
void unsteadyNS::solvesupremizer(word type)
{   
    Info << "SolveSUP" <<endl;
    
   PtrList<dgScalarField> Psup;

   if (type == "modes")
    {
        Info << "Not Supported! " <<endl;
    }
    else if (type == "snapshots")
    {
        // Copy the whole pressure Field: size nPoints x nSnapshots
        for (int i=0; i<Pfield.size(); i++)
        {   
           Psup.append(Pfield[i]);
        }
    }
    else
    {
        std::cout << "You must specify the variable type with either snapshots or modes" << std::endl;
        exit(0);
    }
    // If the supFIeld has already been computed, then read it from file, otherwise  evaluate it
    if (supex == 1)
    {
            readSup("./ITHACAoutput/supfieldDG/");
//         dgVectorField U = _U();
//         dgVectorField Usup
//         (
//             IOobject
//             (
//                 "Usup",
//                 U.time().timeName(),
//                 U.mesh(),
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//                 ),
//             U.mesh(),
//             dimensionedVector("zero", U.dimensions(), vector::zero)
//             );
//         ITHACAstream::read_fields(supfield, Usup, "./ITHACAoutput/supfieldDG/");
    }
    else
    {
       dgVectorField U = _U();

        dgVectorField Usup
        (
            IOobject
            (
                "Usup",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero)
            );

        dimensionedScalar mu_fake
        (
            "mu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
//             dimensionSet(1, -1, -1, 0, 0, 0, 0),
            scalar(1)
            );

        Vector<double> v(0, 0, 0);
        for (label i = 0; i < Usup.boundaryField().size(); i++)
        {
            changeBCtype(Usup, "fixedValue", i);
            assignBC(Usup, i, v);
            assignIF(Usup, v);
        }
        if (type == "snapshots")
        {
            for (label i = 0; i < Psup.size(); i++)
            {
                dg::solveEquation(   -mu_fake*dgm::laplacian( Usup) == dgc::grad(Psup[i]) );
                supfield.append(Usup);
                exportSolution(Usup, name(i + 1), "./ITHACAoutput/supfieldDG/");
            }
            int systemRet = system("ln -s ../../constant ./ITHACAoutput/supfieldDG/constant");
            systemRet += system("ln -s ../../0 ./ITHACAoutput/supfieldDG/0");
            systemRet += system("ln -s ../../system ./ITHACAoutput/supfieldDG/system");
            if (systemRet < 0)
            {
                Info << "System Command Failed in unsteadyNS_DG.C" << endl;
                exit(0);
            }
        }
        else
        {
                exit(0);
        }        
    }
}

void unsteadyNS::readSup(fileName casename, label first_snap, label n_snap)
{
    fileName rootpath(".");
    Foam::Time runTime3(Foam::Time::controlDictName, rootpath, casename);
    label last_s;
    if (first_snap >= runTime3.times().size())
    {
        Info << "Error the index of the first snapshot must be smaller than the number of snapshots" << endl;
        exit(0);
    }
    
    _meshReadSup = autoPtr<dgMesh>
    (
        new dgMesh
        (
            Foam::IOobject
            (
                Foam::dgMesh::defaultRegion,
                runTime3.timeName(),
                runTime3,
                Foam::IOobject::MUST_READ
            )
        )
    );
    Foam::dgMesh& meshReadSup = _meshReadSup();
    
    if (n_snap == 0)
    {
        last_s = runTime3.times().size();
    }
    else
    {
        last_s = min(runTime3.times().size(), n_snap + 1);
    }
    // Read the snapshots for the variables U and p
    for (label i = 2 + first_snap; i < last_s; i++)
    {
        Info << "Reading Snapshot nr " << i - 1 << endl;
        
        Info << "Reading Usup" << endl;
        dgVectorField Usup_tmp(
            IOobject
            (
                "Usup",
                 runTime3.times()[i].name(),
                meshReadSup,
                IOobject::MUST_READ
            ),
            meshReadSup
        );

        // Create the Usup field putting together all the snapshots
        supfield.append(Usup_tmp);
    }
}


/*------------------------------------------------------------------------------------------
    Check if a point in which the solution should be written 
    has been reached
 --------------------------------------------------------------------------------------------*/
bool unsteadyNS::checkWrite(Time& timeObject)
{
    scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));

    scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) - timeObject.deltaTValue());
    if ( diffnow < diffnext)
    {
        return true;
    }
    else
    {
        return false;
    }
}

#include<fstream>

/*------------------------------------------------------------------------------------------
    Write the initialConditions file in ./0
 --------------------------------------------------------------------------------------------*/
void unsteadyNS::writeBCfiles(double rho_bc, double U_bc,  double V_bc )
{
    std::ofstream os("./0/initialConditions", std::ofstream::out);
    os << "U_bc " <<  U_bc  <<  "; \n";
    os<< " V_bc " <<  V_bc  <<  "; \n";
    os << "rhoU_bc " <<  rho_bc*U_bc  <<  "; \n";
    os<< " rhoV_bc " <<  rho_bc*V_bc  <<  "; \n";
    os.close();
}

void unsteadyNS::writeTransportfile(double nu)
{
    std::ofstream os("./constant/initialConditions", std::ofstream::out);
    os << "nu " <<  nu  <<  "; \n";
    os.close();
}
