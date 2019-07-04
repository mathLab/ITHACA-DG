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

#include "unsteadyEuler.H"

// Constructor
unsteadyEuler::unsteadyEuler() {}

// Initialization of the case
unsteadyEuler::unsteadyEuler(int argc, char *argv[])
{

	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
 	#include "createFields.H"

// #include "createFvOptions.H"
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


/*------------------------------------------------------------------------------------
    Perform a full-order simulation
--------------------------------------------------------------------------------------- */

void unsteadyEuler::truthSolve(List<scalar> mu)
{
//         #include "setRootCase.H"

//      #include "createTime.H"
//      #include "createMesh.H"     


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
    //     #include "createTimeControls.H"    
    #include "createFields.H"

    
    //
    // Definition of variables 
    //
    
    // Conservative variables
    //     dgScalarField rho = _rho();
    //     dgVectorField rhoU = _rhoU();
    //     dgVectorField Ener = _Ener;
    

    // Auxiliary variables
    dgScalarField rho1("rho1", rho);
    dgVectorField rhoU1("rhoU1", rhoU);
    dgScalarField Ener1("Ener1", Ener); 
    
    
    //
    // Perform the algorithm for the solution of the 
    //  unsteady  Euler equations
    //

     while(runTime.run())
    {
        runTime++;

        // 2-order SSP Runge-kutta

        //SSP RK stage 1.
        rho1 = rho;
        rhoU1 = rhoU;
        Ener1 = Ener;
//         setBoundaryValues(rho1, rhoU1, Ener1, gamma, runTime - runTime.deltaT());

        rho1.storeOldTime();
        rhoU1.storeOldTime();
        Ener1.storeOldTime();
        rho1.updateGaussField();
        rhoU1.updateGaussField();
        Ener1.updateGaussField();

        gther_U = rhoU1.gaussField()/rho1.gaussField();
        gther_p = (gamma - dimensionedScalar("one",gamma.dimensions(),1.0))*(Ener1.gaussField() - 0.5*(rho1.gaussField()*magSqr(gther_U)));

        Godunov.update(gther_U, gther_p, gamma.value(), rho1.gaussField(), rhoU1.gaussField(), Ener1.gaussField());

        dg::solveEquation(dgm::ddt(rho1) + dgc::div(gther_U, rho1, Godunov.fluxRho()));

        dg::solveEquation(dgm::ddt(rhoU1) + dgc::div(gther_U, rhoU1, Godunov.fluxRhoU()) + dgc::grad(gther_p));

        dg::solveEquation(dgm::ddt(Ener1) + dgc::div(gther_U, Ener1, Godunov.fluxEner())+ dgc::div(gther_U, gther_p));

        Godunov.limite(rho1,rhoU1,Ener1);

        rho1.storeOldTime();
        rhoU1.storeOldTime();
        Ener1.storeOldTime();

//         SSP RK Stage 2.
        rho1.updateGaussField();
        rhoU1.updateGaussField();
        Ener1.updateGaussField();

        gther_U = rhoU1.gaussField()/rho1.gaussField();
        gther_p = (gamma - dimensionedScalar("one",gamma.dimensions(),1.0))*(Ener1.gaussField() - 0.5*rho1.gaussField()*magSqr(gther_U));

        Godunov.update(gther_U, gther_p, gamma.value(), rho1.gaussField(), rhoU1.gaussField(), Ener1.gaussField());

        dg::solveEquation(dgm::ddt(rho1) + dgc::div(gther_U, rho1, Godunov.fluxRho()));

        dg::solveEquation(dgm::ddt(rhoU1) + dgc::div(gther_U, rhoU1, Godunov.fluxRhoU()) + dgc::grad(gther_p));

        dg::solveEquation(dgm::ddt(Ener1) + dgc::div(gther_U, Ener1, Godunov.fluxEner())+ dgc::div(gther_U, gther_p));

         rho = 0.5*rho      + 0.5*rho1;
         rhoU = 0.5*rhoU + 0.5*rhoU1;
         Ener = 0.5*Ener  + 0.5*Ener1;

       Godunov.limite(rho,rhoU,Ener);
        
        rho.correctBoundaryConditions();
        rhoU.correctBoundaryConditions();
        Ener.correctBoundaryConditions();

        runTime.write();

        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
            
            // Write the solution if the target time  is reached
             if (checkWrite(runTime))
             {
                std::ofstream of("./ITHACAoutput/Offline/"+name(counter)+"/"+runTime.timeName());                
        Info<< "Time xxx= " << runTime.timeName() << nl << endl;
                //
                // Evaluation of the primitive variables U and p
                forAll(U , ind)
                {
                                U[ind].x() = rhoU[ind].x()/rho[ind];
                                U[ind].y() = rhoU[ind].y()/rho[ind];
                }

                p = (gamma - dimensionedScalar("one",gamma.dimensions(),1.0))*(Ener - 0.5*(rho*magSqr(U)));
                
                //
                // Write all the output variables
                //
                // State variables
                exportSolution(rhoU, name(counter), "./ITHACAoutput/Offline");
                exportSolution(rho, name(counter), "./ITHACAoutput/Offline");
                exportSolution(Ener, name(counter), "./ITHACAoutput/Offline");

                // Primitive variables 
                exportSolution(p, name(counter), "./ITHACAoutput/Offline");
                exportSolution(U, name(counter), "./ITHACAoutput/Offline");

                rhoU_field.append(rhoU);
                rho_field.append(rho);
                Ener_field.append(Ener);
                p_field.append(p);
                U_field.append(U);
                
                counter++;
                nextWrite += writeEvery;
                writeMu(mu);
             }
    }
}


// Method to solve the supremizer problem
void unsteadyEuler::solvesupremizer(word type)
{   
    
    
    Info << "SolveSUP" <<endl;
    
   PtrList<dgScalarField> P_sup;
Info << "HERE: write type " << type <<endl;
   if (type == "modes")
    {
        Info << "Not Supported! " <<endl;
    }
    else if (type == "snapshots")
    {
        // Copy the whole pressure Field: size nPoints x nSnapshots
        for (int i=0; i<p_field.size(); i++)
        {   
           P_sup.append(p_field[i]);
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
        ITHACAstreamDG::read_fields(sup_field, Usup, "./ITHACAoutput/supfield/");
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
//             dimensionSet(0, 2, -1, 0, 0, 0, 0),
            dimensionSet(1, -1, -1, 0, 0, 0, 0),
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
            for (label i = 0; i < P_sup.size(); i++)
            {

                dg::solveEquation(   -mu_fake*dgm::laplacian( Usup) == dgc::grad(P_sup[i]) );
            
            // ORIGINAL  FV-CODE:obsolete
//                 dgVectorMatrix u_sup_eqn
//                 (
//                     - fvm::laplacian(nu_fake, Usup)
//                     );
//                 solve
//                 (
//                     u_sup_eqn == fvc::grad(P_sup[i])
//                     );
                sup_field.append(Usup);
                exportSolution(Usup, name(i + 1), "./ITHACAoutput/supfield/");
            }
            int systemRet = system("ln -s ../../constant ./ITHACAoutput/supfield/constant");
            systemRet += system("ln -s ../../0 ./ITHACAoutput/supfield/0");
            systemRet += system("ln -s ../../system ./ITHACAoutput/supfield/system");
            if (systemRet < 0)
            {
                Info << "System Command Failed in steadyNS.C" << endl;
                exit(0);
            }
        }
        else
        {
              // NOT YET SUPPORTED!
                exit(0);
//             for (label i = 0; i < Pmodes.size(); i++)
//             {
// 
//                 fvVectorMatrix u_sup_eqn
//                 (
//                     - fvm::laplacian(nu_fake, Usup)
//                     );
//                 solve
//                 (
//                     u_sup_eqn == fvc::grad(Pmodes[i])
//                     );
//                 supmodes.append(Usup);
//                 exportSolution(Usup, name(i+1), "./ITHACAoutput/supremizer/");
//             }
//             int systemRet = system("ln -s ../../constant ./ITHACAoutput/supremizer/constant");
//             systemRet += system("ln -s ../../0 ./ITHACAoutput/supremizer/0");
//             systemRet += system("ln -s ../../system ./ITHACAoutput/supremizer/system");
//             if (systemRet < 0)
//             {
//                 Info << "System Command Failed in steadyNS.C" << endl;
//                 exit(0);
//             }
        }
//         
    }
}




/*------------------------------------------------------------------------------------------
    Compute the lifting function
 --------------------------------------------------------------------------------------------*/
void unsteadyEuler::liftSolve()
{
    for (label k = 0; k < inletIndex.rows(); k++)
    {
    }
}



/*------------------------------------------------------------------------------------------
    Check if a point in which the solution should be written 
    has been reached
 --------------------------------------------------------------------------------------------*/
bool unsteadyEuler::checkWrite(Time& timeObject)
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
void unsteadyEuler::writeBCfiles(double rho_bc, double U_bc,  double V_bc )
{
    std::ofstream os("./0/initialConditions", std::ofstream::out);
    os << "U_bc " <<  U_bc  <<  "; \n";
    os<< " V_bc " <<  V_bc  <<  "; \n";
    os << "rhoU_bc " <<  rho_bc*U_bc  <<  "; \n";
    os<< " rhoV_bc " <<  rho_bc*V_bc  <<  "; \n";
    os.close();
}
