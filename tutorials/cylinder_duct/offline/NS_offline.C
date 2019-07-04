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

Description
    Example of an unsteady Reduction Problem
    
\*---------------------------------------------------------------------------*/

#include "reductionProblem.H"
#include "unsteadyNS_DG.H"
#include "elemReduction.H"
#include "ITHACAstreamDG.H"
#include<iomanip>
#include<fstream>

//
// Definition of the class
class tutorialNS: public unsteadyNS
{
 public:
 
    //Constructor
    explicit tutorialNS(int argc, char *argv[])
    :
    unsteadyNS(argc, argv),
    p(_p()),
    U(_U())
    {};

    // Definition of fields
    dgScalarField& p;
    dgVectorField& U;

    void offlineSolve()
    {
      // 	Vector<double> inl(0, 0, 0);
      // 	label BCind = 1;
      	List<scalar> mu_now(1);
        // Check if the offline solution have already been computed
        if (offline)
        {
           // Read the files with the saved snapshots
            Info<< "------------------------------------------------------------ " << endl;
            Info<< "                   READING FIELDS " <<endl;
            Info<< "------------------------------------------------------------ "  << endl;

            readSolution("./ITHACAoutput/OfflineDG/");    
        }
        else
        {
            // Compute the solutions  of the truth problem
            for (label i = 0; i < mu.rows(); i++)
            {
                for (label j=0;  j < mu.cols();  j++)
                {
                    mu_now[0] = mu(i, j);

                    writeTransportfile(mu(i, j));

                    Info<< "------------------------------------------------------------ " << endl;
                    Info<< "                   SOLVING FOR mu = " << mu( i , j ) <<endl;
                    Info<< "------------------------------------------------------------ "  << endl;
                    // Change value of the parameter and call the truthSolver
                    truthSolve( mu_now );
                }
            }
        }
    }
};


/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

	// Construct the tutorialEuler object
    tutorialNS example(argc, argv);


    // Read parameters from ITHACAdict file
//     ITHACAparameters para;

     /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 15;

    /// Set the parameters infos
    example.setParameters();

    // Set the parameter ranges
    // The first index is associated to a parameter, the second one to the values
    example.mu_range(0, 0) = 0.001;
    example.mu_range(0, 1) = 0.01;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();

   
    // Time parameters
    // startTime: time at which the first snapshot is saved
    // writeEvery: time Interval between two snapshots
    // finalTime: final time of the simulation.
    example.startTime = 0;
    example.finalTime  = 5;
    example.timeStep   =  0.0002;;
    example.writeEvery = 0.2;

    //  Info<< example.mu(0,2) << endl;

    // Perform The Offline Solve
    example.offlineSolve();

    Info << "-----------------------------" << endl;
    Info << "            podex "                 << example.podex << endl;
    Info << "-----------------------------" << endl;

    // Extraction of the modes starting from the primitive variables
    // If the modes have been generated then podex = 1, else podex  = 0
    elemReduction::reduce(example.Pfield);
    elemReduction::reduce(example.Ufield);

    Info<< "Reached final time!" << endl;

    exit 0;
}



/// \dir 04unsteadyNS Folder of the turorial 4
/// \file 
/// \brief Implementation of tutorial 4 for an unsteady Navier-Stokes problem

/// \example 01unsteadyNS.C
/// \section intro_unsreadyNS Introduction to tutorial 4
/// In this tutorial we implement a parametrized unsteady Navier-Stokes 2D problem where the parameter is the kinematic viscosity.
/// The physical problem represents an incompressible flow passing around a very long cylinder. The simulation domain is rectangular
/// with spatial bounds of [-0.2, 2], and [-0.2, 0.21] in the X and Y directions, respectively. The cylinder has a radius of
/// 0.05 unit length and its center is located in the origin of the reference system.  On the upper and lower bounds, no slip boundary conditions are applied.
///The system has a prescribed uniform inlet velocity of 1 m/s which is constant through the whole simulation time.
///
///
/// \section code01 A detailed look into the code
/// 
/// \subsection header ITHACA-DG header files
///
/// First of all let's have a look at the header files that need to be included and what they are responsible for.
/// 
/// The header files of ITHACA-DG necessary for this tutorial are: <unsteadyNS.H> for the full order unsteady NS problem, <elemReduction.H> for the decrease of Computational points,
/// and finally <ITHACAstream.H> for some ITHACA input-output operations.
/// 
/// \dontinclude 04unsteadyNS.C
/// \skip unsteadyNS
/// \until ITHACAstream
/// 
/// \subsection classtutorial01 Definition of the tutorial01 class
/// 
/// We define the tutorial04 class as a child of the unsteadyNS class.
/// The constructor is defined with members that are the fields need to be manipulated
/// during the resolution of the full order problem using pimpleFoam. Such fields are
/// also initialized with the same initial conditions in the solver.
/// \skipline example
/// \until {}
/// 
/// Inside the example class we define the offlineSolve method according to the
/// specific parametrized problem that needs to be solved. 
///It loops over all the parameters, changes the system viscosity with the iterable parameter
/// then performs the offline solve.
/// 
/// \skipline offlineSolve
/// \until }
/// \skipline else
/// \until }
/// \skipline }
/// \skipline }
/// 
/// We note that in the commented line we show that it is possible to parametrize the boundary conditions.
/// For further details we refer to the classes: reductionProblem, and unsteadyNS.
///
/// \subsection main Definition of the main function
/// 
/// In this section we show the definition of the main function.
/// First we construct the object "example" of type tutorial01:
/// 
/// \skipline example
///
/// Now we would like to perform 15 parametrized simulations where the kinematic viscosity
/// is the sole parameter to change, and it lies in the range of {0.01, 0.001} m^2/s equispaced.
/// The inlet velocity and the domain geometry are both kept fixed through all
/// simulations.
///
/// In our implementation, the parameter (viscosity) can be defined by specifying that
/// Nparameters=1, Nsamples=15, and the parameter ranges from 0.01 to 0.001 equispaced, i.e.
///
/// \skipline example.Pnumber
/// \until example.genEquiPar()
///
/// After that we set the inlet boundaries where we have the non homogeneous BC:
///
/// \skipline example.inlet
/// \until example.inletIndex(0, 1) = 0;
///
/// And we set the parameters for the time integration, so as to simulate 20 seconds for each
/// simulation, with a step size = 0.0002 seconds, and the data are dumped every 1.0 seconds, i.e.
///
/// \skipline example.startTime
/// \until example.writeEvery
///
/// Now we are ready to perform the offline stage:
///
/// \skipline Solve()
