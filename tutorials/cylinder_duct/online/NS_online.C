/*------------------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗        ██████╗    ██████╗  
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗       ██╔═══██╗ ██╔════╝  
     ██║   ██║   ███████║███████║██║     ███████║█████╗ ██║   ██║ ██║  ██╗     
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝ ██║   ██║ ██║   ██╗ 
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║       ██████╔═╝ ╚██████╔╝ 
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝       ╚═════╝    ╚═════╝  
                                                                                 
 * In real Time Highly Advanced Computational Applications for Discontinuous Galerkin 
 * Copyright (C) 2018 by the ITHACA-DG authors         
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
    Example of an unsteady Euler Reduction Problem
SourceFiles
    Euler.C
    
\*---------------------------------------------------------------------------*/

#include "reductionProblem.H"
#include "unsteadyNS.H"

#include "ITHACAPOD.H"
#include "reducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include<iomanip>
#include<fstream>


//
// Definition of the class
//
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
    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    
    void offlineSolve()
    {
      // 	Vector<double> inl(0, 0, 0);
      // 	label BCind = 1;
//       	List<scalar> mu_now(1);
        // Check if the offline solution have already been computed
        if (offline)
        {
           // Read the files with the saved snapshots
            ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
            ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
        }
        else
        {
            Info << "NO OFFLINE SOLUTION IS PRESENT" << endl;
            exit(0);
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
    example.finalTime = 5;
    example.timeStep  =  0.0002;;
    example.writeEvery = 0.2;
    
    
    ITHACAparameters para;
    int NmodesUout = para.ITHACAdict->lookupOrDefault<int>("NmodesUout", 12);
    int NmodesPout = para.ITHACAdict->lookupOrDefault<int>("NmodesPout", 12);
    int NmodesSUPout = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPout", 15);

    int NmodesUproj = para.ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);
    int NmodesSUPproj = para.ITHACAdict->lookupOrDefault<int>("NmodesSUPproj", 10);


    // 1 parametrized boundary condition along the x axis(=0)
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 2;
    example.inletIndex(0, 1) = 0;

    // Perform The Offline Solve;
    example.offlineSolve();

    // Solve the supremizer problem
    example.solvesupremizer();

    // Extraction of the modes starting from the primitive variables
    // If the modes have been generated then podex = 1, else podex  = 0
    ITHACAPOD::getModes(example.Ufield, example.Umodes, example.podex, 0, 0, NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example.podex, 0, 0, NmodesPout);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex, example.supex, 1, NmodesSUPout);

    Info << example.podex <<endl;

    example.projectPPE("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
//     example.projectPPE("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    
    reducedUnsteadyNS reduced(example);

    // Set values of the reduced stuff
    reduced.nu = 0.0052;
    reduced.tstart = 0.0;
    reduced.finalTime =  3.5;
    reduced.dt = 0.0002;
    
    // Set the online velocity
    Eigen::MatrixXd vel_now(1,1);
    vel_now(0, 0) =  0.165;
    
    reduced.solveOnline_ch(vel_now);

    // Reconstruct the solution and export it
    reduced.reconstruct_ch("./ITHACAoutput/Reconstruction_CH", 2500);


    Info<< "Final Time Reached\n" << endl;

    return 0;
}


