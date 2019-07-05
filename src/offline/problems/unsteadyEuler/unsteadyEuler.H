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
    unsteadyEuler

Description
    Reduction class for a unsteady Euler problem.

SourceFiles
    unsteadyEuler.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the steadyNS class.
/// \dir
/// Directory containing the header and source files for the steadyNS class.

#ifndef unsteadyEuler_H
#define unsteadyEuler_H
#include "reductionProblem.H"
#include "dgCFD.H"
#include "ITHACAstreamDG.H"
#include "dgFields.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNS Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> steady NS problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a steady NS problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class unsteadyEuler: public reductionProblem
{
	public:

	// Constructors    
	/// Null constructor
	unsteadyEuler();

	/// Construct with argc and argv
	unsteadyEuler(int argc, char *argv[]);
	
	//Destructor
	virtual ~unsteadyEuler(){};


    /// Counter used for the output of the snapshots of the Offline solution
    label counter=1;
    
    /// List of pointers used to form the density snapshots matrix
    PtrList<dgScalarField> rho_field;

    /// List of pointers used to form the velocity snapshots matrix
    PtrList<dgVectorField> rhoU_field;

    /// List of pointers used to form the Energy snapshots matrix
    PtrList<dgScalarField> Ener_field;

    /// List of pointers used to form the velocity snapshots matrix
    PtrList<dgVectorField> U_field;

    /// List of pointers used to form the pressure snapshots matrix
    PtrList<dgScalarField> p_field;
    
    /// List of pointers used to form the supremizer snapshots matrix
    PtrList<dgVectorField> sup_field;

    /// List of pointers used for the density modes
    PtrList<dgVectorField> rho_modes;
    
    /// List of pointers used to form the pressure modes
    PtrList<dgScalarField> p_modes;

    /// List of pointers used to form the velocity modes
    PtrList<dgVectorField> U_modes;

    
	/// List of pointers used to form the supremizer modes
	PtrList<dgVectorField> sup_modes;

	/// List of pointers used to form the list of lifting functions
	PtrList<dgVectorField> liftfield;
    
/*
	/// List of pointers used to form the homogeneous velocity snapshots
	PtrList<volVectorField> Uomfield;
*/
	/// Tolerance for the residual of the stationary problems, there is the same tolerance for velocity and pressure
	scalar tolerance;

	/// Number of maximum iterations to be done for the computation of the truth solution
	scalar maxIter;

	// Specific variable for the unstationary case
	/// Start Time (initial time to start storing the snapshots)
	scalar startTime;

	/// Final time (final time of the simulation and consequently of the acquisition of the snapshots)
	scalar finalTime;

        /// Time step of the simulation
	scalar timeStep;

        /// Time step of the writing procedure
	scalar writeEvery = timeStep;

        /// Auxiliary variable to store the next writing instant
	scalar nextWrite;

	/// adjustTimeStep
	bool adjustTimeStep;

	/// maxCourant
	scalar maxCo;

	/// maxDeltaT
	scalar maxDeltaT;
    
	// Other Variables
	/// Boolean variable to check the existence of the supremizer modes
	bool supex;

    // Dummy variables to transform simplefoam into a class
	/// Pressure field
	autoPtr<dgScalarField> _rho;

    /// Velocity field
    autoPtr<dgVectorField> _rhoU;

    /// Energy Field
    autoPtr<dgScalarField> _Ener;

    autoPtr<dgScalarField> _p;

    autoPtr<dgVectorField> _U;

	/// Mesh
	mutable autoPtr<dgMesh> _mesh;
// 	autoPtr<dgMesh> _mesh;

    /*
	/// simpleControl
	autoPtr<simpleControl> _simple;
*/
	/// Time
	autoPtr<Time> _runTime;	
	

	///
	/// solve the supremizer either with the use of the pressure snaphots or the pressure modes
	///
	/// @param[in]  type  The type of the supremizer approach, either done on the pressure snapshots or on pressure modes.
	///
	void solvesupremizer(word type="snapshots");
    
	/// Perform a lift solve
	void liftSolve();

    //
    // Functions
    //
    
    //--------------------------------------------------------------------------        
    /// @brief      Perform a truthsolve
    ///
    /// @param[in]  mu_now  The actual value of the parameter, used only to write it 
    /// in the Offline folder together with the snapshot.
    ///
    /// @param      mu  value of the parameter
    ///
    void truthSolve(List<scalar> mu);


    //--------------------------------------------------------------------------        
    /// Function to check if the solution must be exported.
    ///
    /// @param      timeObject  The time object of OpenFOAM.
    ///
    /// @return     1 if the time target has been reached, 0 elsewhere.
    ///
    bool checkWrite(Time& timeObject);
    
    //--------------------------------------------------------------------------        
    ///     Write and create the initialConditions file using the given input data.
    ///
    /// @param      rho_bc: velocity component along the x-axis 
    ///
    /// @param      U_bc: velocity component along the x-axis 
    ///
    /// @param      V_bc: velocity component along the y-axis
    ///    
    void writeBCfiles(double rho_bc, double U_bc,  double V_bc );

};

#endif
