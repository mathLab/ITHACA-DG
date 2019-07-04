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


/// \file
/// Source file of the reductionProblem class.


#include "reductionProblem.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reductionProblem::reductionProblem()
{
	offline = ITHACAutilitiesDG::check_off();
// 	podex = ITHACAutilitiesDG::check_pod();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set the parameters
void reductionProblem::setParameters()
{
	mu.resize(Pnumber, Tnumber);
	mu_range.resize(Pnumber, 2);
}

/*
// Generate Random Parameters
void reductionProblem::genRandPar()
{
	for (int k = 0; k < Pnumber; k++)
	{
		for (int n = 0; n < Tnumber; ++n)
		{
			mu(k, n) = mu_range(k, 0) + static_cast <float> (rand()) / static_cast <float> (RAND_MAX / (mu_range(k, 1) - mu_range(k, 0))); // generate numbers
		}
	}
}

void reductionProblem::genRandPar(int Tsize)
{
	mu.resize(Tsize, Pnumber);
	std::srand(std::time(0));
	auto rand = Eigen::MatrixXd::Random(Tsize, Pnumber);
	auto dx = mu_range.col(1) - mu_range.col(0);
	Eigen::MatrixXd k;
	k = (rand.array() + double(1)) / 2;

	for (int i = 0; i < Pnumber; i ++)
	{
		k.col(i) = (k.col(i) * dx(i)).array() + mu_range(i, 0);

	}
	mu = k;
}
*/

// Generate Equidistributed Parameters
void reductionProblem::genEquiPar()
{
	Eigen::VectorXd vec;

	for (int k = 0; k < Pnumber; k++)
	{
		mu.row(k) = vec.LinSpaced(Tnumber, mu_range(k, 0), mu_range(k, 1));
	}
}


// Perform a TruthSolve (To be overridden)
void reductionProblem::truthSolve()
{
    Info << "reductionProblem::truthSolve() is a Method to be overridden -> Exiting the code" << endl;
    exit(0);
}



// Change type of BC
void reductionProblem::changeBCtype(dgVectorField& field, word BCtype, label BC_ind)
{
    field.boundaryFieldRef().set(BC_ind, Foam::dgPatchField<vector>::New(BCtype, field.mesh().boundary()[BC_ind], field));
}

void reductionProblem::changeBCtype(dgScalarField& field, word BCtype, label BC_ind)
{
    field.boundaryFieldRef().set(BC_ind, Foam::dgPatchField<scalar>::New(BCtype, field.mesh().boundary()[BC_ind], field));
}




// Assign a BC for a vector field
void reductionProblem::assignBC(dgScalarField& s, label BC_ind, double& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "zeroGradient")
    {
        zeroGradientDgPatchScalarField& Tpatch = refCast<zeroGradientDgPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
        // NB not all the BC are supported in HopeFOAM
//     else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
//     {
//         fixedGradientFvPatchScalarField& Tpatch = refCast<fixedGradientFvPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
//         scalarField& gradTpatch = Tpatch.gradient();
//         forAll(gradTpatch, faceI)
//         {
//             gradTpatch[faceI] = value;
//         }
//     }
//     else if (s.boundaryField()[BC_ind].type() == "fixedFluxPressure")
//     {
//         for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
//         {
//             s.boundaryFieldRef()[BC_ind][i] = value;
//         }
//     }
//     else if (s.boundaryField()[BC_ind].type() == "freestream")
//     {
//         for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
//         {
//             s.boundaryFieldRef()[BC_ind][i] = value;
//         }
//         freestreamFvPatchField<scalar>& Tpatch = refCast<freestreamFvPatchField<scalar> >(s.boundaryFieldRef()[BC_ind]);
//         scalarField& gradTpatch = Tpatch.freestreamValue();
//         forAll(gradTpatch, faceI)
//         {
//             gradTpatch[faceI] = value;
//         }
// 
//     }
    else if (s.boundaryField()[BC_ind].type() == "calculated")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "empty")
    {

    }
}

// Assign a BC for a scalar field
void reductionProblem::assignBC(dgVectorField& s, label BC_ind, Vector<double>& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "zeroGradient")
    {
        Info << "This Feature is not implemented for this boundary condition" << endl;
        exit(0);
    }
//     else if (s.boundaryField()[BC_ind].type() == "freestream")
//     {
//         for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
//         {
//             s.boundaryFieldRef()[BC_ind][i] = value;
//         }
//         freestreamFvPatchField<vector>& Tpatch = refCast<freestreamFvPatchField<vector> >(s.boundaryFieldRef()[BC_ind]);
//         vectorField& gradTpatch = Tpatch.freestreamValue();
//         forAll(gradTpatch, faceI)
//         {
//             gradTpatch[faceI] = value;
//         }
// 
//     }
}

void reductionProblem::writeMu(List<scalar> mu_now)
{
	mkDir("./ITHACAoutput/Parameters");
	std::ofstream ofs;
	ofs.open ("./ITHACAoutput/Parameters/par", std::ofstream::out | std::ofstream::app);
	forAll(mu_now, i)
	{
		ofs << mu_now[i] << ' ';
	}
	ofs << "\n";
	ofs.close();
}


// ************************************************************************* //

