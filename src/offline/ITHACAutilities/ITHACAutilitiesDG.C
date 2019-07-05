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

#include "ITHACAutilitiesDG.H"

/// \file
/// Source file of the ITHACAutilitiesDG class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Check if the modes exists
bool ITHACAutilitiesDG::check_pod()
{
    struct stat sb;
    bool pod_exist;
    int c = 0;
    if (stat("./ITHACAoutput/POD", &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        pod_exist = true;
        Info << "POD data already exist, reading existing modes" << endl;
    }
    else
    {
        pod_exist = false;
        Info << "POD don't exist, performing a POD decomposition" << endl;
        mkDir("./ITHACAoutput/POD");
        c += abs(system("ln -s ../../constant ./ITHACAoutput/POD/constant"));
        c += abs(system("ln -s ../../0 ./ITHACAoutput/POD/0"));
        c += abs(system("ln -s ../../system ./ITHACAoutput/POD/system"));
    }
    if (c > 0)
    {
        Info << "Error in system operation" << endl;
        exit(0);
    }
    return pod_exist;
}

// Check if the offline data exist
bool ITHACAutilitiesDG::check_off()
{
    struct stat sb;
    bool off_exist;
    int c = 0;
    if (stat("./ITHACAoutput/OfflineDG", &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        off_exist = true;
        Info << "Offline data already exist, reading existing data" << endl;
    }
    else
    {
        off_exist = false;
        Info << "Offline don't exist, performing the Offline Solve" << endl;
//         mkDir("./ITHACAoutput/Offline");
//         c += abs(system("ln -s ../../constant ./ITHACAoutput/Offline/constant"));
//         c += abs(system("ln -s ../../0 ./ITHACAoutput/Offline/0"));
//         c += abs(system("ln -s ../../system ./ITHACAoutput/Offline/system"));
        
        mkDir("./ITHACAoutput/OfflineDG");
        c += abs(system("ln -s ../../constant ./ITHACAoutput/OfflineDG/constant"));
        c += abs(system("ln -s ../../0 ./ITHACAoutput/OfflineDG/0"));
        c += abs(system("ln -s ../../system ./ITHACAoutput/OfflineDG/system"));
    }
    if (c > 0)
    {
        Info << "Error in system operation" << endl;
        exit(0);
    }
    else
    {
    }
return off_exist;
}


// Check if the supremizer data exist
bool ITHACAutilitiesDG::check_sup()
{
    struct stat sb;
    bool sup_exist;
    int c = 0;
    if (stat("./ITHACAoutput/supremizer", &sb) == 0 && S_ISDIR(sb.st_mode))
    {
        sup_exist = true;
        Info << "Supremizer data already exist, reading existing data" << endl;
    }
    else
    {
        sup_exist = false;
        Info << "Supremizers don't exist, performing a POD decomposition" << endl;
        mkDir("./ITHACAoutput/supremizer");
        c += abs(system("ln -s ../../constant ./ITHACAoutput/supremizer/constant"));
        c += abs(system("ln -s ../../0 ./ITHACAoutput/supremizer/0"));
        c += abs(system("ln -s ../../system ./ITHACAoutput/supremizer/system"));
    }
    if (c > 0)
    {
        Info << "Error in system operation" << endl;
        exit(0);
    }
    return sup_exist;
}


// Assign a BC for a vector field
void ITHACAutilitiesDG::assignBC(dgScalarField& s, label BC_ind, double& value)
{
    if (s.boundaryField()[BC_ind].type() == "fixedValue")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
//     else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
//     {
//         fixedGradientDgPatchScalarField& Tpatch = refCast<fixedGradientDgPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
//         scalarField& gradTpatch = Tpatch.gradient();
//         forAll(gradTpatch, faceI)
//         {
//             gradTpatch[faceI] = value;
//         }
//     }
    else if (s.boundaryField()[BC_ind].type() == "zeroGradient")
    {
        zeroGradientDgPatchScalarField& Tpatch = refCast<zeroGradientDgPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            gradTpatch[faceI] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "fixedFluxPressure")
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
//     else if (s.boundaryField()[BC_ind].type() == "freestream")
//     {
//         for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
//         {
//             s.boundaryFieldRef()[BC_ind][i] = value;
//         }
//         freestreamDgPatchField<scalar>& Tpatch = refCast<freestreamDgPatchField<scalar> >(s.boundaryFieldRef()[BC_ind]);
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


void ITHACAutilitiesDG::assignBC(dgScalarField& s, label BC_ind, Eigen::MatrixXd valueVec)
{
    word typeBC = s.boundaryField()[BC_ind].type();
    if (typeBC == "fixedValue" || typeBC == "calculated" || typeBC == "fixedFluxPressure" )
    {
        for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
        {
            double value = valueVec(i);
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "zeroGradient")
    {
        zeroGradientDgPatchScalarField& Tpatch = refCast<zeroGradientDgPatchScalarField>(s.boundaryFieldRef()[BC_ind]);
        scalarField& gradTpatch = Tpatch.gradient();
        forAll(gradTpatch, faceI)
        {
            double value = valueVec(faceI);
            gradTpatch[faceI] = value;
        }
    }

    else if (s.boundaryField()[BC_ind].type() == "empty")
    {

    }
}

// Assign a BC for a vector field
void ITHACAutilitiesDG::assignBC(dgVectorField& s, label BC_ind, Eigen::MatrixXd valueVec)
{
    word typeBC = s.boundaryField()[BC_ind].type();
    int sizeBC = s.boundaryField()[BC_ind].size();
    if (typeBC == "fixedValue" || typeBC == "calculated")
    {
        for (label i = 0; i < sizeBC; i++)
        {
            Vector<double> value(valueVec(i), valueVec(i + sizeBC), valueVec(i + sizeBC * 2));
            s.boundaryFieldRef()[BC_ind][i] = value;
        }
    }
    else if (s.boundaryField()[BC_ind].type() == "zeroGradient")
    {
        Info << "This Feature is not implemented for this boundary condition" << endl;
        exit(0);
    }
// NB. Not all the OpenFOAM  boundary conditons are implemented in HopeFOAM
//     else if (s.boundaryField()[BC_ind].type() == "fixedGradient")
//     {
//         Info << "This Feature is not implemented for this boundary condition" << endl;
//         exit(0);
//     }
//     else if (s.boundaryField()[BC_ind].type() == "freestream")
//     {
//         for (label i = 0; i < s.boundaryField()[BC_ind].size(); i++)
//         {
//             Vector<double> value(valueVec(i), valueVec(i + sizeBC), valueVec(i + sizeBC * 2));
//             s.boundaryFieldRef()[BC_ind][i] = value;
//         }
//         freestreamDgPatchField<vector>& Tpatch = refCast<freestreamDgPatchField<vector> >(s.boundaryFieldRef()[BC_ind]);
//         vectorField& gradTpatch = Tpatch.freestreamValue();
//         forAll(gradTpatch, faceI)
//         {
//             Vector<double> value(valueVec(faceI), valueVec(faceI + sizeBC), valueVec(faceI + sizeBC * 2));
//             gradTpatch[faceI] = value;
//         }
//     }
}


// * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //

