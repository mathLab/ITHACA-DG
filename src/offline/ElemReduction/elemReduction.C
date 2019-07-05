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
/// source file for the elemReduction class

#include "elemReduction.H"
// #include "EigenFunctions.H"


void elemReduction::reduce(PtrList<dgVectorField>& snapshotsU, bool sup)
{
//         ITHACAparameters para;

        // Cycle on the snapshots
        int nDof;
        int dofStart;
        int nElems = snapshotsU[0].mesh().cellElementsTree().size();
          
        //
        // From DG to Vol
        // 
        for (int iSnap =0;  iSnap< snapshotsU.size();  iSnap++)
        {
            //
            dgVectorField U_tmp(snapshotsU[0].name(), snapshotsU[iSnap]);
            U_tmp.resize(nElems);
            // Reset the inner field
            U_tmp *= scalar(0);
            
            // Cycle on the cells of the domain 
            for(dgTree<physicalCellElement>::leafIterator iter = snapshotsU[iSnap].mesh().cellElementsTree().leafBegin(); iter != snapshotsU[iSnap].mesh().cellElementsTree().end(); ++iter)
            {
                dofStart = iter()->value().dofStart();
                nDof = iter()->value().nDof();
                label nQuadPoint = iter()->value().cellVandermonde().size()/nDof;
                label nPhyPoints = nDof;
                
                Eigen::MatrixXd cellV = Eigen::MatrixXd::Zero(nQuadPoint, nPhyPoints );
                Eigen::VectorXd CellValues = Eigen::VectorXd::Zero(nPhyPoints );
                Eigen::VectorXd CellQuadValues = Eigen::VectorXd::Zero(nQuadPoint );
                // Convert the cellVandermonde matrix to Eigen format
                for(int nRow=0; nRow < nQuadPoint; nRow++)
                {
                    for(int nCol=0; nCol < nPhyPoints; nCol++)
                    {
                        cellV(nRow, nCol) =  iter()->value().cellVandermonde()[nRow*nDof+nCol];
                    }
                }
                // Cycle on the Variable components
                for (int nVar=0; nVar < 3; nVar++)
                {
                    for(int nCol=0; nCol < nPhyPoints; nCol++)
                    {
                        CellValues(nCol) =  snapshotsU[iSnap][dofStart+nCol][nVar];
                    }
                    CellQuadValues = cellV*CellValues;
                    Eigen::VectorXd JwValues = Eigen::VectorXd::Zero(nQuadPoint );
                    for(int nCol=0; nCol < nQuadPoint; nCol++)
                    {
                        JwValues(nCol)  = iter()->value().jacobianWeights()[nCol];
                    }
                    Eigen::VectorXd averageUvalue = Eigen::VectorXd::Zero(1);
                    Eigen::MatrixXd JwValuesT = JwValues.transpose();
                    averageUvalue = JwValuesT*CellQuadValues;
                    averageUvalue = averageUvalue/JwValuesT.sum();
                    U_tmp[iter()->value().sequenceIndex()][nVar] = averageUvalue(0) ;
                }
            }

            // Output reduced field
            elemReduction::exportReducedSolution(U_tmp, name(iSnap+1), sup);                        
    } 
}


void elemReduction::reduce(PtrList<dgScalarField>& snapshotsP, bool sup)
{

//         ITHACAparameters para;

        // Initialization of the integers and of the hyperreduced matrix
        int nDof;
        int dofStart;
        int nElems = snapshotsP[0].mesh().cellElementsTree().size();

        //
        //Evaluation of the elemental modes
        //
        // Cycle on the snapshots
        for (int iSnap =0;  iSnap< snapshotsP.size();  iSnap++)
        {
            dgScalarField p_tmp(snapshotsP[0].name(), snapshotsP[iSnap]);

            p_tmp.resize(nElems);
            // Reset the inner field
            p_tmp *= scalar(0);

            // Cycle on the cells of the domain 
            for(dgTree<physicalCellElement>::leafIterator iter = snapshotsP[iSnap].mesh().cellElementsTree().leafBegin(); iter != snapshotsP[iSnap].mesh().cellElementsTree().end(); ++iter)
            {
                dofStart = iter()->value().dofStart();
                nDof = iter()->value().nDof();
                label nQuadPoint = iter()->value().cellVandermonde().size()/nDof;
                label nPhyPoints = nDof;

                Eigen::MatrixXd cellV = Eigen::MatrixXd::Zero(nQuadPoint, nPhyPoints );
                Eigen::VectorXd CellValues = Eigen::VectorXd::Zero(nPhyPoints );
                Eigen::VectorXd CellQuadValues = Eigen::VectorXd::Zero(nQuadPoint );
                // Convert the cellVandermonde matrix to Eigen format
                for(int nRow=0; nRow < nQuadPoint; nRow++)
                {
                    for(int nCol=0; nCol < nPhyPoints; nCol++)
                    {
                        cellV(nRow, nCol) =  iter()->value().cellVandermonde()[nRow*nDof+nCol];
                    }
                }
                // Convert local cell values to an Eigen Vector
                for(int nCol=0; nCol < nPhyPoints; nCol++)
                {
                    CellValues(nCol) =  snapshotsP[iSnap][dofStart+nCol];
                }
                CellQuadValues = cellV*CellValues;
                Eigen::VectorXd JwValues = Eigen::VectorXd::Zero(nQuadPoint );
                for(int nCol=0; nCol < nQuadPoint; nCol++)
                {
                     JwValues(nCol)  = iter()->value().jacobianWeights()[nCol];
                }
                Eigen::VectorXd averagePvalue = Eigen::VectorXd::Zero(1);
                Eigen::MatrixXd JwValuesT = JwValues.transpose();
                averagePvalue = JwValuesT*CellQuadValues;
                averagePvalue = averagePvalue/JwValuesT.sum();
                p_tmp[iter()->value().sequenceIndex()] = averagePvalue(0) ;
            }
            // Output reduced field
            elemReduction::exportReducedSolution(p_tmp, name(iSnap+1), sup);
        }
    
}



/// Export the Bases
void elemReduction::exportBases(PtrList<dgVectorField>& s, PtrList<dgVectorField>& _snapshots, bool sup)
{
    if (sup)
    {
        fileName fieldname;
        for (label i = 0; i < s.size(); i++)
        {
            mkDir("./ITHACAoutput/supremizer/" + name(i + 1));
            fieldname = "./ITHACAoutput/supremizer/" + name(i + 1) + "/" + _snapshots[i].name();
            word type = "volVectorField";
            // Open the file
            OFstream os(fieldname);
            // Write the Header
            Foam::IOobject:: writeBanner(os)
                << "FoamFile\n{\n"
                << "    version     " << os.version() << ";\n"
                << "    format      " << os.format() << ";\n"
                << "    class       " << type << ";\n";
            os  << "    location    " << name(i + 1) << ";\n"
                <<"    object      " << _snapshots[i].name() << ";\n"
                << "}" << nl;
            os  <<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" 
            <<nl;
            // Write the data
            s[i].writeData(os);
        }
    }
    else
    {
        fileName fieldname;
        for (label i = 0; i < s.size(); i++)
        {
            mkDir("./ITHACAoutput/POD/" + name(i + 1));
            fieldname = "./ITHACAoutput/POD/" + name(i + 1) + "/" + _snapshots[i].name();
            word type = "volVectorField";
            // Open the file
            OFstream os(fieldname);
            // Write the Header
            Foam::IOobject:: writeBanner(os)
                << "FoamFile\n{\n"
                << "    version     " << os.version() << ";\n"
                << "    format      " << os.format() << ";\n"
                << "    class       " << type << ";\n";
            os  << "    location    " << name(i + 1) << ";\n"
                <<"    object      " << _snapshots[i].name() << ";\n"
                << "}" << nl;
            os  <<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" 
            <<nl;
            // Write the data
            s[i].writeData(os);
        }
    }
}


/// Export the Bases
void elemReduction::exportBases(PtrList<dgScalarField>& s, PtrList<dgScalarField>& _snapshots, bool sup)
{
    if (sup)
    {
        fileName fieldname;
        for (label i = 0; i < s.size(); i++)
        {
            mkDir("./ITHACAoutput/supremizer/" + name(i + 1));
            fieldname = "./ITHACAoutput/supremizer/" + name(i + 1) + "/" + _snapshots[i].name();
            word type = "volScalarField";
            // Open the file
            OFstream os(fieldname);
            // Write the Header
            Foam::IOobject:: writeBanner(os)
                << "FoamFile\n{\n"
                << "    version     " << os.version() << ";\n"
                << "    format      " << os.format() << ";\n"
                << "    class       " << type << ";\n";
            os  << "    location    " << name(i + 1) << ";\n"
                <<"    object      " << _snapshots[i].name() << ";\n"
                << "}" << nl;
            os  <<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" 
            <<nl;
            // Write the data
            s[i].writeData(os);
        }
    }
    else
    {
        fileName fieldname;
        for (label i = 0; i < s.size(); i++)
        {
            mkDir("./ITHACAoutput/POD/" + name(i + 1));
            fieldname = "./ITHACAoutput/POD/" + name(i + 1) + "/" + _snapshots[i].name();
            word type = "volScalarField";
            // Open the file
            OFstream os(fieldname);
            // Write the Header
            Foam::IOobject:: writeBanner(os)
                << "FoamFile\n{\n"
                << "    version     " << os.version() << ";\n"
                << "    format      " << os.format() << ";\n"
                << "    class       " << type << ";\n";
            os  << "    location    " << name(i + 1) << ";\n"
                <<"    object      " << _snapshots[i].name() << ";\n"
                << "}" << nl;
            os  <<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" 
            <<nl;
            // Write the data
//             _snapshots[i].writeHeader(os, type);
            s[i].writeData(os);
        }
    }
}


/// Export the reduced solutions
void elemReduction::exportReducedSolution(dgVectorField& s, fileName name, bool sup)
{
    fileName fieldname;
    int c =0;

    if (sup)
    {  
        struct stat sb;
        if (stat("./ITHACAoutput/supfield", &sb) == 0 && S_ISDIR(sb.st_mode))
        {
        }
        else
        {
            mkDir("./ITHACAoutput/supfield");
            c += abs(system("ln -s ../../constant ./ITHACAoutput/supfield/constant"));
            c += abs(system("ln -s ../../0 ./ITHACAoutput/supfield/0"));
            c += abs(system("ln -s ../../system ./ITHACAoutput/supfield/system"));
        } 
        mkDir("./ITHACAoutput/supfield/" + name);
        fieldname = "./ITHACAoutput/supfield/" + name + "/" + s.name();
    }
    else
    {
        struct stat sb;

        if (stat("./ITHACAoutput/Offline", &sb) == 0 && S_ISDIR(sb.st_mode))
        {
//             Info << "./ITHACAoutput/Offline already exist!" << endl;
        }
        else
        {
            mkDir("./ITHACAoutput/Offline");
            c += abs(system("ln -s ../../constant ./ITHACAoutput/Offline/constant"));
            c += abs(system("ln -s ../../0 ./ITHACAoutput/Offline/0"));
            c += abs(system("ln -s ../../system ./ITHACAoutput/Offline/system"));
        }
        mkDir("./ITHACAoutput/Offline/" + name);
        fieldname = "./ITHACAoutput/Offline/" + name + "/" + s.name();
    }
    word type = "volVectorField";
    // Open the file
    OFstream os(fieldname);
    //     s.writeHeader(os, type);   
    // Write the Header
    Foam::IOobject:: writeBanner(os)
        << "FoamFile\n{\n"
        << "    version     " << os.version() << ";\n"
        << "    format      " << os.format() << ";\n"
        << "    class       " << type << ";\n";
    os  << "    location    " << name<< ";\n"
        <<"    object      " << s.name() << ";\n"
        << "}" << nl;
    os  <<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" 
    <<nl;
    // Write the data
    s.writeData(os);
}


/// Export the reduced solutions
void elemReduction::exportReducedSolution(dgScalarField& s, fileName name, bool sup)
{
    #include "IOobject.H"
    #include "objectRegistry.H"

    fileName fieldname;
    struct stat sb;
    int c =0;
    if (stat("./ITHACAoutput/Offline", &sb) == 0 && S_ISDIR(sb.st_mode))
    {
//         Info << "./ITHACAoutput/Offline already exist!" << endl;
    }
    else
    {
        mkDir("./ITHACAoutput/Offline");
        c += abs(system("ln -s ../../constant ./ITHACAoutput/Offline/constant"));
        c += abs(system("ln -s ../../0 ./ITHACAoutput/Offline/0"));
        c += abs(system("ln -s ../../system ./ITHACAoutput/Offline/system"));
    }
    mkDir("./ITHACAoutput/Offline/" + name);
    fieldname = "./ITHACAoutput/Offline/" + name + "/" + s.name();
    word type = "volScalarField";
    // Open the file
    OFstream os(fieldname);
//     s.writeHeader(os, type);   
    // Write the Header
    Foam::IOobject:: writeBanner(os)
        << "FoamFile\n{\n"
        << "    version     " << os.version() << ";\n"
        << "    format      " << os.format() << ";\n"
        << "    class       " << type << ";\n";
    os  << "    location    " << name<< ";\n"
        <<"    object      " << s.name() << ";\n"
        << "}" << nl;
    os  <<"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n" 
    <<nl;
    // Write the data
    s.writeData(os);
}



// * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * * //







