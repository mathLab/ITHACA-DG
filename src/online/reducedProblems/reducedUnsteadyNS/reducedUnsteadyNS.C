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
/// Source file of the reducedUnsteadyNS class


#include "reducedUnsteadyNS.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor initialization
reducedUnsteadyNS::reducedUnsteadyNS()
{
}

reducedUnsteadyNS::reducedUnsteadyNS(unsteadyNS& FOMproblem)
{
    problem = &FOMproblem;

    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();

    // Create locally the velocity modes
    for (label k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append(problem->liftfield[k]);
    }
    for (label k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append(problem->Umodes[k]);
    }
    for (label k = 0; k < problem->NSUPmodes; k++)
    {
        Umodes.append(problem->supmodes[k]);
    }

    // Create locally the pressure modes
    for (label k = 0; k < problem->NPmodes; k++)
    {
        Pmodes.append(problem->Pmodes[k]);
    }

    // Store locally the snapshots for projections
    for (label k = 0; k < problem->Ufield.size(); k++)
    {
        Usnapshots.append(problem->Ufield[k]);
        Psnapshots.append(problem->Pfield[k]);
    }
    newton_object_sup = newton_unsteadyNS_sup(Nphi_u + Nphi_p , Nphi_u + Nphi_p, FOMproblem);
    newton_object_PPE = newton_unsteadyNS_PPE(Nphi_u + Nphi_p , Nphi_u + Nphi_p, FOMproblem);
}

// * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //

int newton_unsteadyNS_sup::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
     Eigen::VectorXd a_tmp_old(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    a_tmp_old = y_old.head(Nphi_u);
//         a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    a_dot = (1.5*x.head(Nphi_u) - 2.0* y_old.head(Nphi_u) + 0.5*y_old_old.head(Nphi_u)) / dt;
    Eigen::MatrixXd gg(1, 1);

    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->P_matrix * a_tmp;
    Eigen::VectorXd M6 = problem->D_matrix * b_tmp;

    for (label i = 0; i < Nphi_u; i++)
    {
//         cc = a_tmp.transpose() * problem->C_matrix[i] * a_tmp;
        cc = 2.0* a_tmp.transpose() * problem->C_matrix[i] * a_tmp- a_tmp_old.transpose() * problem->C_matrix[i] * a_tmp_old;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M2(i);
//         fvec(i) = - M5(i) + M1(i) - cc(0, 0) ; // - M2(i);
    }
    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        fvec(k) =  M3(j) ;
    }
    for (label j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - BC(j);
    }
    return 0;
}



// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyNS_sup::df(const Eigen::VectorXd &x,  Eigen::MatrixXd &fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyNS_sup> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}



int newton_unsteadyNS_PPE::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd a_tmp_old(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    a_tmp_old = y_old.head(Nphi_u);

    // Convective terms
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->D_matrix * b_tmp;

    // BC PPE
    Eigen::VectorXd M6 = problem->BC1_matrix * a_tmp * nu;

    // BC PPE
    Eigen::VectorXd M7 = problem->BC3_matrix * a_tmp * nu;

    for (label i = 0; i < Nphi_u; i++)
    {
//         cc = a_tmp.transpose() * problem->C_matrix[i] * a_tmp;
        cc = 2.0* a_tmp.transpose() * problem->C_matrix[i] * a_tmp- a_tmp_old.transpose() * problem->C_matrix[i] * a_tmp_old;
        fvec(i) = - M5(i) + M1(i) - cc(0, 0) - M2(i);
    }
    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        gg = a_tmp.transpose() * problem->G_matrix[j] * a_tmp;
        bb = a_tmp.transpose() * problem->BC2_matrix[j] * a_tmp;
        fvec(k) = M3(j, 0) + gg(0, 0) - M7(j, 0);
    }
    for (label j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - BC(j);
    }
    return 0;
}


// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyNS_PPE::df(const Eigen::VectorXd &x,  Eigen::MatrixXd &fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyNS_PPE> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}


// * * * * * * * * * * * * * Solve Functions supremizer * * * * * * * * * * * //

void reducedUnsteadyNS::solveOnline_sup(Eigen::MatrixXd& vel_now, label startSnap)
{
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    // y.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[startSnap], Umodes);
    // y.tail(Nphi_p) = ITHACAutilities::get_coeffs(problem->Pfield[startSnap], Pmodes);
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(Usnapshots[startSnap], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs(Psnapshots[startSnap], Pmodes);

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = vel_now(j, 0);
    }

    // Set some properties of the newton object
    newton_object_sup.nu = nu;
    newton_object_sup.y_old_old = y;
    newton_object_sup.y_old = y;

//     newton_object_sup.y_oldOld = y;
    newton_object_sup.dt = dt;
    newton_object_sup.BC.resize(N_BC);

    for (label j = 0; j < N_BC; j++)
    {
        newton_object_sup.BC(j) = vel_now(j, 0);
    }


    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    online_solution.resize(Ntsteps);

    // Set the initial time
    time = tstart;

    // Counting variable
    int counter = 0;

    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    if(time!=0)
    {
        online_solution[counter] = tmp_sol;
        counter ++;
    }

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNS_sup> hnls(newton_object_sup);

    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    
    time = time + dt;
    while (time < finalTime)
    {
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);
        for (label j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
        newton_object_sup.operator()(y, res);
        newton_object_sup.y_old_old = newton_object_sup.y_old;
        newton_object_sup.y_old = y;


        std::cout << "################## Online solve N° " << count_online_solve << " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << vel_now << std::endl;
        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        count_online_solve += 1;
        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;
        if (counter >= online_solution.size())
        {
            online_solution.append(tmp_sol);
        }
        else
        {
            online_solution[counter] = tmp_sol;
        }
        counter ++;
        time = time + dt;
    }
    // Export the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");
    count_online_solve += 1;
}



void reducedUnsteadyNS::solveOnline_ch(Eigen::MatrixXd& vel_now, label startSnap)
{
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    y.head(Nphi_u) = ITHACAutilities::get_coeffs(Usnapshots[startSnap], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs(Psnapshots[startSnap], Pmodes);
    
     // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = vel_now(j, 0);
    }

     // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    online_solution.resize(Ntsteps);

     // Set the initial time
    time = tstart;
// 
//     // Counting variable
    int counter = 0;

//     // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    
    if(time!=0)
    {
        online_solution[counter] = tmp_sol;
        counter ++;
    }

    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    
    // Intialization
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd a_tmp_old(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    Eigen::VectorXd y_old;
    Eigen::VectorXd y_old_old ;

    Eigen::MatrixXd gg(1, 1);

    // Convective term
    Eigen::MatrixXd cc(1, 1);
    
    Eigen::VectorXd UT(Nphi_u);
    Eigen::VectorXd MUT(Nphi_u);
    Eigen::VectorXd UTT(Nphi_u);
    Eigen::VectorXd MUTT(Nphi_u);
    Eigen::VectorXd PT(Nphi_p);
    Eigen::VectorXd MPT(Nphi_p);
    Eigen::VectorXd U(Nphi_u);

     // Parameters
    scalar gamma0 = 0;
    scalar a0 =  1.0;
    scalar a1 =  0.0;
    scalar b0 =  1.0;
    scalar b1 =  0.0;
    
    // Intial values for the y_old vector
    y_old = y;

    time = time + dt;
    while (time <= finalTime)
    {
        a_tmp = y.head(Nphi_u) ;
        b_tmp = y.tail(Nphi_p) ;
        
        // Update the val,ues
        y_old_old = y_old;
        y_old = y;
 
        a_tmp_old = y_old.head(Nphi_u);
        
        a_dot = (- a0* y_old.head(Nphi_u) -a1*y_old_old.head(Nphi_u)) / dt;
            
        // Mom Term
        Eigen::VectorXd M1 = problem->B_matrix * a_tmp * nu;
        // Gradient of pressure
        Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
        // Mass Term
        Eigen::VectorXd M5 = problem->M_matrix * a_dot;
        
        // Pressure Term
        Eigen::VectorXd M3 = problem->P_matrix * a_tmp;
        Eigen::VectorXd M6 = problem->D_matrix * b_tmp;

        //Laplacian velocity term
//         Eigen::MatrixXd Diag_Contr= 1.5/dt/nu*Eigen::MatrixXd::Identity(Nphi_u, Nphi_u);
        Eigen::MatrixXd Diag_Contr= 1.5/dt/nu*problem->M_matrix;
        Eigen::MatrixXd LD_matrix = -problem->L_matrix+Diag_Contr;
//         Eigen::VectorXd M8 = LD_matrix * a_tmp;
        
        if(gamma0==0.0)
        {
            gamma0 = 1.5;
            a0 =  2.0;
            a1 = -0.5;
            b0 =  2.0;
            b1 = -1.0;
        }

        // Step 1: velocity guess
        for (label i = 0; i < Nphi_u; i++)
        {
            cc = b0* a_tmp.transpose() * problem->C_matrix[i] * a_tmp + b1* a_tmp_old.transpose() * problem->C_matrix[i] * a_tmp_old;
//             UT(i) = ( ( a0*a_tmp_old(i) + a1*a_tmp_old(i) ) - dt*cc(0, 0) ) /gamma0;
                MUT(i) = (-M5(i) - cc(0, 0) ) /gamma0*dt;
        }

        UT = problem->M_matrix.inverse() *MUT;
        // Step 2: evaluation of the pressure
        for (label i = 0; i < Nphi_p; i++)
        {
            gg = UT.transpose() * problem->G_matrix[i] * UT;
//         gg = UT.transpose() * problem->C_matrix[i] * UT;
            
            Eigen::VectorXd M7 = problem->BC3_matrix * a_tmp * nu* gamma0;
    
//             MPT(i) = - gamma0/dt*gg(0,0);
                MPT(i) =  -gg(0,0) + M7(i, 0);

        }
//         
        PT = problem->D_matrix.inverse()*MPT;

        // Step 3: velocity correction
        for (label i = 0; i < Nphi_u; i++)
        {
            UT(i) = UT(i)*gamma0/dt;
        }
        MUTT = problem->M_matrix*(UT) - (problem->K_matrix*PT);
//          MUTT = problem->M_matrix*(UT) - (problem->K_matrix*b_tmp);

    
        for (label i = 0; i < Nphi_u; i++)
        {
            MUTT(i) = MUTT(i)*dt/gamma0;
        }
        UTT = problem->M_matrix.inverse() *MUTT;

        // Step 4: computation of the final velocity  
        UTT = UTT*gamma0/dt/nu;
        U = LD_matrix.inverse() * UTT;

        y.head(Nphi_u) = U;
        y.tail(Nphi_p) = PT;

        for (label j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }

        std::cout << "################## Online solve N° " << count_online_solve << " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << vel_now << std::endl;

        count_online_solve += 1;
        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;
        if (counter >= online_solution.size())
        {
            online_solution.append(tmp_sol);
        }
        else
        {
            online_solution[counter] = tmp_sol;
        }
        counter ++;
        time = time + dt;
    }

    count_online_solve += 1;
}




// * * * * * * * * * * * * * * * Solve Functions PPE * * * * * * * * * * * * * //

void reducedUnsteadyNS::solveOnline_PPE(Eigen::MatrixXd& vel_now, label startSnap)
{
    // Create and resize the solution vector
    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    // Set Initial Conditions
    // y.head(Nphi_u) = ITHACAutilities::get_coeffs(problem->Ufield[startSnap], Umodes);
    // y.tail(Nphi_p) = ITHACAutilities::get_coeffs(problem->Pfield[startSnap], Pmodes);
    y.head(Nphi_u) = ITHACAutilities::get_coeffs(Usnapshots[startSnap], Umodes);
    y.tail(Nphi_p) = ITHACAutilities::get_coeffs(Psnapshots[startSnap], Pmodes);

    // Change initial condition for the lifting function
    for (label j = 0; j < N_BC; j++)
    {
        y(j) = vel_now(j, 0);
    }

    // Set some properties of the newton object
    newton_object_PPE.nu = nu;
    newton_object_PPE.y_old = y;
    newton_object_PPE.y_old_old = y;
    newton_object_PPE.dt = dt;
    newton_object_PPE.BC.resize(N_BC);

    for (label j = 0; j < N_BC; j++)
    {
        newton_object_PPE.BC(j) = vel_now(j, 0);
    }

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    online_solution.resize(Ntsteps);

    // Set the initial time
    time = tstart;

    // Counting variable
    int counter = 0;

    // Create vectpr to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_p + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;
    if(time!=0)
    {
        online_solution[counter] = tmp_sol;
        counter ++;
    }

    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_unsteadyNS_PPE> hnls(newton_object_PPE);

    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    time = time + dt;

    // Start the time loop
    while (time < finalTime + 0.99*dt)
    {
        Info << time <<  "     " << dt << "    "<<finalTime <<endl;

        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);
        for (label j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
        newton_object_PPE.operator()(y, res);
        newton_object_PPE.y_old_old = newton_object_PPE.y_old;
        newton_object_PPE.y_old = y;

        std::cout << "################## Online solve N° " << count_online_solve << " ##################" << std::endl;
        Info << "Time = " << time << endl;
        std::cout << "Solving for the parameter: " << vel_now << std::endl;
        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        count_online_solve += 1;
        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;
        if (counter >= online_solution.size())
        {
            online_solution.append(tmp_sol);
        }
        else
        {
            online_solution[counter] = tmp_sol;
        }
        counter ++;
        time = time + dt;
    }

    // Export the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python", "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab", "./ITHACAoutput/red_coeff");
    count_online_solve += 1;
}

void reducedUnsteadyNS::reconstruct_PPE(fileName folder, int printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);

    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);
            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }
            problem->exportSolution(U_rec,  name(counter2), folder);

            volScalarField P_rec("P_rec", problem->Pmodes[0] * 0);
            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += problem->Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }
            problem->exportSolution(P_rec, name(counter2), folder);
            nextwrite += printevery;
            
            double timenow = online_solution[i](0,0);
            std::ofstream of(folder + name(counter2) + "/" + name(timenow));    

            counter2 ++;

            UREC.append(U_rec);
            PREC.append(P_rec);
        }
        counter++;
    }
}

void reducedUnsteadyNS::reconstruct_sup(fileName folder, int printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);

    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);
            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }
            problem->exportSolution(U_rec,  name(counter2), folder);

            volScalarField P_rec("P_rec", problem->Pmodes[0] * 0);
            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += problem->Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }
            problem->exportSolution(P_rec, name(counter2), folder);
            nextwrite += printevery;

            double timenow = online_solution[i](0,0);
            std::ofstream of(folder + name(counter2) + "/" + name(timenow));   
            
            counter2 ++;

            UREC.append(U_rec);
            PREC.append(P_rec);
        }
        counter++;
    }
}

void reducedUnsteadyNS::reconstruct_ch(fileName folder, int printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);

    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);
            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }
            problem->exportSolution(U_rec,  name(counter2), folder);

            volScalarField P_rec("P_rec", problem->Pmodes[0] * 0);
            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += problem->Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }
            problem->exportSolution(P_rec, name(counter2), folder);
            nextwrite += printevery;

            double timenow = online_solution[i](0,0);
            std::ofstream of(folder + name(counter2) + "/" + name(timenow));   
            
            counter2 ++;

            UREC.append(U_rec);
            PREC.append(P_rec);
        }
        counter++;
    }
}

// ************************************************************************* //
