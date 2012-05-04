
#include "Solvers/abstractSolver.h"



solver_QuadProg::solver_QuadProg()
{

}

solver_QuadProg::~solver_QuadProg()
{

}

//-----------------------------------------------------------------------------------------//
#if defined QUADPROG_IS_AVAILABLE

#include <iostream>
#include <sstream>
#include <string>


void solver_QuadProg::initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality)
{
    AbstractSolver::initSolver(_nProblem, _nEquality, _nInequality);

    P.resize(nProblem, nProblem);
    q.resize(nProblem);
    c = 0;
    
    A.resize(nProblem, nEquality);
    b.resize(nEquality);
    
    G.resize(nProblem, nInequality);
    h.resize(nInequality);
    
    QuadProgx.resize(nProblem);
}


/*! set the QuadProg Problem
 *
 * QuadProg solve the following problem:
 * 
 * (1/2) * X.P.X + q.X + c
 *
 * st: A.X + b  = 0
 *     G.X + h >= 0
 *
 */
void solver_QuadProg::setCostFunction(Eigen::MatrixXd& _P, Eigen::VectorXd& _q, double _c)
{
    P.set(_P.data(), nProblem, nProblem);
    q.set(_q.data(), nProblem);
    c = _c;
}

void solver_QuadProg::setEqualityConstraint(Eigen::MatrixXd& _A, Eigen::VectorXd& _b)
{
    A.set(_A.data(), nProblem, nEquality);
    b.set(Eigen::VectorXd(-_b).data(), nEquality);
}


void solver_QuadProg::setInequalityConstraint(Eigen::MatrixXd& _G, Eigen::VectorXd& _h)
{
    G.set(Eigen::MatrixXd(-_G).data(), nProblem, nInequality);
    h.set(_h.data(), nInequality);
}

#include <stdexcept>

bool solver_QuadProg::solveProblem()
{
    try
    {
        cost = QuadProgPP::solve_quadprog(P, q, A, b, G, h, QuadProgx);
    }
    catch ( std::logic_error logErr )
    {
        std::cout<<"ERROR: QuadProg solver failed! "<<logErr.what()<<std::endl;
        return false;
    }
    catch (std::runtime_error logErr )
    {
        std::cout<<"ERROR: QuadProg solver failed! "<<logErr.what()<<std::endl;
        return false;
    }
    
    for (int i=0; i<nProblem; i++)
        x[i] = QuadProgx[i];

    return true;
    
}
#else

void solver_QuadProg::initSolver(int _nProblem, int _nEquality, int _nInequality)
{

}

void solver_QuadProg::freeSolver()
{

}

void solver_QuadProg::setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h, double _c)
{

}

bool solver_QuadProg::solveProblem()
{
    std::cout<<"QuadProg++ solver has not been implemented..."<<std::endl;
    return false;
}
#endif

