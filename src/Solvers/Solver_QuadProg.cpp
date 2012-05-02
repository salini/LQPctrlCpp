
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
//    nProblem    = _nProblem;
//    nEquality   = _nEquality;
//    nInequality = _nInequality;

    P.resize(nProblem, nProblem);
    q.resize(nProblem);
    
    A.resize(nProblem, nEquality);
    b.resize(nEquality);
    
    G.resize(nProblem, nInequality);
    h.resize(nInequality);
    
    x.resize(nProblem);
}

void solver_QuadProg::freeSolver()
{

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
void solver_QuadProg::setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h)
{
    int i, j;
    
    double _Pcw[nProblem*nProblem];
    double _Acw[nProblem*nProblem];
    double _Gcw[nProblem*nProblem];
    
    for (i=0; i<nProblem; i++)
        for (j=0; j<nProblem; j++)
            _Pcw[j*nProblem+i] = _P[i*nProblem+j];
    
    for (i=0; i<nEquality; i++)
        for (j=0; j<nProblem; j++)
            _Acw[j*nEquality+i] = -_A[i*nProblem+j];

    for (i=0; i<nInequality; i++)
        for (j=0; j<nProblem; j++)
            _Gcw[j*nInequality+i] = -_G[i*nProblem+j];


    P.set(_Pcw, nProblem, nProblem);
    q.set(_q, nProblem);
    
    A.set(_Acw, nProblem, nEquality);
    b.set(_b, nEquality);
    
    G.set(_Gcw, nProblem, nInequality);
    h.set(_h, nInequality);
    
//    std::cout<<"P:"<<std::endl<<P<<std::endl;
//    std::cout<<"q:"<<std::endl<<q<<std::endl;
//    std::cout<<"A:"<<std::endl<<A<<std::endl;
//    std::cout<<"b:"<<std::endl<<b<<std::endl;
//    std::cout<<"G:"<<std::endl<<G<<std::endl;
//    std::cout<<"h:"<<std::endl<<h<<std::endl;
}


void solver_QuadProg::solveProblem()
{
    std::cout << "f: " << solve_quadprog(P, q, A, b, G, h, x) << std::endl;
    std::cout << "x: " << x << std::endl;
}
#else

void solver_QuadProg::initSolver(int _nProblem, int _nEquality, int _nInequality)
{

}

void solver_QuadProg::freeSolver()
{

}

void solver_QuadProg::setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h)
{

}

void solver_QuadProg::solve()
{
    std::cout<<"QuadProg++ solver has not been implemented..."<<std::endl;
}
#endif

