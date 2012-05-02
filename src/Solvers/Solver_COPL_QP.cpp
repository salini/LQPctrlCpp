
#include "Solvers/abstractSolver.h"



solver_COPL_QP::solver_COPL_QP()
{

}

solver_COPL_QP::~solver_COPL_QP()
{

}

//-----------------------------------------------------------------------------------------//
#if defined COPL_QP_IS_AVAILABLE



void solver_COPL_QP::initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality)
{
    AbstractSolver::initSolver(_nProblem, _nEquality, _nInequality);

}

void solver_COPL_QP::freeSolver()
{

}
void solver_COPL_QP::setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h)
{

}


void solver_COPL_QP::solveProblem()
{

}

#else

void solver_COPL_QP::initSolver(int _nProblem, int _nEquality, int _nInequality)
{

}

void solver_COPL_QP::freeSolver()
{

}

void solver_COPL_QP::setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h)
{

}

void solver_COPL_QP::solve()
{
    std::cout<<"QuadProg++ solver has not been implemented..."<<std::endl;
}
#endif

