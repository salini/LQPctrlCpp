
#include "Solvers/abstractSolver.h"


AbstractSolver::AbstractSolver()
{

}

AbstractSolver::~AbstractSolver()
{

}

void AbstractSolver::initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality)
{
    nProblem    = _nProblem;
    nEquality   = _nEquality;
    nInequality = _nInequality;
}

void AbstractSolver::freeSolver()
{

}

void AbstractSolver::setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h)
{

}

void AbstractSolver::solveProblem()
{

}
