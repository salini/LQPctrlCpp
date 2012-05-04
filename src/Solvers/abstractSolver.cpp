
#include "Solvers/abstractSolver.h"


AbstractSolver::AbstractSolver()
{
}

AbstractSolver::~AbstractSolver()
{

}

double AbstractSolver::getCost()
{
    return cost;
}

Eigen::VectorXd AbstractSolver::getSolution()
{
    return x;
}

void AbstractSolver::initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality)
{
    nProblem    = _nProblem;
    nEquality   = _nEquality;
    nInequality = _nInequality;
    
    x.resize(nProblem);
}

void AbstractSolver::setCostFunction(Eigen::MatrixXd& _P, Eigen::VectorXd& _q, double _c)
{

}

void AbstractSolver::setEqualityConstraint(Eigen::MatrixXd& _A, Eigen::VectorXd& _b)
{

}

void AbstractSolver::setInequalityConstraint(Eigen::MatrixXd& _G, Eigen::VectorXd& _h)
{

}

bool AbstractSolver::solveProblem()
{
    return false;
}
