
#include "Solver.h"

#include <iostream>
#include <fstream>


////////////////////////////////////////////////////////////////////////////////


void GenericSolver::initProblem(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality)
{
    nProblem    = _nProblem;
    nEquality   = _nEquality;
    nInequality = _nInequality;

    P.setZero(nProblem, nProblem);       // cost function
    q.setZero(nProblem);                 //
    c = 0;
    
    A.setZero(nEquality, nProblem);      // equality constraint
    b.setZero(nEquality);
    
    G.setZero(nInequality, nProblem);    // inequality constraint
    h.setZero(nInequality);
    
    currentSolver->initSolver(_nProblem, _nEquality, _nInequality);
    problemSolved = false;
}


GenericSolver::GenericSolver()
{
    solverName = "";
    currentSolver=NULL;
    setSolver(new solver_QuadProg);
    initProblem(0,0,0);
}


GenericSolver::~GenericSolver()
{
    if (currentSolver!=NULL)
    {
        delete currentSolver;
    }
}


void GenericSolver::setSolver(std::string _solverName)
{
    if (_solverName == "QuadProg")
    {
        setSolver(new solver_QuadProg);
    }
    else
    {
        std::cout<<_solverName<<" solver is unknown. Keep previous solver."<<std::endl;
    }
}

void GenericSolver::setSolver(AbstractSolver* _newSolver)
{
    if (currentSolver!=NULL)
    {
        delete currentSolver;
    }
    currentSolver = _newSolver;
}

void GenericSolver::setCostFunction(double* _P, double* _q, double _c)
{
    Eigen::Map<Eigen::MatrixXd> _Eigen_P(_P, nProblem, nProblem);
    Eigen::Map<Eigen::VectorXd> _Eigen_q(_q, nProblem);
    
    P = _Eigen_P;
    q = _Eigen_q;
    setCostFunction(P, q, _c);
}

void GenericSolver::setEqualityConstraint(double* _A, double* _b)
{
    Eigen::Map<Eigen::MatrixXd> _Eigen_A(_A, nEquality, nProblem);
    Eigen::Map<Eigen::VectorXd> _Eigen_b(_b, nEquality);
    
    A = _Eigen_A;
    b = _Eigen_b;
    setEqualityConstraint(A, b);
}

void GenericSolver::setInequalityConstraint(double* _G, double* _h)
{
    Eigen::Map<Eigen::MatrixXd> _Eigen_G(_G, nInequality, nProblem);
    Eigen::Map<Eigen::VectorXd> _Eigen_h(_h, nInequality);
    
    G = _Eigen_G;
    h = _Eigen_h;
    setInequalityConstraint(G, h);
}

void GenericSolver::setCostFunction(Eigen::MatrixXd& _P, Eigen::VectorXd& _q, double _c)
{
    P = _P;
    q = _q;
    c = _c;
    
    currentSolver->setCostFunction(P, q, c);
    problemSolved = false;
}

void GenericSolver::setEqualityConstraint(Eigen::MatrixXd& _A, Eigen::VectorXd& _b)
{
    A = _A;
    b = _b;
    
    currentSolver->setEqualityConstraint(A, b);
    problemSolved = false;
}

void GenericSolver::setInequalityConstraint(Eigen::MatrixXd& _G, Eigen::VectorXd& _h)
{
    G = _G;
    h = _h;
    
    currentSolver->setInequalityConstraint(G,h);
    problemSolved = false;
}


void GenericSolver::setProblem(Eigen::MatrixXd& _P, Eigen::VectorXd& _q, Eigen::MatrixXd& _A, Eigen::VectorXd& _b, Eigen::MatrixXd& _G, Eigen::VectorXd& _h, double _c)
{
    setCostFunction(_P, _q, _c);
    setEqualityConstraint(_A, _b);
    setInequalityConstraint(_G, _h);
}


bool GenericSolver::solveProblem()
{
    if (nProblem == 0)
    {
        std::cout<<"nProblem is 0. Nothing to optimize."<<std::endl;
        problemSolved = false;
    }
    else
    {
        problemSolved = currentSolver->solveProblem();
    }
    return problemSolved;
}



void GenericSolver::setProblemFromFile(std::string _fileName)
{
    int i,j;
    std::string output;

    std::ifstream fileIn(_fileName.c_str());

    if(fileIn.is_open())
    {
        while (!fileIn.eof())
        {
            fileIn >> output;
            if (output == "nProblem:")
                fileIn >> nProblem;
            else if (output == "nEquality:")
                fileIn >> nEquality;
            else if (output == "nInequality:")
                fileIn >> nInequality;
        }
    }
    fileIn.close();

    initProblem(nProblem, nEquality, nInequality);

    std::cout<<"IN2"<<std::endl;
    fileIn.open(_fileName.c_str(), std::ifstream::in);
    if(fileIn.is_open())
    {
        while (!fileIn.eof())
        {
            fileIn >> output;
            std::cout<<output<<std::endl;
            if (output == "P:")
            {
                for (i=0; i<nProblem; i++)
                {
                    for (j=0; j<nProblem; j++)
                    {
                        fileIn >> P(i,j);
                    }
                }
            }
            else if (output == "q:")
            {
                for (i=0; i<nProblem; i++)
                {
                    fileIn >> q(i);
                }
            }
            else if (output == "c:")
            {
                fileIn >> c;
            }
            else if (output == "A:")
            {
                for (i=0; i<nEquality; i++)
                {
                    for (j=0; j<nProblem; j++)
                    {
                        fileIn >> A(i,j);
                    }
                }
            }
            else if (output == "b:")
            {
                for (i=0; i<nEquality; i++)
                {
                    fileIn >> b(i);
                }
            }
            else if (output == "G:")
            {
                for (i=0; i<nInequality; i++)
                {
                    for (j=0; j<nProblem; j++)
                    {
                        fileIn >> G(i,j);
                    }
                }
            }
            else if (output == "h:")
            {
                for (i=0; i<nInequality; i++)
                {
                    fileIn >> h(i);
                }
            }
        }
    }
    setProblem(P, q, A, b, G, h, c);
    problemSolved = false;
}

void GenericSolver::printProblem()
{
    std::cout << "nProblem:    " << nProblem   << std::endl
              << "nEquality:   " << nEquality  << std::endl
              << "nInequality: " << nInequality<< std::endl;
    
    std::cout<<"P:\n"<<P<<std::endl;
    std::cout<<"q:\n"<<q<<std::endl;
    std::cout<<"c:\n"<<c<<std::endl;
    
    std::cout<<"A:\n"<<A<<std::endl;
    std::cout<<"b:\n"<<b<<std::endl;
    
    std::cout<<"G:\n"<<G<<std::endl;
    std::cout<<"h:\n"<<h<<std::endl;
}


bool GenericSolver::problemIsSolved()
{
    return problemSolved;
}

double GenericSolver::getCost()
{
    if (!problemIsSolved())
    {
        std::cout<<"WARNING: Problem unsolved. Cost not computed."<<std::endl;
    }
    return currentSolver->getCost();
}

Eigen::VectorXd GenericSolver::getSolution()
{
    if (!problemIsSolved())
    {
        std::cout<<"WARNING: Problem unsolved. Optimal value not computed."<<std::endl;
    }
    return currentSolver->getSolution();
}

unsigned int GenericSolver::getnProblem()
{
    return nProblem;
}

unsigned int GenericSolver::getnEquality()
{
    return nEquality;
}

unsigned int GenericSolver::getnInequality()
{
    return nInequality;
}





