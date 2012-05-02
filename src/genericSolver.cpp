
#include "genericSolver.h"

#include <iostream>
#include <fstream>


////////////////////////////////////////////////////////////////////////////////
void GenericSolver::initProblem()
{
    nProblem    = 0;
    nEquality   = 0;
    nInequality = 0;

    c = 0;        // constant in cost function
    P = NULL; // cost function
    q = NULL;
    
    A = NULL; // equality constraint
    b = NULL;
    
    G = NULL; // inequality constraint
    h = NULL;
    
    solverName = "";
    currentSolver = new AbstractSolver();
}


void GenericSolver::freeProblem()
{
    if (P!=NULL) delete[] P;
    if (q!=NULL) delete[] q;
    
    if (A!=NULL) delete[] A;
    if (b!=NULL) delete[] b;
    
    if (G!=NULL) delete[] G;
    if (h!=NULL) delete[] h;
}


GenericSolver::GenericSolver()
{
    currentSolver=NULL;
    initProblem();
}


GenericSolver::~GenericSolver()
{
    freeProblem();
}


void GenericSolver::setSolver(std::string _solverName)
{
    if (currentSolver!=NULL)
    {
        delete currentSolver;
    }
    
    if (_solverName == "QuadProg")
    {
        currentSolver = new solver_QuadProg;
    }
    else if (_solverName == "CGAL")
    {
        currentSolver = new solver_CGAL;
    }
    else
    {
        std::cout<<_solverName<<" solver is unknown."<<std::endl;
        currentSolver = new AbstractSolver;
    }
    
    currentSolver->initSolver(nProblem, nEquality, nInequality);
    currentSolver->setProblem(P, q, A, b, G, h);
}

void GenericSolver::setSolver(AbstractSolver* _newSolver)
{
    currentSolver = _newSolver;
}





void GenericSolver::solveProblem()
{
    currentSolver->solveProblem();
}

void GenericSolver::setProblem()
{

}

void GenericSolver::setProblemFromFile(std::string _fileName)
{
    freeProblem();
    initProblem();

    std::string output;
    
    std::ifstream fileIn(_fileName.c_str());
    
    if(fileIn.is_open())
    {
        while (!fileIn.eof()) {
            fileIn >> output;
            if (output == "nProblem:")
                fileIn >> nProblem;
            else if (output == "nEquality:")
                fileIn >> nEquality;
            else if (output == "nInequality:")
                fileIn >> nInequality;

            else if (output == "P:")
            {
                P = new double[nProblem*nProblem];
                for (int i=0; i<nProblem*nProblem; i++)
                {
                    fileIn >> P[i];
                }
            }
            else if (output == "q:")
            {
                q = new double[nProblem];
                for (int i=0; i<nProblem; i++)
                {
                    fileIn >> q[i];
                }
            }
            else if (output == "c:")
            {
                fileIn >> c;
            }
            
            else if (output == "A:")
            {
                A = new double[nEquality*nProblem];
                for (int i=0; i<nEquality*nProblem; i++)
                {
                    fileIn >> A[i];
                }
            }
            else if (output == "b:")
            {
                b = new double[nEquality];
                for (int i=0; i<nEquality; i++)
                {
                    fileIn >> b[i];
                }
            }
            else if (output == "G:")
            {
                G = new double[nInequality*nProblem];
                for (int i=0; i<nInequality*nProblem; i++)
                {
                    fileIn >> G[i];
                }
            }
            else if (output == "h:")
            {
                h = new double[nInequality];
                for (int i=0; i<nInequality; i++)
                {
                    fileIn >> h[i];
                }
            }
        }
    }
}

void GenericSolver::printProblem()
{
    std::cout << "nProblem:    " << nProblem   << std::endl
              << "nEquality:   " << nEquality  << std::endl
              << "nInequality: " << nInequality<< std::endl;
    
    std::cout<<"P:"<<std::endl;
    printMatrix(nProblem, nProblem, P);
    std::cout<<std::endl;
    std::cout<<"q:"<<std::endl;
    printVector(nProblem, q);
    std::cout<<std::endl;
    std::cout<<"c:"<<std::endl<<c<<std::endl;
    
    std::cout<<std::endl;
    std::cout<<std::endl;
    
    std::cout<<"A:"<<std::endl;
    printMatrix(nEquality, nProblem, A);
    std::cout<<std::endl;
    std::cout<<"b:"<<std::endl;
    printVector(nEquality, b);
    
    std::cout<<std::endl;
    std::cout<<std::endl;
    
    std::cout<<"G:"<<std::endl;
    printMatrix(nInequality, nProblem, G);
    std::cout<<std::endl;
    std::cout<<"h:"<<std::endl;
    printVector(nInequality, h);
}


void GenericSolver::printMatrix(unsigned int nRow, unsigned int nCol, double* mat)
{
    for (int i=0; i<nRow; i++)
    {
        for (int j=0; j<nCol; j++)
        {
            std::cout<<mat[i*nCol+j]<<"  ";
        }
        std::cout<<std::endl;
    }
}

void GenericSolver::printMatrix(unsigned int nRow, unsigned int nCol, double** mat)
{
    for (int i=0; i<nRow; i++)
    {
        for (int j=0; j<nCol; j++)
        {
            std::cout<<mat[i][j]<<"  ";
        }
        std::cout<<std::endl;
    }
}

void GenericSolver::printVector(unsigned int nDim, double* vec)
{
    for (int i=0; i<nDim; i++)
    {
        std::cout<<vec[i]<<"  "<<std::endl;
    }
}
