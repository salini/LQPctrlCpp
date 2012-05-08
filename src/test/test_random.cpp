

#include "Solver.h"

#include <time.h>
#include <stdlib.h>
#include <iostream>


void printSolverSolution(GenericSolver& s)
{
    std::cout<<"cost: "<<s.getCost()<<std::endl;
    std::cout<<"x   : ";
    s.printVector(s.getSolution(), s.getnProblem());
    std::cout<<std::endl;
}


void eye(double* mat, unsigned int n)
{
    for (int i=0; i<n*n; i++)
    {
        mat[i] = 0;
    }
    for (int i=0; i<n; i++)
        mat[i*n+i] = 1;
}

void randVec(double* vec, unsigned int nDim)
{
    for (int i=0; i<nDim; i++)
        vec[i] = (rand() - (RAND_MAX/2.))/(rand()+1);
}

void randMat(double* mat, unsigned int nRow, unsigned int nCol)
{
    for (int i=0; i<nRow*nCol; i++)
        mat[i] = (rand() - (RAND_MAX/2.))/(rand()+1);
}

void createProblem(unsigned int nP, unsigned int nE, unsigned int nI, double* _P, double* _q, double* _A, double* _b, double* _G, double* _h)
{
    eye(_P, nP);
    randVec(_q, nP);
    randMat(_A, nE, nP);
    randVec(_b, nE);
    randMat(_G, nI, nP);
    randVec(_h, nI);
}


int main(int argc, char** argv)
{
    srand ( time(NULL) );
    
    unsigned int nProblem=6;
    unsigned int nEquality=3;
    unsigned int nInequality=2;
    
    if (argc == 4)
    {
        nProblem = atoi(argv[1]);
        nEquality = atoi(argv[2]);
        nInequality = atoi(argv[3]);
    }
    
    double* P = new double[nProblem*nProblem];
    double* q = new double[nProblem];
    double* A = new double[nEquality*nProblem];
    double* b = new double[nEquality];
    double* G = new double[nInequality*nProblem];
    double* h = new double[nInequality];
    
    createProblem(nProblem, nEquality, nInequality, P, q, A, b, G, h);
    
    GenericSolver solver;

    solver.initProblem(nProblem, nEquality, nInequality);
    
//    solver.printMatrix(P, nProblem, nProblem);
    
    solver.setProblem(P, q, A, b, G, h);
    solver.printProblem();
    
    

    std::cout<<"==================================================="<<std::endl;
    std::cout<<"USE QuadProg++:"<<std::endl;
    solver.setSolver("QuadProg");
    solver.solveProblem();
    printSolverSolution(solver);
    
    
//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE CGAL:"<<std::endl;
//    solver.setSolver("CGAL");
//    solver.solveProblem();
//    printSolverSolution(solver);
    
    
//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE OOQP:"<<std::endl;
//    solver.setSolver("OOQP");
//    solver.solveProblem();
//    printSolverSolution(solver);
//    
//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE COPL_QP:"<<std::endl;
//    solver.setSolver("COPL_QP");
//    solver.solveProblem();
//    printSolverSolution(solver);

//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE Mosek:"<<std::endl;
//    solver.setSolver("Mosek");
//    solver.solveProblem();
//    printSolverSolution(solver);


#if defined WIN32
    System("pause");
#endif

    return 0;
}


