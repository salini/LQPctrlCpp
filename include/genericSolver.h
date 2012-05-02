#ifndef __LQPCTRLCPP_GENERICSOLVER__
#define __LQPCTRLCPP_GENERICSOLVER__


#include <string>

#include "Solvers/abstractSolver.h"

//enum LQPSOLVER {
//                QUADPROG,
//                CGAL,
//                MOSEK,
//                COPL_QP,
//                OOQP
//               };


class GenericSolver
{
public:
    GenericSolver();
    virtual ~GenericSolver();
    
    void setSolver(std::string _solverName);
    void setSolver(AbstractSolver* _newSolver);
    
    void initProblem();
    void freeProblem();
    void setProblem();
    void solveProblem();
    
    void setProblemFromFile(std::string _fileName);
    void printProblem();

private:

//    void initSolver();
//    void freeSolver();

    void printMatrix(unsigned int nRow, unsigned int nCol, double* mat); //TODO: row-wise or column-wise
    void printMatrix(unsigned int nRow, unsigned int nCol, double** mat);
    void printVector(unsigned int nDim, double* vec);
    
    std::string       solverName;
    
    unsigned int nProblem, nEquality, nInequality;

    
    double* P; // cost function
    double* q;
    double  c; // constant value of cost function

    double* A; // equality constraint
    double* b;
    
    double* G; // inequality constraint
    double* h;
    
    
    AbstractSolver* currentSolver;
};



#endif
