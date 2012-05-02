#ifndef __LQPCTRLCPP_ABSTRACTSOLVER__
#define __LQPCTRLCPP_ABSTRACTSOLVER__


class AbstractSolver
{
public:
    AbstractSolver();
    virtual ~AbstractSolver();
    
    virtual void initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality);
    virtual void freeSolver();
    virtual void setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h);
    virtual void solveProblem();

protected:
    unsigned int nProblem;
    unsigned int nEquality;
    unsigned int nInequality;
};


#include "QuadProg++.hh"

class solver_QuadProg: public AbstractSolver
{
public:
    solver_QuadProg();
    virtual ~solver_QuadProg();
    
    void initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality);
    void freeSolver();
    void setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h);
    void solveProblem();

private:
    Matrix<double> P;   // cost function
    Vector<double> q;
    
    Matrix<double> A;   // equality constraints
    Vector<double> b;
    
    Matrix<double> G;   // inequality constraints
    Vector<double> h;
    
    Vector<double> x;   // solution vector
};




#include <CGAL/basic.h>
#include <CGAL/QP_models.h>

class solver_CGAL: public AbstractSolver
{
public:
    solver_CGAL();
    virtual ~solver_CGAL();
    
    void initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality);
    void freeSolver();
    void setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h);
    void solveProblem();

private:
    unsigned int nConstraint;

    double** P;     // cost function
    double*  q;
    
    double** AG;    // it concatenates the equality and inequality constraints
    double*  bh;
    
    CGAL::Comparison_result* r;
    
    bool*  fl;
    bool*  fu;
    double*  l;
    double*  u;
    
    double* x;       // solution vector
};



class solver_COPL_QP: public AbstractSolver
{
public:
    solver_COPL_QP();
    virtual ~solver_COPL_QP();
    
    void initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality);
    void freeSolver();
    void setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h);
    void solveProblem();
};


#endif
