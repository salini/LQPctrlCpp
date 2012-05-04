#ifndef __LQPCTRLCPP_ABSTRACTSOLVER__
#define __LQPCTRLCPP_ABSTRACTSOLVER__

#include <Eigen/Dense>

class AbstractSolver
{
public:
    AbstractSolver();
    virtual ~AbstractSolver();

    virtual void initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality);
    virtual bool solveProblem();
    virtual void setCostFunction(Eigen::MatrixXd& _P, Eigen::VectorXd& _q, double _c=0);
    virtual void setEqualityConstraint(Eigen::MatrixXd& _A, Eigen::VectorXd& _b);
    virtual void setInequalityConstraint(Eigen::MatrixXd& _G, Eigen::VectorXd& _h);

    double          getCost();
    Eigen::VectorXd getSolution();

protected:
    unsigned int nProblem;
    unsigned int nEquality;
    unsigned int nInequality;
    
    double  cost;
    Eigen::VectorXd x;
};




#include "QuadProg++.hh"

class solver_QuadProg: public AbstractSolver
{
public:
    solver_QuadProg();
    virtual ~solver_QuadProg();
    
    void initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality);
    bool solveProblem();
    void setCostFunction(Eigen::MatrixXd& _P, Eigen::VectorXd& _q, double _c=0);
    void setEqualityConstraint(Eigen::MatrixXd& _A, Eigen::VectorXd& _b);
    void setInequalityConstraint(Eigen::MatrixXd& _G, Eigen::VectorXd& _h);
    

private:
    QuadProgPP::Matrix<double> P;   // cost function
    QuadProgPP::Vector<double> q;
    double                     c;
    
    QuadProgPP::Matrix<double> A;   // equality constraints
    QuadProgPP::Vector<double> b;
    
    QuadProgPP::Matrix<double> G;   // inequality constraints
    QuadProgPP::Vector<double> h;
    
    QuadProgPP::Vector<double> QuadProgx;   // solution vector
};


class solver_Mosek: public AbstractSolver
{
public:
    solver_Mosek();
    virtual ~solver_Mosek();
    
    void initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality);
    bool solveProblem();
    void setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h, double _c=0);
    
};

#endif
