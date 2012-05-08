#ifndef __LQPCTRLCPP_GENERICSOLVER__
#define __LQPCTRLCPP_GENERICSOLVER__


#include <string>

#include <Eigen/Dense>

#include "Solvers/abstractSolver.h"



class GenericSolver
{
public:
    GenericSolver();
    virtual ~GenericSolver();
    
    void setSolver(std::string _solverName);
    void setSolver(AbstractSolver* _newSolver);
    
    void initProblem(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality);
    
    void setCostFunction(double* _P, double* _q, double _c=0);
    void setEqualityConstraint(double* _A, double* _b);
    void setInequalityConstraint(double* _G, double* _h);
    void setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h, double _c=0);

    void setCostFunction(Eigen::MatrixXd& _P, Eigen::VectorXd& _q, double _c=0);
    void setEqualityConstraint(Eigen::MatrixXd& _A, Eigen::VectorXd& _b);
    void setInequalityConstraint(Eigen::MatrixXd& _G, Eigen::VectorXd& _h);
    void setProblem(Eigen::MatrixXd& _P, Eigen::VectorXd& _q, Eigen::MatrixXd& _A, Eigen::VectorXd& _b, Eigen::MatrixXd& _G, Eigen::VectorXd& _h, double _c=0);
    
    bool solveProblem();
    
    void setProblemFromFile(std::string _fileName);
    void printProblem();

    unsigned int getnProblem();
    unsigned int getnEquality();
    unsigned int getnInequality();

    double  getCost();
    Eigen::VectorXd getSolution();
    
    bool    problemIsSolved();

private:

    bool problemSolved;

    std::string solverName;
    
    unsigned int nProblem, nEquality, nInequality;

    Eigen::MatrixXd P; // cost function
    Eigen::VectorXd q; //
    double          c; //
    
    Eigen::MatrixXd A; // Equality constraints
    Eigen::VectorXd b; //
    
    Eigen::MatrixXd G; // Inequality constraints
    Eigen::VectorXd h; //
    
    AbstractSolver* currentSolver;
};



#endif
