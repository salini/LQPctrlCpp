
#include "Solvers/abstractSolver.h"



solver_CGAL::solver_CGAL()
{

}

solver_CGAL::~solver_CGAL()
{

}

//-----------------------------------------------------------------------------------------//
#if defined CGAL_IS_AVAILABLE


#include <cassert>
#include <CGAL/QP_functions.h>
//#include <CGAL/QP_solution.h>

// choose exact integral type
//#if defined USE_GMP_WIH_CGAL
//    #if defined CGAL_USE_GMP
//        #include <CGAL/Gmpz.h>
//        typedef CGAL::Gmpz ET;
//    #else
//        #include <CGAL/MP_Float.h>
//        typedef CGAL::MP_Float ET;
//    #endif
//#else
    #include <CGAL/MP_Float.h>
    typedef CGAL::MP_Float ET;
//#endif


typedef CGAL::Quadratic_program_from_iterators
<double**,                                                // for A
 double*,                                                 // for b
 CGAL::Comparison_result*,                                // for r
 bool*,                                                   // for fl
 double*,                                                 // for l
 bool*,                                                   // for fu
 double*,                                                 // for u
 double**,                                                // for D
 double*>                                                 // for c 
Program;


typedef CGAL::Quadratic_program_solution<ET> Solution;











void solver_CGAL::initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality)
{
    
    AbstractSolver::initSolver(_nProblem, _nEquality, _nInequality);
    nConstraint = nEquality + nInequality;

    unsigned int i;
    
    P = new double*[nProblem];
    for (i=0; i<nProblem; i++)
        P[i] = new double[i+1];
    q = new double[nProblem];

    AG = new double*[nProblem];
    for (i=0; i<nProblem; i++)
        AG[i] = new double[nConstraint];
    bh = new double[nConstraint];
    
    r = new CGAL::Comparison_result[nConstraint];
    for (i=0; i<nEquality; i++)
        r[i] = CGAL::EQUAL;
    for (i=nEquality; i<nConstraint; i++)
        r[i] = CGAL::SMALLER;

    fl = new bool[nProblem];
    fu = new bool[nProblem];
    l  = new double[nProblem];
    u  = new double[nProblem];
    for (i=0; i<nProblem; i++)
    {
        fl[i] = false;
        fu[i] = false;
        l[i]  = 0;
        u[i]  = 0;
    }
    
    x = new double[nProblem];
}

void solver_CGAL::freeSolver()
{

}


void solver_CGAL::setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h)
{
    unsigned int i,j;
    
    for (i=0; i<nProblem; i++)
        for (j=0; j<i+1; j++)
            P[i][j] = _P[i*nProblem+j]; // DO NOT forget the 2 constant!!!

    for (i=0; i<nProblem; i++)
        q[i] = _q[i];             //TODO: maybe optimization by replacing the pointer, not copying all the value?!?

    for (i=0; i<nProblem; i++)
        for (j=0; j<nEquality; j++)
            AG[i][j] = _A[j*nProblem+i];
    for (i=0; i<nProblem; i++)
        for (j=0; j<nInequality; j++)
            AG[i][nEquality+j] = _G[j*nProblem+i];

    for (i=0; i<nEquality; i++)
        bh[i] = _b[i];
    for (i=0; i<nInequality; i++)
        bh[nEquality+i] = _h[i];

////////////////////////////////////////////////////////////////////////////////
    std::cout<<"P:"<<std::endl;
    for (i=0; i<nProblem; i++)
    {
        for (j=0; j<i+1; j++)
        {
            std::cout<<P[i][j]<<" ";
        }
        std::cout<<std::endl;
    }

    std::cout<<"q:"<<std::endl;
    for (i=0; i<nProblem; i++)
        std::cout<<q[i]<<" ";
    std::cout<<std::endl;

    std::cout<<"AG:"<<std::endl;
    for (i=0; i<nProblem; i++)
    {
        for (j=0; j<nConstraint; j++)
        {
            std::cout<<AG[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    
    std::cout<<"bh:"<<std::endl;
    for (i=0; i<nConstraint; i++)
        std::cout<<bh[i]<<" ";
    std::cout<<std::endl;
    
}




void solver_CGAL::solveProblem()
{
    std::cout<<"solve Problem"<<std::endl;
    // now construct the quadratic program; the first two parameters are
    // the number of variables and the number of constraints (rows of A)
    Program qp (nProblem, nConstraint, AG, bh, r, fl, l, fu, u, P, q, 0); // the last term is c0, the constant term

    // solve the program, using ET as the exact type
    Solution s = CGAL::solve_quadratic_program(qp, ET());

    // output solution
    //std::cout << s;
    std::cout<<to_double(s.objective_value())<<std::endl;
    
    Solution::Variable_value_iterator it  = s.variable_values_begin();
    Solution::Variable_value_iterator end = s.variable_values_end();
    for (; it != end; ++it)
        std::cout << to_double(*it) << " ";
    std::cout << std::endl;
}
#else

void solver_CGAL::initSolver(unsigned int _nProblem, unsigned int _nEquality, unsigned int _nInequality)
{

}

void solver_CGAL::freeSolver()
{

}

void solver_CGAL::setProblem(double* _P, double* _q, double* _A, double* _b, double* _G, double* _h)
{

}

void solver_CGAL::solve()
{
    std::cout<<"CGAL solver has not been implemented..."<<std::endl;
}
#endif


