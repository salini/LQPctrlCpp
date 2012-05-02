

#include "genericSolver.h"

//#include "Solvers/abstractSolver.h"

#include <iostream>

int main(int argc, char** argv)
{
    int exitVal=0;
    
    GenericSolver solver;
    
//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE QUADPROG:"<<std::endl;
//    solver.setSolver("QuadProg");
//    std::cout<<"solve"<<std::endl;
//    solver.solveProblem();

//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE CGAL:"<<std::endl;
//    solver.setSolver("CGAL");
//    solver.solveProblem();




//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE MOSEK:"<<std::endl;
//    solver.setSolver("MOSEK");
//    solver.solveProblem();

//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE COPL_QP:"<<std::endl;
//    solver.setSolver("COPL_QP");
//    solver.solveProblem();

//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE OOQP:"<<std::endl;
//    solver.setSolver("OOQP");
//    solver.solveProblem();

////////////////////////////////////////////////////////////////////////////////

//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE CGAL:"<<std::endl;
//    solver_QuadProg sQP;
//    sQP.solveProblem();

//    std::cout<<"==================================================="<<std::endl;
//    std::cout<<"USE CGAL:"<<std::endl;
//    solver_CGAL sCGAL;
//    sCGAL.solveProblem();



#if defined WIN32
    System("pause");
#endif

    return 0;
}


