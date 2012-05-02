

#include "genericSolver.h"

#include <iostream>



int main(int argc, char** argv)
{
    int exitVal=0;
    
    GenericSolver solver;

    solver.setProblemFromFile("/home/joe/myProblem.txt");
    solver.printProblem();
    
    std::cout<<"==================================================="<<std::endl;
    std::cout<<"USE QuadProg++:"<<std::endl;
    solver.setSolver("QuadProg");
    solver.solveProblem();
    
    std::cout<<"==================================================="<<std::endl;
    std::cout<<"USE CGAL:"<<std::endl;
    solver.setSolver("CGAL");
    solver.solveProblem();
    
    std::cout<<"==================================================="<<std::endl;
    std::cout<<"USE COPL_QP:"<<std::endl;
    solver.setSolver("COPL_QP");
    solver.solveProblem();


#if defined WIN32
    System("pause");
#endif

    return 0;
}


