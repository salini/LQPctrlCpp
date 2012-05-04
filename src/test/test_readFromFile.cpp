

#include "genericSolver.h"

#include <iostream>

void showResult(GenericSolver& s)
{
    std::cout<<"Cost:\n"<<s.getCost()<<std::endl<<std::endl
             <<"X:\n"<<s.getSolution()<<std::endl;
}



int main(int argc, char** argv)
{
    std::cout<<"------------------ test set problem from file ------------------"<<std::endl;
    GenericSolver solver;
    
    if (argc==1)
    {
        std::cout<<"No argument file. Quit program."<<std::endl;
        return 1;
    }
    solver.setProblemFromFile(argv[1]);
    solver.printProblem();
    
    solver.solveProblem();
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);


#if defined WIN32
    System("pause");
#endif

    return 0;
}


