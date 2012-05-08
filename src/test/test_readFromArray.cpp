

#include "Solver.h"


void showResult(GenericSolver& s)
{
    std::cout<<"Cost:\n"<<s.getCost()<<std::endl<<std::endl
             <<"X:\n"<<s.getSolution()<<std::endl;
}


int main(int argc, char** argv)
{
    GenericSolver solver;
    
    // WARNING: ARRAY MUST BE INSTANTIATE COLUMN-WISE //
    
    std::cout<<"------------------ test non initialized problem ------------------"<<std::endl;
    solver.printProblem();
    solver.solveProblem();  // should fail without crash
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);
    
    std::cout<<"------------------ test void problem ------------------"<<std::endl;
    solver.initProblem(3,2,1);
    solver.printProblem();
    solver.solveProblem();  // should fail without crash
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);
    
    std::cout<<"------------------ test unconstrained problem ------------------"<<std::endl;
    solver.initProblem(3,0,0);
    double P[]={1,0,0, 0,1,0, 0,0,1};
    double q[]={1,2,3};
    solver.setCostFunction(P, q);
    solver.printProblem();
    solver.solveProblem();
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);
    
    std::cout<<"------------------ test problem  with equality constraint ------------------"<<std::endl;
    solver.initProblem(3,2,0);
    double A[]={1,0, 0,1, 0,1};
    double b[]={7,8};
    solver.setCostFunction(P, q);
    solver.setEqualityConstraint(A, b);
    solver.printProblem();
    solver.solveProblem();
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);

    std::cout<<"------------------ test problem  with inequality constraint ------------------"<<std::endl;
    solver.initProblem(3,0,2);
    double G[]={1,0, 0,-1, 0,-1};
    double h[]={-7,8};
    solver.setCostFunction(P, q);
    solver.setInequalityConstraint(G, h);
    solver.printProblem();
    solver.solveProblem();
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);

    std::cout<<"------------------ test problem  with inequality constraint ------------------"<<std::endl;
    solver.initProblem(3,1,2);
    double A2[]={-1,1,-1};
    double b2[]={-5};
    solver.setCostFunction(P, q);
    solver.setEqualityConstraint(A2, b2);
    solver.setInequalityConstraint(G, h);
    solver.printProblem();
    solver.solveProblem();
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);

    std::cout<<"------------------ test void problem (again) ------------------"<<std::endl;
    solver.initProblem(4,3,3);
    solver.printProblem();
    solver.solveProblem();  // should fail without crash
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);

#if defined WIN32
    System("pause");
#endif

    return 0;
}


