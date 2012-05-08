

#include "Solver.h"

#include <Eigen/Dense>



void showResult(GenericSolver& s)
{
    std::cout<<"Cost:\n"<<s.getCost()<<std::endl<<std::endl
             <<"X:\n"<<s.getSolution()<<std::endl;
}


int main(int argc, char** argv)
{

    Eigen::MatrixXd P, A, G;
    Eigen::VectorXd q, b, h;
    
    GenericSolver solver;
    
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
    P.resize(3,3);
    q.resize(3);
    P << 1,0,0, 0,1,0, 0,0,1;
    q << 1,2,3;
    solver.setCostFunction(P, q);
    solver.printProblem();
    solver.solveProblem();
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);
    
    std::cout<<"------------------ test problem  with equality constraint ------------------"<<std::endl;
    solver.initProblem(3,2,0);
    A.resize(2,3);
    b.resize(2);
    A << 1,0,0, 0,1,1;
    b << 7,8;
    solver.setCostFunction(P, q);
    solver.setEqualityConstraint(A, b);
    solver.printProblem();
    solver.solveProblem();
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);

    std::cout<<"------------------ test problem  with inequality constraint ------------------"<<std::endl;
    solver.initProblem(3,0,2);
    G.resize(2,3);
    h.resize(2);
    G << 1,0,0, 0,-1,-1;
    h << -7,8;
    solver.setCostFunction(P, q);
    solver.setInequalityConstraint(G, h);
    solver.printProblem();
    solver.solveProblem();
    std::cout<<"Problem solved? "<<solver.problemIsSolved()<<std::endl;
    showResult(solver);

    std::cout<<"------------------ test problem  with inequality constraint ------------------"<<std::endl;
    solver.initProblem(3,1,2);
    A.resize(1,3);
    b.resize(1);
    A << -1,1,-1;
    b << -5;
    solver.setCostFunction(P, q);
    solver.setEqualityConstraint(A, b);
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


