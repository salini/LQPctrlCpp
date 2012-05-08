
#include "LQPctrl.h"

LQPctrl::LQPctrl()
{
}

LQPctrl::~LQPctrl()
{
}

void LQPctrl::init()
{
}

void LQPctrl::update(double dt)
{
    updateRobotState(dt);
    updateTasks(dt);
    updateConstraints(dt);
}

void LQPctrl::solve()
{

    setProblemConstraints();  // compute A, b, G, h

    sortedListOfTasks = getSortedListOfTasks();
    for (int i=0; i<sortedListOfTasks.size(); i++)
    {
        setCostFunction(sortedListOfTasks[i]); // compute current E, f and consequently P, q
        
        solver->solve(P, q, A, b, G, h);

        if i<levelMax-1:
        {
            updateProblemConstraints(E, x);
        }
    }
}

Eigen::VectorXd LQPctrl::getSolution()
{
    return solver->getSolution();
}