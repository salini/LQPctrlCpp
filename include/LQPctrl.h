#ifndef __LQPCTRLCPP_LQPCTRL__
#define __LQPCTRLCPP_LQPCTRL__


#include "Task.h"
#include "Constraint.h"
#include "Solver.h"


struct RobotStateContainer


class LQPctrl
{
public:
    LQPctrl();
    ~LQPctrl();

    void init();
    void update(dt);
    void solve();

    Eigen::VectorXd getSolution();
private:    
    void updateRobotState(double dt);
    void updateTasks(double dt);
    void updateConstraints(double dt);

    RobotStateContainer robotState;
    std::List<Task*> listOfTasks;

    GenericSolver* solver;
};




#endif