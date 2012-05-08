#ifndef __LQPCTRLCPP_TASK__
#define __LQPCTRLCPP_TASK__

#include <Eigen/Dense>

#include "LQPctrl.h"


enum
{
    COST_NORMAL
} CostType;

enum
{
    NORM_NORMAL;
} NormType;




class Task
{
public:
    Task();
    ~Task();
    
protected:
    Eigen::MatrixXd E;
    Eigen::VectorXd f;

    Eigen::MatrixXd J;
    double          error;

    unsigned int* cdof;
    double*       weight;
    double*       level;
    bool          isActive;

    CostType      cost;
    NormType      norm;
    FormalismType formalism;

//    self._inv_lambda    = zeros((0, 0))
//    self._lambda        = zeros((0, 0))
//    self._inv_ellipsoid = zeros((0, 0))
//    self._ellipsoid     = zeros((0, 0))
//    self._L_T           = zeros((0, 0))

//    self._inv_lambda_is_already_computed    = False
//    self._lambda_is_already_computed        = False
//    self._inv_ellipsoid_is_already_computed = False
//    self._ellipsoid_is_already_computed     = False
};


//////////////////////////////////////////////////////
class AccelerationTask: public Task
{
};

class ForceTask: public Task
{
};


//////////////////////////////////////////////////////
class JointTask: public AccelerationTask
{
};

class MultiJointTask: public AccelerationTask
{
};

class FrameTask: public AccelerationTask
{
};


class comTask: public AccelerationTask
{
};


class LQPcomTask: public AccelerationTask
{
};



//////////////////////////////////////////////////////
class TorqueTask: public ForceTask
{
};

class MultiTorqueTask: public ForceTask
{
};

class WrenchTask: public ForceTask
{
};




//////////////////////////////////////////////////////
class MultiTask: public Task
{
};





#endif


