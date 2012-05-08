#ifndef __LQPCTRLCPP_CONSTRAINT__
#define __LQPCTRLCPP_CONSTRAINT__

#include <Eigen/Dense>

class Constraint
{
public:
    Constraint();
    virtual ~Constraint();

    virtual void init(unsigned int _nDof, unsigned int _nGForce, unsigned int _nFc, unsigned int _nConstraintMax, bool _useDependentFormalism=true);
    
    //void setIndependentFormalism(Eigen::MatrixXd& _invM_Jchi_T, Eigen::VectorXd& _invM_G_N);
    
    void setFormalism(bool useDependentFormalism);
    unsigned int getnConstraintMax();
    unsigned int getnConstraint();
    unsigned int getnProblem();

protected:
    unsigned int nDof;
    unsigned int nGForce;
    unsigned int nFc;
    unsigned int nProblem;
    unsigned int nConstraint;
    unsigned int nConstraintMax;

    bool useDependentFormalism;

    Eigen::MatrixXd* invM_Jchi_T;
    Eigen::VectorXd* invM_G_N;
};



class EqualityConstraint: public Constraint
{
public:
    EqualityConstraint();
    virtual ~EqualityConstraint();

    virtual void init(unsigned int _nDof, unsigned int _nGForce, unsigned int _nFc, unsigned int _nConstraintMax, bool _useDependentFormalism=true);
protected:
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
};



class InequalityConstraint: public Constraint
{
public:
    InequalityConstraint();
    virtual ~InequalityConstraint();

    virtual void init(unsigned int _nDof, unsigned int _nGForce, unsigned int _nFc, unsigned int _nConstraintMax, bool _useDependentFormalism=true);

protected:
    Eigen::MatrixXd G;
    Eigen::VectorXd h;
};




////////////////////////////////////////////////
class EquationOfMotion: public EqualityConstraint
{
public:
    EquationOfMotion();
    ~EquationOfMotion();

    void setEquationOfMotion(Eigen::MatrixXd& M, Eigen::MatrixXd& Jchi_T, Eigen::VectorXd& G_minus_N);
};





////////////////////////////////////////////////
class ContactAcceleration: public EqualityConstraint
{
public:
    ContactAcceleration();
    ~ContactAcceleration();

    void setContactAcceleration(Eigen::MatrixXd& Jc, Eigen::VectorXd& dJc_dot_gvel, bool* constActivity);
};




class GForceMax: public InequalityConstraint
{
    GForceMax();
    ~GForceMax();

    void setGForceMax(Eigen::VectorXd& _gForceMax, Eigen::VectorXd& _dgForceMax, double dt, Eigen::VectorXd& previousGForce);
};


class FrictionContact: public InequalityConstraint
{
    FrictionContact();
    ~FrictionContact();

    void setFrictionContact(Eigen::VectorXd& mus, bool* constActivity);
};


class JointLimits: public InequalityConstraint
{
    JointLimits();
    ~JointLimits();

    void setJointLimits(Eigen::VectorXd& qLim, Eigen::VectorXd& vLim, Eigen::VectorXd& aLim, double hPos, double hVel);
};


class CollisionAvoidance: public InequalityConstraint
{
    CollisionAvoidance();
    ~CollisionAvoidance();

    void setCollisionAvoidance(Eigen::VectorXd& sDist, Eigen::VectorXd& sVel, Eigen::MatrixXd& Jca, Eigen::MatrixXd& dJca, double hPos);
};



#endif

