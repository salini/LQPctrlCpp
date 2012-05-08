
#include "Constraint.h"


Constraint::Constraint()
{

}

Constraint::~Constraint()
{

}

void Constraint::init(unsigned int _nDof, unsigned int _nGForce, unsigned int _nFc, unsigned int _nConstraintMax, bool _useDependentFormalism)
{
    nDof = _nDof;
    nGForce = _nGForce;
    nFc = _nFc;
    nConstraintMax = _nConstraintMax;
    setFormalism(_useDependentFormalism);
    
    if (useDependentFormalism)
    {
        nProblem = nDof + nGForce + nFc;
    }
    else
    {
        nProblem = nGForce + nFc;
    }
}


void Constraint::setFormalism(bool _useDependentFormalism)
{
    useDependentFormalism = _useDependentFormalism;
}

unsigned int Constraint::getnConstraintMax()
{
    return nConstraintMax;
}

unsigned int Constraint::getnConstraint()
{
    return nConstraint;
}

unsigned int Constraint::getnProblem()
{
    return nProblem;
}



//////////////////////////////////////////////////
EqualityConstraint::EqualityConstraint()
{
}

EqualityConstraint::~EqualityConstraint()
{
}

void EqualityConstraint::init(unsigned int _nDof, unsigned int _nGForce, unsigned int _nFc, unsigned int _nConstraintMax, bool _useDependentFormalism)
{
    Constraint::init(_nDof, _nGForce, _nFc, _nConstraintMax, _useDependentFormalism);

    A.setZero(nConstraintMax, nProblem);
    b.setZero(nConstraintMax);
}




//////////////////////////////////////////////////
InequalityConstraint::InequalityConstraint()
{
}

InequalityConstraint::~InequalityConstraint()
{
}

void InequalityConstraint::init(unsigned int _nDof, unsigned int _nGForce, unsigned int _nFc, unsigned int _nConstraintMax, bool _useDependentFormalism)
{
    Constraint::init(_nDof, _nGForce, _nFc, _nConstraintMax, _useDependentFormalism);

    G.setZero(nConstraint, nProblem);
    h.setZero(nConstraint);
}





//////////////////////////////////////////////////
EquationOfMotion::EquationOfMotion()
{

}

EquationOfMotion::~EquationOfMotion()
{

}

void EquationOfMotion::setEquationOfMotion(Eigen::MatrixXd& M, Eigen::MatrixXd& Jchi_T, Eigen::VectorXd& G_minus_N)
{
//    if ( useDependentFormalism )
//    {
//        A = M; //hstack(M, Jchi_T);  //TODO: change
//        b = G_minus_N;
//    }
//    else
//    {
//        // Nothing to do: A and b should be void matrix, or should not be added during constraints concatenation
//    }
}



//////////////////////////////////////////////////
ContactAcceleration::ContactAcceleration()
{

}

ContactAcceleration::~ContactAcceleration()
{

}



void ContactAcceleration::setContactAcceleration(Eigen::MatrixXd& Jc, Eigen::VectorXd& dJc_dot_gvel, bool* constActivity)
{

}



//////////////////////////////////////////////////
GForceMax::GForceMax()
{
}

GForceMax::~GForceMax()
{
}

void GForceMax::setGForceMax(Eigen::VectorXd& _gForceMax, Eigen::VectorXd& _dgForceMax, double dt, Eigen::VectorXd& previousGForce)
{
}




//////////////////////////////////////////////////
FrictionContact::FrictionContact()
{
}

FrictionContact::~FrictionContact()
{
}

void FrictionContact::setFrictionContact(Eigen::VectorXd& mus, bool* constActivity)
{
}

//////////////////////////////////////////////////
JointLimits::JointLimits()
{
}

JointLimits::~JointLimits()
{
}

void JointLimits::setJointLimits(Eigen::VectorXd& qLim, Eigen::VectorXd& vLim, Eigen::VectorXd& aLim, double hPos, double hVel)
{
}




//////////////////////////////////////////////////
CollisionAvoidance::CollisionAvoidance()
{
}

CollisionAvoidance::~CollisionAvoidance()
{
}

void CollisionAvoidance::setCollisionAvoidance(Eigen::VectorXd& sDist, Eigen::VectorXd& sVel, Eigen::MatrixXd& Jca, Eigen::MatrixXd& dJca, double hPos)
{
}





