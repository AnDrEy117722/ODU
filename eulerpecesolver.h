#ifndef EULERPECESOLVER_H
#define EULERPECESOLVER_H

#include <Eigen/Eigen>

class EulerPECESolver
{
public:

    EulerPECESolver();

    void InitValues (Eigen::VectorXd const & initVals, double initTime)
    {
        m_values = initVals;
        m_time = initTime;
    }
    void CalcStep();
    double t() const {return m_time;}
    Eigen::VectorXd vals() const {return m_values;}
    double m_step = 0.01;

protected:

    virtual Eigen::VectorXd RecalcSystem(double time, Eigen::VectorXd& val) = 0;
    Eigen::VectorXd m_values, steps;
    double m_time;
    // double m_step = 0.1;
    Eigen::Matrix<double, 4, 4> A;
    Eigen::Matrix<double, 1, 4> b1, b2;

    double m_max = 0.001;
    double m_min = 0.0001;
};


#endif // EULERPECESOLVER_H
