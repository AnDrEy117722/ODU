#include "eulerpecesolver.h"
#include "iostream"
using namespace std;
EulerPECESolver::EulerPECESolver() {

    steps.resize(4);
    A.resize(4, 4);
    b1.resize(4);
    b2.resize(4);
    A << 0, 0, 0, 0,
        1.0/2, 0, 0, 0,
        0, 1.0/2, 0, 0,
        0, 0, 1, 0;
    b1 << 1.0/6, 1.0/3, 1.0/3, 1.0/6;
    b2 << 0, 0, 0, 0;
    steps << 0, 1.0/2, 1.0/2, 1;


}

void EulerPECESolver::CalcStep()
{

    // Eigen::VectorXd predictor = m_values + RecalcSystem(m_time, m_values)*m_step;
    // Eigen::VectorXd corrector = m_values + (RecalcSystem(m_time, m_values)
    //                                         + RecalcSystem(m_time+m_step, predictor))*(m_step/2);
    // m_values = corrector;
    // m_time += m_step;

    Eigen::MatrixXd tmp(4, m_values.size());
    Eigen::VectorXd x1, x2;
    double diff = 1;
    do
    {
        Eigen::VectorXd val = m_values;
        tmp.row(0) = RecalcSystem(m_time, val);
        for (int i = 1; i < 4; i++)
        {
            Eigen::VectorXd val = m_values;
            for (int j = 0; j < i; j++)
                val += tmp.row(j)*A(i,j)*m_step;
            tmp.row(i) = RecalcSystem(m_time + steps(i)*m_step, val);
        }
        x1 = b1 * tmp;
        x2 = b2 * tmp;
        diff = (x1-x2).cwiseAbs().maxCoeff()*m_step;
        if (diff > m_max)
            m_step /= 2;
        if (diff < m_min)
            m_step *=2;
    }
    while (diff>m_max);
    m_values = m_values + x1*m_step;
    m_time += m_step;



}
