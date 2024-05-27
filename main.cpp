#include <iostream>
#include <cmath>
#include <vector>

#include "eulerpecesolver.h"

#define USE_MATH_DEFINES

using namespace std;

const int a = -2;
const float b = 0.2;
const int d = 1;

double f (double x, double y)
{
    return (a*x - b*y);
}

double runge4(double t, double y0, double step)
{
    double k1, k2, k3, k4;
    k1 = f(t, y0);
    k2 = f(t+step/2, y0+step*k1/2);
    k3 = f(t+step/2, y0+step*k2/2);
    k4 = f(t+step, y0+step*k3);
    double y = y0 +(k1+2*k2+2*k3+k4)*step/6;
    return y;
}

class PECESolver : public EulerPECESolver
{
    virtual Eigen::VectorXd RecalcSystem(double time, Eigen::VectorXd& val) override
    {
        Eigen::VectorXd ret(val.size());
        double y1 = val(0);
        double y2 = val(1);
        ret(0) = 9*y1 + 24*y2 + 5*cos(time) - sin(time)/3;
        ret(1) = -24*y1 - 51*y2 - 9*cos(time) + sin(time)/3;
        return ret;
    }
};


int main()
{
    //task1

    // double initVal = d;
    // double step = 0.0045;
    // double C = d + a/(pow(b, 2));
    // double maxErr = 0;
    // vector<double> time, y;

    // time.push_back(0);
    // y.push_back(initVal);
    // for(double t = 0; t <= 1; t+=step)
    // {
    //     double dy = runge4(time.back(), y.back(), step);
    //     time.push_back(t);
    //     y.push_back(dy);
    // }

    // for (size_t i=0; i<time.size(); i++)
    // {
    //     double u = a/b*(time[i] - 1/b) + C*exp(-b*time[i]);
    //     double err = fabs(u - y[i]);
    //     cout << "u(t) = " << u << ", y(t) = " << y[i] << ", error: " <<  err << endl;
    //     if (err >= maxErr)
    //         maxErr = err;
    // }
    // cout << "Max error: " << maxErr << endl;

    //task2

    PECESolver system;
    Eigen::VectorXd init(2);
    init << 4.0/3, 2.0/3;
    system.InitValues(init, 0);
    double max_Err = 0;
    double min_Err = 10;
    while (system.t() < 5)
    {
        system.CalcStep();
        Eigen::VectorXd vals = system.vals();

        double u1 = 2*exp(-3*system.t()) - exp(-39*system.t()) + cos(system.t())/3;
        double u2 = -exp(-3*system.t()) + 2*exp(-39*system.t()) - cos(system.t())/3;
        double err1 = fabs(u1 - vals(0));
        double err2 = fabs(u2 - vals(1));

        if (err1 > max_Err)
            max_Err = err1;
        if (err2 > max_Err)
            max_Err = err2;
        if (err1 < min_Err)
            min_Err = err1;
        if (err2 < min_Err)
            min_Err = err2;

        cout << "U1 = " << vals(0) << " err1 = " << err1 << "\tU2 = " <<
            vals(1) << " err2 = " << err2  << "\tT = " << system.t() << endl;

    }
    cout << "Max error: " << max_Err << endl;
    cout << "Min error: " << min_Err << endl;
    cout << "step: " << system.m_step << endl;


    return 0;
}
