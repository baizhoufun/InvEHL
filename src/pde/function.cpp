#include "pde/function.hpp"

namespace invEHL
{
namespace pde
{
Eigen::VectorXd Function::ehd(double (*fp)(double x, double y), const Eigen::VectorXd &h, const Eigen::VectorXd &f)
{
    Eigen::VectorXd tmp = h;
    tmp.setZero();
    for (int i = 0; i < h.size(); i++)
    {
        tmp(i) = (*fp)(h(i), f(i)); // functional pointer takes two variables h and f
    }
    return tmp;
}

// overload of Mesh::ehd method
Eigen::VectorXd Function::ehd(double (*fp)(double x), const Eigen::VectorXd &h)
{
    Eigen::VectorXd tmp = h;
    tmp.setZero();
    for (int i = 0; i < h.size(); i++)
    {
        tmp(i) = (*fp)(h(i)); // functional pointer only takes one variable h
    }
    return tmp;
}
double Function::h3(double h)
{
    return h * h * h;
}; // mobility = H^3
double Function::dh3dh(double h)
{
    return 3. * h * h;
}; // d mobility / d H = 3 * H^2

// Pi = 1/(D - H)^2 eqn (3.18) where D = 1+ f here :
double Function::PI(double h, double f)
{
    return 1. / pow(1. + f - h, 2.0);
};
// Potental = integral of Pi in H
double Function::intPIdh(double h, double f)
{
    return 1. / (1. + f - h);
};
double Function::dPIdh(double h, double f)
{
    return 2. / pow(1. + f - h, 3.0);
}; // partial Pi / partial  H
// partial Pi / partial  D = partial Pi / partial  f since D = 1 + f
double Function::dPIdd(double h, double f)
{
    return -2. / pow(1. + f - h, 3.0);
};

//static Eigen::VectorXd dPIdh(double(*fp)(double x, double y), const Eigen::VectorXd &h, const Eigen::VectorXd &f);
} // namespace pde
} // namespace invEHL
