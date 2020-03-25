#ifndef FUNCTION_H
#define FUNCTION_H

#include <eigen3/Eigen/Core>

namespace invEHL
{
namespace pde
{
struct Function
{
public:
    static Eigen::VectorXd ehd(double (*fp)(double x, double y), const Eigen::VectorXd &h, const Eigen::VectorXd &f);
    static Eigen::VectorXd ehd(double (*fp)(double x), const Eigen::VectorXd &h);
    static double h3(double h);
    static double dh3dh(double h);
    static double PI(double h, double f);
    static double intPIdh(double h, double f);
    static double dPIdh(double h, double f);
    static double dPIdd(double h, double f);
};

} // namespace pde
} // namespace invEHL
#endif