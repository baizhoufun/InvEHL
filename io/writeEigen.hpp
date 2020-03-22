#ifndef WRITEEIGEN_HPP
#define WRITEEIGEN_HPP

#include <eigen3/Eigen/Core>
#include <string>

namespace invEHL
{
namespace io
{
class IOEigen
{
public:
    static void write(const std::string &fileName, const Eigen::VectorXd &f);

    static void write(const std::string &fileName, const Eigen::MatrixXd &f);
};
} // namespace io
} // namespace invEHL

#endif