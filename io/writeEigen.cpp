
#include <eigen3/Eigen/Core>
#include <fstream>

#include "writeEigen.hpp"

namespace invEHL
{
namespace io
{
void IOEigen::write(const std::string &fileName, const Eigen::VectorXd &f)
{
    std::ofstream file;
    file.open(fileName);
    file << f.format(Eigen::FullPrecision);
    file.close();
    printf("Write to %s\n", fileName.c_str());
}

void IOEigen::write(const std::string &fileName, const Eigen::MatrixXd &f)
{
    std::ofstream file;
    file.open(fileName);
    file << f.format(Eigen::FullPrecision);
    file.close();
    printf("Write to %s\n", fileName.c_str());
}
} // namespace io
} // namespace invEHL
