
#include <eigen3/Eigen/Core>
#include <fstream>

#include "ioEigen.hpp"

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

void IOEigen::img2Mat(const cv::Mat &img, Eigen::VectorXd &b)
{
    Eigen::size_t row = img.rows - 1;
    Eigen::size_t col = img.cols - 1;
    b.resize(row * col);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            b(i * col + j) = static_cast<double>(img.at<float>(i, j));
        }
    }
}
// namespace io
} // namespace io
} // namespace invEHL
