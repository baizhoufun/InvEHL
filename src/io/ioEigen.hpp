#ifndef WRITEEIGEN_HPP
#define WRITEEIGEN_HPP

#include <opencv2/opencv.hpp>
#include <eigen3/Eigen/Core>
#include <string>

namespace invEHL
{
namespace io
{
class IOEigen
{
public:
    static Eigen::MatrixXd readMatrix(const char *filename, int MAXBUFFSIZE);

    static void write(const std::string &fileName, const Eigen::VectorXd &f);

    static void write(const std::string &fileName, const Eigen::MatrixXd &f);

    static void write(const std::string &fileName, const std::vector<Eigen::VectorXd> &fContainer, int k = 1);

    static void img2Mat(const cv::Mat &img, Eigen::VectorXd &b);
};
} // namespace io
} // namespace invEHL

#endif