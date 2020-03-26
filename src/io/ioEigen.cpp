#include <eigen3/Eigen/Core>
#include <fstream>
#include <omp.h>

#include "io/ioEigen.hpp"

#define MAXBUFSIZE ((int)1e6)

namespace invEHL
{
namespace io
{
void IOEigen::waterMark()
{
    std::cout << "  ________   ________   _________   ___   _____ ______    ________   ___           \n";
    std::cout << " |\\   __  \\ |\\   __  \\ |\\___   ___\\|\\  \\ |\\   _ \\  _   \\ |\\   __  \\ |\\  \\          \n";
    std::cout << " \\ \\  \\|\\  \\\\ \\  \\|\\  \\\\|___ \\  \\_|\\ \\  \\\\ \\  \\\\\\__\\ \\  \\\\ \\  \\|\\  \\\\ \\  \\         \n";
    std::cout << "  \\ \\  \\\\\\  \\\\ \\   ____\\    \\ \\  \\  \\ \\  \\\\ \\  \\\\|__| \\  \\\\ \\   __  \\\\ \\  \\        \n";
    std::cout << "   \\ \\  \\\\\\  \\\\ \\  \\___|     \\ \\  \\  \\ \\  \\\\ \\  \\    \\ \\  \\\\ \\  \\ \\  \\\\ \\  \\____   \n";
    std::cout << "    \\ \\_______\\\\ \\__\\         \\ \\__\\  \\ \\__\\\\ \\__\\    \\ \\__\\\\ \\__\\ \\__\\\\ \\_______\\ \n";
    std::cout << "     \\|_______| \\|__|          \\|__|   \\|__| \\|__|     \\|__| \\|__|\\|__| \\|_______| \n";
    std::cout << "                               _______    ___  ___   ___                        \n";
    std::cout << "                              |\\  ___ \\  |\\  \\|\\  \\ |\\  \\                       \n";
    std::cout << "                              \\ \\   __/| \\ \\  \\\\\\  \\\\ \\  \\                      \n";
    std::cout << "                               \\ \\  \\_|/__\\ \\   __  \\\\ \\  \\                     \n";
    std::cout << "                                \\ \\  \\_|\\ \\\\ \\  \\ \\  \\\\ \\  \\____                \n";
    std::cout << "                                 \\ \\_______\\\\ \\__\\ \\__\\\\ \\_______\\              \n";
    std::cout << "                                  \\|_______| \\|__|\\|__| \\|_______|              \n";
    std::cout << std::endl;
}
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

void IOEigen::img2Mat(const cv::Mat &inputImg, Eigen::VectorXd &b)
{
    int dep = CV_32FC1;
    cv::Mat img = cv::Mat::zeros(inputImg.rows, inputImg.cols, dep);
    cvtColor(inputImg, img, CV_BGR2GRAY);

    Eigen::size_t row = img.rows - 1;
    Eigen::size_t col = img.cols - 1;
    b.resize(row * col);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            b(i * col + j) = static_cast<float>(img.at<uchar>(i, j)) / 255.0;
        }
    }
}
void IOEigen::mat2Img(const Eigen::VectorXd &b, int row, int col, float aspectRatio, float bmin, float bmax)
{
    int dep = CV_8UC1;
    cv::Mat img(row, col, dep);

//Eigen::size_t row = img.rows - 1;
//Eigen::size_t col = img.cols - 1;
//b.resize(row * col);
#pragma omp parallel for num_threads(4)
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        { //
            img.at<uchar>(i, j) = 255.0 * (b(i * col + j) - bmin) / (bmax - bmin);
            //= static_cast<float>() / 255.0;
        }
    }

    cv::resize(img, img, cv::Size(500, 500 * aspectRatio));

    cv::imshow("", img);
    cv::waitKey(10);
};

void IOEigen::write(const std::string &fileName, const std::vector<Eigen::VectorXd> &fContainer, int k)
{
    int dataLength = fContainer.size();
    unsigned int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        dataLength /= 10;
    } while (dataLength);

#pragma omp parallel for num_threads(4) std::cout << data.control().maxCoeff();
    std::cout << data.control().minCoeff();
    for (int i = 0; i < fContainer.size(); i += k)
    {
        int number_of_zeros = number_of_digits - std::to_string(i).length();
        std::ofstream localfile(fileName + std::string(number_of_zeros, '0').append(std::to_string(i)) + ".txt");
        localfile << fContainer[i].format(Eigen::FullPrecision);
        localfile.close();
    }
    printf("Write to %s %d --%d\n", fileName.c_str(), 0, fContainer.size());
}

Eigen::MatrixXd IOEigen::readMatrix(const char *filename, int MAXBUFFSIZE)
{
    int cols = 0, rows = 0;
    double *buff = new double[MAXBUFFSIZE];
    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(filename);
    while (!infile.eof())
    {
        std::string line;
        std::getline(infile, line);

        int temp_cols = 0;
        std::stringstream stream(line);
        while (!stream.eof())
            stream >> buff[cols * rows + temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    infile.close();
    //rows--;
    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i, j) = buff[cols * i + j];

    delete[] buff;
    return result;
};
// namespace io
} // namespace io
} // namespace invEHL
