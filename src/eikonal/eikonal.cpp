#include <iostream>
#include <opencv2/opencv.hpp>
#include "eikonal.hpp"

using namespace std;
using namespace cv;

namespace invEHL
{
namespace image
{

Eikonal::Eikonal()
{
}

Eikonal::Eikonal(cv::Mat inputImg, bool FLIP, Flag flag)
{
    rec_ = cv::Rect(0, 0, inputImg.rows, inputImg.cols);
    flip = FLIP;
    initializePhi(inputImg, rec_, flag);
}

Eikonal::~Eikonal()
{
}

void Eikonal::initializePhi(cv::Mat inputImg, cv::Rect inputPhi, Flag flag)
{
    cvtColor(inputImg, img, CV_BGR2GRAY);

    col = inputImg.cols;
    row = inputImg.rows;
    dep = CV_32FC1;

    phiInit_ = Mat::zeros(row, col, dep);

    float phiIn = 0.0f;
    float phiOut = 1.0f;
    if (flip)
    {
        phiIn = 1.0f;
        phiOut = 0.0f;
    }

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if (j < inputPhi.y || j > inputPhi.y + inputPhi.height || i < inputPhi.x || i > inputPhi.x + inputPhi.width)
            {
                phiInit_.at<float>(i, j) = 0;
            }
            else
            {
                switch (flag)
                {
                case Flag::INITIAL_ORIGINAL:
                {
                    phiInit_.at<float>(i, j) = static_cast<float>((int)img.at<uchar>(i, j));
                    break;
                }
                case Flag::INITIAL_BINARY:
                {
                    if ((int)img.at<uchar>(i, j) <= 130)
                        phiInit_.at<float>(i, j) = phiIn;
                    else
                        phiInit_.at<float>(i, j) = phiOut;
                    break;
                }

                default:
                    break;
                }
            }
        }
    }
    phi = phiInit_.clone();
}

void Eikonal::evolution(int iterNum, float dt, float c1, float c2)
{
    Mat grad_x, grad_y, lap_xy;
    phi = phiInit_.clone();

    for (int a = 0; a < iterNum; a++)
    {

        Sobel(phi, grad_y, dep, 0, 1, 3);
        Sobel(phi, grad_x, dep, 1, 0, 3);
        Laplacian(phi, lap_xy, dep, 1);
        //GaussianBlur(phi, lap_xy, cv::Size(5, 5), 1);

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                if (j < rec_.y || j > rec_.y + rec_.height || i < rec_.x || i > rec_.x + rec_.width)
                {
                    phi.at<float>(i, j) = 0.0;
                }
                else
                {
                    float dx = grad_x.at<float>(i, j);
                    float dy = grad_y.at<float>(i, j);
                    float p = phi.at<float>(i, j);
                    float dxdy = lap_xy.at<float>(i, j);
                    phi.at<float>(i, j) -= dt * (p - p0) * (sqrt(c1 + dx * dx + dy * dy) - 1.0f);
                    phi.at<float>(i, j) += dt * c2 * dxdy;
                }
            }
        }
    }
    //rescaleMinMax(0.0, 255.0);
    //GaussianBlur(img, phi, cv::Size(13, 13), 5);
}

void Eikonal::gaussianBlur(int kGF, float sigmaGF)
{
    GaussianBlur(phiInit(), phi, cv::Size(kGF, kGF), sigmaGF);
}; // namespace image
} // namespace image
} // namespace invEHL
