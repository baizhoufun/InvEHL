#include <vector>
#include <opencv2/opencv.hpp>

namespace invEHL
{
namespace image
{
class Eikonal
{
public:
    Eikonal();
    Eikonal(cv::Mat inputImg, bool FLIP = false);
    ~Eikonal();

    void initializePhi(cv::Mat inputImg, cv::Rect inputPhi);
    void rescaleMinMax(float phiMin = 0.0f, float phiMax = 1.0f)
    {
        cv::normalize(phi, phi, phiMax, phiMin, cv::NORM_MINMAX);
    };
    void displayPhi()
    {
    }

    const cv::Mat &phiInit() { return phiInit_; };

    //  input phi, iteration Max, dt, c1, c2 (lap smoothing)
    void evolution(int iterMax = 100, float dt = 0.01f, float c1 = 0.05f, float c2 = 0.01f);
    cv::Mat phi;

private:
    cv::Rect rec_;
    cv::Mat phiInit_;
    cv::Mat img;
    int col;
    int row;
    int dep;
    bool flip;
    float p0 = 0.0001;
};
} // namespace image
} // namespace invEHL
