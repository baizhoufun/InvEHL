#include <iostream>
#include <eigen3/Eigen/Core>
#include <thread>

#include "eikonal/eikonal.hpp"
//#include "io/iniReader.hpp"
#include "io/ioEigen.hpp"
#include "io/utilities.hpp"
//#include "pde/mesh.hpp"
#include "pde/tfe.hpp"

using namespace cv;

class fk
{
public:
    void tk(int a){};
};

int main(int argc, char **argv)
{
    std::string confPath = "../resources/in.ini";
    invEHL::io::INIReader reader(confPath);

    if (reader.ParseError() < 0)
    {
        std::cout << "Can't load 'test.ini'\n";
        return 1;
    };

    if (argc > 1)
        confPath = argv[1];

    Mat img = imread(reader.GetString("eikonal", "input"));
    cv::resize(img, img,
               cv::Size(reader.GetInteger("mesh", "nx") * 2 + 1,
                        reader.GetInteger("mesh", "ny") * 2 + 1));

    invEHL::image::Eikonal ls(img, reader.GetBoolean("eikonal", "flip"));

    ls.evolution(reader.GetInteger("eikonal", "iterLS"),
                 reader.GetReal("eikonal", "dtLS"),
                 reader.GetReal("eikonal", "c1"),
                 reader.GetReal("eikonal", "c2"));

    ls.rescaleMinMax(0, 255);
    cv::imwrite(reader.GetString("eikonal", "outputLS"), ls.phi);

    // GaussianBlur(ls.phiInit(), ls.phi,
    //              cv::Size(reader.GetInteger("eikonal", "kGF"), reader.GetInteger("eikonal", "kGF")),
    //              reader.GetReal("eikonal", "sigmaGF"));

    // ls.rescaleMinMax(0, 255);
    // cv::imwrite(reader.GetString("eikonal", "outputGF"), ls.phi);

    ls.rescaleMinMax();
    invEHL::pde::TFE tfe("../resources/in.ini");

    //invEHL::io::Utilities::waterMark();
    //Eigen::VectorXd b;
    //invEHL::io::IOEigen::img2Mat(ls.phi, b);

    return 0;
}
