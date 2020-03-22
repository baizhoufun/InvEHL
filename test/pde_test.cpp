#include <iostream>
#include <eigen3/Eigen/Core>

#include "../eikonal/eikonal.hpp"
#include "../iniReader/iniReader.hpp"

using namespace cv;

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

    Mat img = imread(reader.GetString("eikonal", "input", "wtf"));
    cv::resize(img, img, cv::Size(200, 200), 0.5, 0.5);
    invEHL::image::Eikonal ls(img, reader.GetBoolean("eikonal", "flip", 0));

    ls.evolution(reader.GetInteger("eikonal", "iterLS", -1),
                 reader.GetReal("eikonal", "dtLS", -1),
                 reader.GetReal("eikonal", "c1", -1),
                 reader.GetReal("eikonal", "c2", -1));
    ls.rescaleMinMax(0, 255);
    cv::imwrite(reader.GetString("eikonal", "outputLS", "wtf"), ls.phi);

    GaussianBlur(ls.phiInit(), ls.phi,
                 cv::Size(reader.GetInteger("eikonal", "kGF", -1), reader.GetInteger("eikonal", "kGF", -1)),
                 reader.GetReal("eikonal", "sigmaGF", -1));

    ls.rescaleMinMax(0, 255);
    cv::imwrite(reader.GetString("eikonal", "outputGF", "wtf"), ls.phi);

    //std::cout << "test Int: " << reader.GetInteger("PDE", "nx", -1) << std::endl;
    //std::cout << "test real: " << reader.GetReal("PDE", "dt", -1) << std::endl;
    //std::cout << "test bool: " << reader.GetBoolean("PDE", "oS", 0) << std::endl;
    //std::cout << "test string " << reader.GetString("path", "rF", "wtf") << std::endl;
    //std::cout << "test section: " << reader.HasSection("PDE") << std::endl;
    //std::cout << "test value: " << reader.HasValue("path", "rF") << std::endl;

    return 0;
}
