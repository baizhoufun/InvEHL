#include <iostream>

#include "eikonal/eikonal.hpp"
#include "io/iniReader.hpp"

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

    Mat img = imread(reader.GetString("eikonal", "input"));
    cv::resize(img, img, cv::Size(121, 151), 0.5, 0.5);
    invEHL::image::Eikonal ls(img, reader.GetBoolean("eikonal", "flip"));

    ls.evolution(reader.GetInteger("eikonal", "iterLS"),
                 reader.GetReal("eikonal", "dtLS"),
                 reader.GetReal("eikonal", "c1"),
                 reader.GetReal("eikonal", "c2"));
    ls.rescaleMinMax(0, 255);
    cv::imwrite(reader.GetString("eikonal", "outputLS"), ls.phi);

    GaussianBlur(ls.phiInit(), ls.phi,
                 cv::Size(reader.GetInteger("eikonal", "kGF"), reader.GetInteger("eikonal", "kGF")),
                 reader.GetReal("eikonal", "sigmaGF"));

    ls.rescaleMinMax(0, 255);
    cv::imwrite(reader.GetString("eikonal", "outputGF"), ls.phi);

    //std::cout << "test Int: " << reader.GetInteger("PDE", "nx", -1) << std::endl;
    //std::cout << "test real: " << reader.GetReal("PDE", "dt", -1) << std::endl;
    //std::cout << "test bool: " << reader.GetBoolean("PDE", "oS", 0) << std::endl;
    //std::cout << "test string " << reader.GetString("path", "rF", "wtf") << std::endl;
    //std::cout << "test section: " << reader.HasSection("PDE") << std::endl;
    //std::cout << "test value: " << reader.HasValue("path", "rF") << std::endl;

    return 0;
}
