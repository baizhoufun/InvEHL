#include <iostream>
#include <eigen3/Eigen/Core>

#include "../eikonal/eikonal.hpp"
#include "../io/iniReader.hpp"
#include "../io/writeEigen.hpp"
#include "../pde/mesh.hpp"

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
    cv::resize(img, img, cv::Size(121, 121), 0.5, 0.5);
    invEHL::image::Eikonal ls(img, reader.GetBoolean("eikonal", "flip"));

    ls.evolution(reader.GetInteger("eikonal", "iterLS"),
                 reader.GetReal("eikonal", "dtLS"),
                 reader.GetReal("eikonal", "c1"),
                 reader.GetReal("eikonal", "c2"));
    ls.rescaleMinMax(0, 255);
    cv::imwrite(reader.GetString("eikonal", "outputLS"), ls.phi);
    ls.rescaleMinMax();

    invEHL::pde::Mesh mesh(
        reader.GetReal("mesh", "lx"),
        reader.GetReal("mesh", "ly"),
        reader.GetInteger("mesh", "nx"),
        reader.GetInteger("mesh", "nx"), 1, 1);

    mesh.initNode();
    mesh.initElement();
    mesh.assembleMass();
    mesh.assembleStiff();

    mesh.outputMesh(
        reader.GetString("mesh", "outputElement"),
        reader.GetString("mesh", "outputNode"));

    Eigen::MatrixXd a(120, 120);
    Eigen::VectorXd b(mesh.node.size());
    std::cout << b.size();
    for (int i = 0; i < 120; i++)
    {
        for (int j = 0; j < 120; j++)
        {
            a(i, j) = ls.phi.at<float>(i, j);
            b(i * 120 + j) = ls.phi.at<float>(i, j);
        }
    }

    Eigen::VectorXd c = mesh.lumpedLaplaceMatrix * b;

    invEHL::io::IOEigen::write("../resources/abc.txt", b);

    return 0;
}
