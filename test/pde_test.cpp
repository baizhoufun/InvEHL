#include <iostream>
#include <eigen3/Eigen/Core>
#include <thread>

//#include "eikonal/eikonal.hpp"
//#include "io/iniReader.hpp"
//#include "io/ioEigen.hpp"
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

    // ls.rescaleMinMax(0, 255);
    // cv::imwrite(reader.GetString("eikonal", "outputGF"), ls.phi);

    invEHL::io::Utilities::waterMark();

    invEHL::pde::TFE tfe("../resources/in.ini");
    tfe.rescale(0.5, 0.0, tfe.data.control());
    invEHL::io::IOEigen::write("../resources/control.txt", tfe.data.control());
    tfe.setFunction(tfe.data.state()[0], 0.125);

    for (int i = 1; i < tfe.param.tStep; ++i)
    {
        tfe.BDF(tfe.data.state()[i - 1], tfe.data.state()[i - 1], tfe.param.dt, tfe.param.dt, tfe.data.state()[i], invEHL::pde::TFE::Flag::BDFINFO_ON);
    }
    invEHL::io::IOEigen::write("../resources/state/", tfe.data.state());

    //Eigen::VectorXd b;
    //invEHL::io::IOEigen::img2Mat(ls.phi, b);

    return 0;
}
