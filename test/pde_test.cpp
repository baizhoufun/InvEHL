#include <iostream>
#include <eigen3/Eigen/Core>
#include <thread>

//#include "eikonal/eikonal.hpp"
//#include "io/iniReader.hpp"
#include "io/ioEigen.hpp"
#include "io/utilities.hpp"
//#include "pde/mesh.hpp"
#include "pde/tfe.hpp"

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

    //invEHL::io::Utilities::waterMark();

    invEHL::pde::TFE tfe(confPath);
    tfe.rescale(0.5, 0.0, tfe.data.control());
    invEHL::io::IOEigen::write("../resources/control.txt", tfe.data.control());
    tfe.setFunction(tfe.data.state()[0], 0.15);

    double dt = tfe.param.dt;
    double dtLast = dt;
    double tNow = 0;
    int q = 0;
    double g = 0.1;

    while (tfe.BDF(tfe.data.state()[0], tfe.data.state()[0], tfe.data.state()[1], dt, dt, 1, invEHL::pde::TFE::Flag::BDFINFO_ON) != 0)
    {
        q = q + 1;
        dt = dt / pow(1.0 + g, q);
        std::cout << "REDUCE TIME STEP ...\n";
        std::cout << tNow << " | " << dtLast << " | " << dt;
    }

    tNow += dt;

    for (int i = 2; i < tfe.param.tStep; ++i)
    {
        const Eigen::VectorXd &stateLast = tfe.data.state()[i - 1];
        invEHL::io::IOEigen::mat2Img(stateLast,
                                     tfe.mesh.info.ny * 2, tfe.mesh.info.nx * 2, tfe.mesh.info.ly / tfe.mesh.info.lx,
                                     0, 0.5);
        std::cout << tNow << " | " << dtLast << " | " << dt;
        q = 0;
        while (tfe.BDF(tfe.data.state()[i - 2], tfe.data.state()[i - 1], tfe.data.state()[i], dtLast, dt, tfe.param.bdf, invEHL::pde::TFE::Flag::BDFINFO_ON) != 0)
        {
            q = q + 1;
            dt = dt / pow(1.0 + g, q);
            std::cout << "REDUCE TIME STEP ...\n";
            std::cout << tNow << " | " << dtLast << " | " << dt;
            //tfe.BDF(tfe.data.state()[i - 2], tfe.data.state()[i - 1], tfe.data.state()[i], dtLast, dt, 2, invEHL::pde::TFE::Flag::BDFINFO_ON);
        }
        tNow += dt;
        dtLast = dt;
        //tfe.BDF(tfe.data.state()[i - 1], tfe.data.state()[i - 1], tfe.data.state()[i], tfe.param.dt, tfe.param.dt, 1, invEHL::pde::TFE::Flag::BDFINFO_ON);
    }
    invEHL::io::IOEigen::write("../resources/state/", tfe.data.state());
    return 0;
}
