#include <iostream>
#include "../io/iniReader.hpp"

int main(int argc, char **argv)
{
    std::string confPath = "../resources/in.ini";

    if (argc > 1)
        confPath = argv[1];

    invEHL::io::INIReader reader(confPath);

    if (reader.ParseError() < 0)
    {
        std::cout << "Can't load 'test.ini'\n";
        return 1;
    };

    std::cout << reader.GetInteger("PDE", "xnx", -1) << std::endl;
    std::cout << reader.GetReal("PDE", "c0", -1) << std::endl;
    std::cout << reader.HasSection("PDE") << std::endl;
    std::cout << reader.HasValue("PDE", "xnx") << std::endl;
    std::cout << reader.GetBoolean("PDE", "oS", 0) << std::endl;
    std::cout << std::endl;
    return 0;
}