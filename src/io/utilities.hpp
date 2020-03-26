#ifndef UTILITIES_HPP
#define UTILITIES_HPP

namespace invEHL
{
namespace io
{
class Utilities
{
public:
    static double startTime;
    static double endTime;
    static void waterMark();
    static void tic(bool output = false);
    static void toc(bool output = false);
    static double tictoc(bool output = false);
};
} // namespace io
} // namespace invEHL

#endif