#ifndef UTILITIES_HPP
#define UTILITIES_HPP

namespace invEHL
{
namespace io
{
struct Utilities
{
    static double startTime;
    static double endTime;
    static void waterMark();
    static void tic(bool output = false);
    static void toc(bool output = false);
    static double tictoc(bool output = false);
    template <typename T>
    T clamp(const T &value, const T &low, const T &high)
    {
        return value < low ? low : (value > high ? high : value);
    }
};

} // namespace io
} // namespace invEHL
#endif