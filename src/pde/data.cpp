#include "data.hpp"

namespace invEHL
{
namespace pde
{
void Data::clear()
{
	for (auto &it : state())
	{
		it.setZero();
	}
	for (auto &it : adjoint())
	{
		it.setZero();
	}
	for (auto &it : constraint())
	{
		it.setZero();
	}
	rawGradient().setZero();
};

} // namespace pde
} // namespace invEHL