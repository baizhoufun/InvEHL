#ifndef DATA_H
#define DATA_H
#include <vector>
#include <eigen3/Eigen/Core>

namespace invEHL
{
namespace pde
{
class Data
{
public:
	Data(){};

	const std::vector<Eigen::VectorXd> &state() const { return _state; };
	const std::vector<Eigen::VectorXd> &adjoint() const { return _adjoint; };
	const std::vector<Eigen::VectorXd> &constraint() const { return _constraint; };
	const Eigen::VectorXd &rawGradient() const { return _rawGradient; };
	std::vector<Eigen::VectorXd> &state() { return _state; };
	std::vector<Eigen::VectorXd> &adjoint() { return _adjoint; };
	std::vector<Eigen::VectorXd> &constraint() { return _constraint; };
	Eigen::VectorXd &rawGradient() { return _rawGradient; };

	const Eigen::VectorXd &control() const { return _control; };
	const Eigen::VectorXd &target() const { return _target; };
	const std::vector<double> &time() const { return _time; };

	Eigen::VectorXd &control() { return _control; };
	Eigen::VectorXd &target() { return _target; };
	std::vector<double> &time() { return _time; };
	void clear();

private:
	std::vector<Eigen::VectorXd> _state;
	std::vector<Eigen::VectorXd> _adjoint;
	std::vector<Eigen::VectorXd> _constraint;
	Eigen::VectorXd _control;
	Eigen::VectorXd _rawGradient;
	Eigen::VectorXd _target;
	std::vector<double> _time;
};

//static Eigen::VectorXd dPIdh(double(*fp)(double x, double y), const Eigen::VectorXd &h, const Eigen::VectorXd &f);
} // namespace pde
} // namespace invEHL
#endif