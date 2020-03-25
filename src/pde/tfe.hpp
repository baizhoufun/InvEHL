#ifndef TFE_H
#define TFE_H

#include <eigen3/Eigen/Dense>
#include <fstream>
#include <string>

#include "eikonal/eikonal.hpp"
#include "pde/mesh.hpp"
#include "pde/data.hpp"
#include "pde/function.hpp"
#include "io/iniReader.hpp"
#include "io/ioEigen.hpp"

namespace invEHL
{
namespace pde
{
class TFE
{
public:
	struct Param
	{
		double tolFRes = 1e-13;
		double tolAdjointSolver = 1e-4;
		double tolStateSolver = 1e-6;
		double dt = 0.01;
		double tMax = 1;
		int tStep = 1;
		bool breakAlarm = false;
		int bdf = 1;
		std::string rootDir;
	} param;

	enum class Flag : unsigned int
	{
		BDFINFO_ON = 1,
		BDFINFO_OFF = 2
	};

	TFE(){};
	TFE(const std::string &iniFIleName) { initialization(iniFIleName); };
	io::INIReader iniReader;
	Mesh mesh;
	Data data;
	image::Eikonal eikonal;

	void initialization(const std::string &iniFIleName);
	void resetData();
	void setFunction(double (*fp)(double x, double y), Eigen::VectorXd &f, double f0 = -1.) const;
	void setFunction(Eigen::VectorXd &f, double f0) const;
	void setFunction(const char *filename, Eigen::VectorXd &f, double f0 = -1.);
	void rescale(double zScale, double zAvg, Eigen::VectorXd &h0) const;
	void BDF(const Eigen::VectorXd &h0, const Eigen::VectorXd &h1, double dt0, double dt1, Eigen::VectorXd &h2, Flag flag = Flag::BDFINFO_OFF);

	const Eigen::VectorXd &one() const { return one_; }

	// ================= DATA MEMBERS ================== //

	struct solverInfo
	{
		int newtonIter = 0;	  // how many inner newton step used for each time step
		double deltaTime = 0; // size of each time step
	};
	enum options
	{
		FORWARD_ONLY,
		FORWARD_AND_BACKWARD,
		FORWARD_AND_BACKWARD_AND_CONSTRAINTFORCE,
		INFO_OFF,
		INFO_ON
	};

	std::vector<solverInfo> stepMonitor;

	// ================= DYNAMICS AND OPTIMIZATION RELATED MEMBER FUNCTIONS ================== //
	solverInfo BDF(int oBDF, int tk, Eigen::VectorXd &hNext);
	void backW(int tk, Eigen::VectorXd &lk, options OPTION = FORWARD_AND_BACKWARD);
	void forward(options OPTION, options SOLVERINFO = INFO_OFF);
	void processRawGradient();
	void updateInvHessian(const Eigen::VectorXd &pk, const Eigen::VectorXd &sk, Eigen::MatrixXd &invHessian);
	double objective(double c0);
	Eigen::MatrixXd computeFreeEnergy();

	// ================= INLINE STATIC FUNCTIONS FOR EHL MODEL ================== //

	// Apply initialized nodal vector h with a function fp evaluated at each FEM nodal point

private:
	double alpha[4] = {0, 0, 0, 0};
	Eigen::SparseMatrix<double> W, dW;
	Eigen::SparseMatrix<double, Eigen::RowMajor> J;
	Eigen::VectorXd F, one_;

	static void bdf(int o, double *alpha);
	static void bdf(int o, double *alpha, const Eigen::VectorXd &h0, const Eigen::VectorXd &h1, Eigen::VectorXd &hBDF);
};
} // namespace pde
} // namespace invEHL
#endif
