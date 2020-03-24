#ifndef TFE_H
#define TFE_H

#include <eigen3/Eigen/Dense>
#include <fstream>
#include <string>

#include "pde/mesh.hpp"
#include "pde/data.hpp"
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
	const Eigen::VectorXd &one() const { return one_; }

	void initialization(const std::string &iniFIleName);
	void resetData();
	void setFunction(double (*fp)(double x, double y), Eigen::VectorXd &f, double f0 = -1.);
	void setFunction(Eigen::VectorXd &f, double f0);
	void setFunction(const char *filename, Eigen::VectorXd &f, double f0 = -1.);
	void BDF(const Eigen::VectorXd &h0, const Eigen::VectorXd &h1, double dt0, double dt1, Eigen::VectorXd &h2, Flag flag = Flag::BDFINFO_OFF);

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

	// ================= DATA MEMBERS ================== //

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
	void initVector(double (*fp)(double x, double y), Eigen::VectorXd &h) const
	{
		for (size_t i = 0; i < mesh.node.size(); i++)
		{
			h(i) = (*fp)(mesh.node[i][0], mesh.node[i][1]);
		}
	};

	static Eigen::VectorXd ehd(double (*fp)(double x, double y), const Eigen::VectorXd &h, const Eigen::VectorXd &f);
	// ehd function : evaluate function fp which only depends function h, e.g. mobility h^3
	static Eigen::VectorXd ehd(double (*fp)(double x), const Eigen::VectorXd &h);
	static double h3(double h)
	{
		return h * h * h;
	}; // mobility = H^3
	static double dh3dh(double h)
	{
		return 3. * h * h;
	}; // d mobility / d H = 3 * H^2
	// Pi = 1/(D - H)^2 eqn (3.18) where D = 1+ f here :
	static double PI(double h, double f)
	{
		return 1. / pow(1. + f - h, 2.0);
	};
	// Potental = integral of Pi in H
	static double intPIdh(double h, double f)
	{
		return 1. / (1. + f - h);
	};
	static double dPIdh(double h, double f)
	{
		return 2. / pow(1. + f - h, 3.0);
	}; // partial Pi / partial  H
	// partial Pi / partial  D = partial Pi / partial  f since D = 1 + f
	static double dPIdd(double h, double f)
	{
		return -2. / pow(1. + f - h, 3.0);
	};
	double alpha[4] = {0, 0, 0, 0};

private:
	Eigen::SparseMatrix<double> W, dW;
	Eigen::SparseMatrix<double, Eigen::RowMajor> J;
	Eigen::VectorXd F, one_;

	static void bdf(int o, double *alpha);
	static void bdf(int o, double *alpha, const Eigen::VectorXd &h0, const Eigen::VectorXd &h1, Eigen::VectorXd &hBDF);
};
} // namespace pde
} // namespace invEHL
#endif
