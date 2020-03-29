#ifndef TFE_H
#define TFE_H

#include <eigen3/Eigen/Dense>
#include <fstream>
#include <string>

#include "eikonal/eikonal.hpp"
#include "pde/mesh.hpp"
#include "pde/data.hpp"
#include "io/iniReader.hpp"

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
		double tolStateError = 1e-10;
		double dt = 0.01;
		double tMax = 1;
		double h0 = 0.125;
		int tStep = 1;
		bool breakAlarm = false;
		int bdf = 1;
		int maxNewtonIter = 6;
		std::string rootDir;
	} param;

	enum class Flag : unsigned int
	{
		BDFINFO_ON = 1,
		BDFINFO_OFF = 3
	};

	TFE(){};
	TFE(const std::string &iniFIleName) { initialization(iniFIleName); };
	io::INIReader iniReader;
	Mesh mesh;
	Data data;
	image::Eikonal eikonal;

	void initialization(const std::string &iniFIleName);
	void resetData();
	//void setFunction(const double (*fp)(double x, double y), Eigen::VectorXd &f, double f0 = -1.) const;
	void setFunction(Eigen::VectorXd &f, double f0) const;
	void setFunction(const char *filename, Eigen::VectorXd &f, double f0 = -1.);
	void rescale(double zScale, double zAvg, Eigen::VectorXd &h0) const;
	int BDF(const Eigen::VectorXd &h0, const Eigen::VectorXd &h1, Eigen::VectorXd &h2, double dt0, double dt1, int bdfOrder, Flag flag = Flag::BDFINFO_OFF);
	double FNewton(double gamma, const Eigen::VectorXd &hBDF, const Eigen::VectorXd &h2, Eigen::VectorXd &F);
	double JNewton(double gamma, const Eigen::VectorXd &h1, Eigen::SparseMatrix<double, Eigen::RowMajor> &J);
	const Eigen::VectorXd &one() const { return one_; }

	// ================= DATA MEMBERS ================== //
	//virtual double fAnalytic(double x, double y) = 0;

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
	void setParam();
	double alpha[4] = {0, 0, 0, 0};
	Eigen::SparseMatrix<double> W_, dW_, WLap_;
	Eigen::SparseMatrix<double, Eigen::RowMajor> J_;
	Eigen::VectorXd F_, one_;

	//static void bdf(int o, double *alpha);
	static void bdf(int o, double *alpha, double delta1 = 1.0);
	static void bdf(int o, double *alpha, const Eigen::VectorXd &h0, const Eigen::VectorXd &h1, Eigen::VectorXd &hBDF);
};
} // namespace pde
} // namespace invEHL
#endif
