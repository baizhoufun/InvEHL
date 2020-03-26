#ifndef MESH_H
#define MESH_H
#include <vector>
#include <eigen3/Eigen/Sparse>

namespace invEHL
{
namespace pde
{
class Mesh
{
public:
	struct properties
	{
		double lx = 1, ly = 1;			// domain size x and y
		double joc = 1, jx = 1, jy = 1; //coordiante jacobian or seperately in x and y direction
		int nx = 1, ny = 1, dof = 1;	// number of elements in x and y and total nodal degree of freedom
	} info;

	std::vector<std::vector<double>> node;			  // index container of nodal points
	std::vector<std::vector<int>> element;			  // index container of elements
	Eigen::SparseMatrix<double> consistentMassMatrix; // not-lumped mass matrix
	Eigen::SparseMatrix<double> lumpedMassMatrix;	  // lumped mass matrix
	Eigen::SparseMatrix<double> inverseMassMatrix;	  // inverse of mass matrix
	Eigen::SparseMatrix<double> stiffnessMatrix;	  // stiffness matrix
	Eigen::SparseMatrix<double> lumpedLaplaceMatrix;  // laplacian with lumped mass matrix

	Mesh(){};
	Mesh(double lx_, double ly_, int nx_, int ny_)
	{
		info.lx = lx_;
		info.ly = ly_;
		info.nx = nx_;
		info.ny = ny_;
		//info.dt = dt_;
		//info.tStep = tStep_;
	}
	void initNode();
	void initElement();

	const double &oneMassOne() const { return oneMassOne_; } // one vector dot mass matrix dot one vector = total area of mesh
	void assembleMass();
	void assembleStiff();
	// just assemble weighted stiff matrix
	void assembleWeightedStiff(Eigen::SparseMatrix<double> &matrix, const Eigen::VectorXd &h);
	// assemble weighted stiff then contract it with another one-dimensional vector
	void assembleWeightedStiff(
		Eigen::SparseMatrix<double> &W, const Eigen::VectorXd &H3,
		Eigen::SparseMatrix<double> &dW, const Eigen::VectorXd &dH3dH, const Eigen::VectorXd &P);
	// write mesh to a path
	void outputMesh(const std::string &elementPath, const std::string &nodePath) const;
	// ehd function : evaluate function fp which depends on function f and function h, e.g. Pi = 1/ (1+f-h)
	void initVector(double (*fp)(double x, double y), Eigen::VectorXd &h) const;
	void allocSparseStiff(Eigen::SparseMatrix<double> &W) const;

private:
	const static double localMassMatrix[9][9];		// equation (3.111)
	const static double gaussQuadrature09[9][3];	//
	const static double localStiffnessMatrix[9][9]; // equation (3.111)
	const static double basisAtQuadrature[9][9];	// evluate 9 Lagrange basis at each of the 9 quadrature abscissa
	const static double stiffAtQuadrature[9][9][9];
	const static double stiffDxAtQuadrature[9][9][9];
	const static double stiffDyAtQuadrature[9][9][9];
	const static double weightedStiffDx[9][9][9]; // tensor N1ijk equation (3.114)
	const static double weightedStiffDy[9][9][9]; // tensor N2ijk equation (3.114)
	double oneMassOne_ = 0.0;
};

//static Eigen::VectorXd dPIdh(double(*fp)(double x, double y), const Eigen::VectorXd &h, const Eigen::VectorXd &f);
} // namespace pde
} // namespace invEHL
#endif