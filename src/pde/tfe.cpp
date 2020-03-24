#include "tfe.hpp"
#include <iostream>
namespace invEHL
{
namespace pde
{
// overload of Mesh::ehd method
Eigen::VectorXd TFE::ehd(double (*fp)(double x, double y), const Eigen::VectorXd &h, const Eigen::VectorXd &f)
{
    Eigen::VectorXd tmp = h;
    tmp.setZero();
    for (int i = 0; i < h.size(); i++)
    {
        tmp(i) = (*fp)(h(i), f(i)); // functional pointer takes two variables h and f
    }
    return tmp;
}

// overload of Mesh::ehd method
Eigen::VectorXd TFE::ehd(double (*fp)(double x), const Eigen::VectorXd &h)
{
    Eigen::VectorXd tmp = h;
    tmp.setZero();
    for (int i = 0; i < h.size(); i++)
    {
        tmp(i) = (*fp)(h(i)); // functional pointer only takes one variable h
    }
    return tmp;
}

void TFE::resetData()
{
    const int &tStep = param.tStep;
    const int &dof = mesh.info.dof;

    data.time().setZero(tStep);

    data.state().resize(tStep);
    for (auto &it : data.state())
    {
        it.setZero(dof);
    }
    data.adjoint().resize(tStep);
    for (auto &it : data.adjoint())
    {
        it.setZero(dof);
    }
    data.constraint().resize(tStep);
    for (auto &it : data.adjoint())
    {
        it.setZero(dof);
    }
    data.control().setZero(dof);
    data.target().setZero(dof);
    data.rawGradient().setZero(dof);
    //F.setZero(fem.info.dof);
    one_.setOnes(dof);
}
void TFE::initialization(const std::string &iniFIleName)
{
    iniReader.Read(iniFIleName);
    while (iniReader.ParseError() < 0)
    {
        std::cout << "Can't load.\n";
        return;
    }

    mesh.info.lx = iniReader.GetReal("mesh", "lx");
    mesh.info.ly = iniReader.GetReal("mesh", "ly");
    mesh.info.nx = iniReader.GetInteger("mesh", "nx");
    mesh.info.ny = iniReader.GetInteger("mesh", "ny"),
    param.dt = iniReader.GetReal("pde", "dt");
    param.tStep = iniReader.GetInteger("pde", "tStep");
    param.rootDir = iniReader.GetString("pde", "rootDir");
    param.tolFRes = iniReader.GetReal("pde", "tolFRes");
    param.tolAdjointSolver = iniReader.GetReal("pde", "tolAdjointSolver");
    param.tolStateSolver = iniReader.GetReal("pde", "tolStateSolver");

    mesh.initNode();
    mesh.initElement();
    mesh.assembleMass();
    mesh.assembleStiff();
    mesh.outputMesh(iniReader.GetString("mesh", "outputElement"), iniReader.GetString("mesh", "outputNode"));
    resetData();
}
void TFE::setFunction(double (*fp)(double x, double y), Eigen::VectorXd &f, double f0)
{
    initVector(fp, f);
    if (f0 > -0.5)
    {
        double favg = one().dot(mesh.lumpedMassMatrix * f) / mesh.oneMassOne();
        f = f + (f0 - favg) * one();
        std::cout << "Analytic input loaded; avg from " << favg << " to "
                  << one().dot(mesh.lumpedMassMatrix * f) / mesh.oneMassOne() << " !\n";
    }
};

void TFE::setFunction(Eigen::VectorXd &f, double f0) { f.setConstant(f0); };

void TFE::setFunction(const char *filename, Eigen::VectorXd &f, double f0)
{
    const int &dof = mesh.info.dof;
    Eigen::MatrixXd tt = io::IOEigen::readMatrix(filename, 5 * dof);
    while (tt.rows() != dof)
    {
        printf("Input DOF (%d) and Mesh DOF (%d) do not match! Press enter to reload ...", tt.rows(), dof);
        std::cin.get();
        tt.resize(0, 0);
        tt = io::IOEigen::readMatrix(filename, 5 * dof);
    }
    f = Eigen::VectorXd(tt.col(0));

    if (f0 > -0.5)
    {
        double favg = one().dot(mesh.lumpedMassMatrix * f) / mesh.oneMassOne();
        f = f + (f0 - favg) * one();
        std::cout << "User input loaded; avg from " << favg << " to "
                  << one().dot(mesh.lumpedMassMatrix * f) / mesh.oneMassOne() << " !\n";
    }
    else
    {
        std::cout << "User input loaded and unmodified ! \n";
    }
};

void TFE::BDF(const Eigen::VectorXd &h0, const Eigen::VectorXd &h1, double dt0, double dt1, Eigen::VectorXd &h2, Flag flag)
{
    const Eigen::VectorXd &f = data.control();
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> sol;
    double F_res = 1;
    int maxIter = 15;
    const int &dof = mesh.info.dof;
    const double &tolStateSolver = param.tolStateSolver;
    const double &tolFRes = param.tolFRes;

    W.setZero();
    mesh.assembleWeightedStiff(W, ehd(h3, h1));
    Eigen::SparseMatrix<double> dPI(dof, dof);
    dPI.setIdentity();
    Eigen::SparseMatrix<double> tmp(dof, dof);
    tmp = W * mesh.lumpedLaplaceMatrix;

    h2 = h1;
    dPI.diagonal() = ehd(dPIdh, h1, f);
    Eigen::VectorXd hBDF(h1.size());
    bdf(1, alpha, h0, h1, hBDF);

    sol.compute(-J);
    sol.setTolerance(tolStateSolver);
    int iter = 0;

    // nonlinear equation (3.44) that we need to solve
    F = mesh.lumpedMassMatrix * (h2 - hBDF) - dt1 * alpha[2] * (tmp * h2 + W * ehd(PI, h2, f));
    F_res = F.cwiseAbs().maxCoeff();
    while (F_res > tolFRes && iter < maxIter)
    {
        iter++;
        h2 += sol.solve(F);
        if (sol.info() != Eigen::Success)
        {
            printf("State solver failure ...\n");
            //breakAlarm = true;
        }
        F = mesh.lumpedMassMatrix * (h2 - hBDF) - dt1 * alpha[2] * (tmp * h2 + W * ehd(PI, h2, f));
        F_res = F.cwiseAbs().maxCoeff();
    }
    if (flag == Flag::BDFINFO_ON)
    {
        std::cout << "BDF nonlinear iter step = " << iter;
    }
}
void TFE::bdf(int o, double *alpha)
{
    switch (o)
    {
    case 1:
    {
        alpha[2] = 1.0;
        alpha[1] = 1.0;
        alpha[0] = 0;

        //hBDF = h_0;
        break;
    }
    case 2:
    {
        alpha[2] = 2. / 3.;
        alpha[1] = 4. / 3.;
        alpha[0] = -1. / 3.;
        //hBDF = a_zero * h1 + a_minus * h0;
        break;
    }
    default:
        break;
    }
};
void TFE::bdf(int o, double *alpha, const Eigen::VectorXd &h0, const Eigen::VectorXd &h1, Eigen::VectorXd &hBDF)
{
    bdf(1, alpha);
    switch (o)
    {
    case 1:
    {
        hBDF = h1;
        break;
    }
    case 2:
    {
        alpha[2] = 2. / 3.;
        alpha[1] = 4. / 3.;
        alpha[0] = -1. / 3.;
        hBDF = alpha[1] * h1 + alpha[0] * h0;
        break;
    }
    default:
        break;
    }
};

// void thinFilm::initMesh(double lx_, double ly_, int nx_, int ny_, double dt_, int tStep_)
// {
// 	//waterMark();
// 	// init info
// 	fem.info.nx = nx_;
// 	fem.info.ny = ny_;
// 	fem.info.lx = lx_;
// 	fem.info.ly = ly_;
// 	fem.info.dt = dt_;
// 	fem.info.tStep = tStep_;
// 	// init Mesh
// 	fem.initNode();
// 	fem.initElement();
// 	// assemble state-independent matrix
// 	fem.assembleMass();
// 	fem.assembleStiff();
// 	// allocate state container
// 	state.resize(fem.info.tStep);
// 	for (auto &it : state)
// 	{
// 		it.setZero(fem.info.dof);
// 	}
// 	// allocate control vector
// 	control.setZero(fem.info.dof);
// 	F.setZero(fem.info.dof);
// 	// allocate stepMonitor
// 	timeStamp.setZero(fem.info.tStep);
// 	stepMonitor.resize(fem.info.tStep);
// 	stepMonitor[0].deltaTime = fem.info.dt;
// 	// allocate adjoint container
// 	adjoint.resize(fem.info.tStep);
// 	for (auto &it : adjoint)
// 	{reader
// 		it.setZero(fem.info.dof);
// 	}
// 	// allocate constraint force container
// 	constraintforce.resize(fem.info.tStep);
// 	for (auto &it : constraintforce)
// 	{
// 		it.setZero(fem.info.dof);
// 	}
// 	// allocate target vector
// 	target.setZero(fem.info.dof);
// 	// allocate gradient dJdD
// 	rawGradient.setZero(fem.info.dof);
// 	one.setOnes(fem.info.dof);
// }

// //set function f from functional pointer fp and normalized it to value f0
// void thinFilm::setFunction(double (*fp)(double x, double y), Eigen::VectorXd &f, double f0)
// {
// 	fem.initVector(fp, f);
// 	if (f0 > -0.5)
// 	{
// 		double favg = one.dot(fem.lumpedMassMatrix * f) / Mesh::oneMassOne();
// 		f = f + (f0 - favg) * one;
// 		std::cout << "Analytic input loaded; avg from " << favg << " to "
// 				  << one.dot(fem.lumpedMassMatrix * f) / Mesh::oneMassOne() << " !\n";
// 	}
// };

// //set function f to const value f0
// void thinFilm::setFunction(double f0, Eigen::VectorXd &f) { f.setConstant(f0); };
//setFun
// //set function f from text file input and normalize to value f0
// void thinFilm::setFunction(const char *filename, Eigen::VectorXd &f, double f0)
// {
// 	Eigen::MatrixXd tt = readMatrix(filename, 5 * fem.info.dof);
// 	while (tt.rows() != fem.info.dof)
// 	{
// 		printf("Input DOF (%d) and Mesh DOF (%d) do not match! Press enter to reload ...", tt.rows(), fem.info.dof);
// 		std::cin.get();
// 		tt.resize(0, 0);
// 		tt = readMatrix(filename, 5 * fem.info.dof);
// 	}
// 	f = Eigen::VectorXd(tt.col(0));

// 	if (f0 > -0.5)
// 	{
// 		double favg = one.dot(fem.lumpedMassMatrix * f) / Mesh::oneMassOne();
// 		f = f + (f0 - favg) * one;
// 		std::cout << "User input loaded; avg from " << favg << " to "
// 				  << one.dot(fem.lumpedMassMatrix * f) / Mesh::oneMassOne() << " !\n";
// 	}
// 	else
// 	{
// 		std::cout << "User input loaded! \n";
// 	}
// };

// thinFilm::solverInfo thinFilm::BDF(int oBDF, int tk, Eigen::VectorXd &hNext)
// {
// 	const Eigen::VectorXd &h_0 = state[tk];
// 	const Eigen::VectorXd &f = control;
// 	const double &dtLast = stepMonitor[((tk - 1) < 0 ? 0 : (tk - 1))].deltaTime;
// 	double dt = dtLast;
// 	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> sol;
// 	double F_res = 1;
// 	int maxIteration = 15;

// 	// ---------------------------
// 	W.setZero();
// 	fem.assembleWeightedStiff(W, Mesh::ehd(h3, h_0));
// 	Eigen::SparseMatrix<double> dPI(fem.info.dof, fem.info.dof);
// 	dPI.setIdentity();
// 	Eigen::SparseMatrix<double> tmp(fem.info.dof, fem.info.dof);
// 	tmp = W * fem.lumpedLaplaceMatrix;
// 	// ---------------------------

// 	hNext = h_0;
// 	dPI.diagonal() = Mesh::ehd(dPIdh, h_0, f);
// 	Eigen::VectorXd hBDF(h_0.size());

// 	//BDF1 AND BDF2 coefficients - table 3.1
// 	switch (oBDF)
// 	{
// 	case 1:
// 	{
// 		a_plus = 1.;
// 		a_zero = 1.;
// 		a_minus = 0.;
// 		hBDF = h_0;
// 		break;
// 	}
// 	case 2:
// 	{
// 		a_plus = 2. / 3.;
// 		a_zero = 4. / 3.;
// 		a_minus = -1. / 3.;
// 		hBDF = a_zero * h_0 + a_minus * state[tk - 1];
// 		break;
// 	}
// 	default:
// 		break;
// 	}

// 	// newton iteration jacobian - equation (3.47)
// 	J = fem.lumpedMassMatrix - dt * a_plus * (tmp + W * dPI);
// 	sol.compute(-J);
// 	sol.setTolerance(tolStateSolver);
// 	int iteration = 0;
// 	// nonlinear equation (3.44) that we need to solve
// 	F = fem.lumpedMassMatrix * (hNext - hBDF) - dt * a_plus * (tmp * hNext + W * Mesh::ehd(PI, hNext, f));
// 	F_res = F.cwiseAbs().maxCoeff();
// 	// newton iteration at each time step - equation (3.48)
// 	while (F_res > tolFRes && iteration < maxIteration)
// 	{
// 		iteration++;
// 		//hNext += sol.solveWithGuess(F,hNext);
// 		hNext += sol.solve(F);
// 		if (sol.info() != Eigen::Success)
// 		{
// 			printf("State solver failure ...\n");
// 			breakAlarm = true;
// 		}
// 		F = fem.lumpedMassMatrix * (hNext - hBDF) - dt * a_plus * (tmp * hNext + W * Mesh::ehd(PI, hNext, f));
// 		F_res = F.cwiseAbs().maxCoeff();
// 	}

// 	thinFilm::solverInfo solverInfoNow;
// 	solverInfoNow.deltaTime = dt;
// 	solverInfoNow.newtonIter = iteration;
// 	timeStamp(tk + 1) = timeStamp(tk) + dt;

// 	return solverInfoNow;
// };
// 	hNext = h_0;
// 	dPI.diagonal() = Mesh::ehd(dPIdh, h_0, f);
// 	Eigen::VectorXd hBDF(h_0.size());

// 	//BDF1 AND BDF2 coefficients - table 3.1
// 	switch (oBDF)
// 	{
// 	case 1:
// 	{
// 		a_plus = 1.;
// 		a_zero = 1.;
// 		a_minus = 0.;
// 		hBDF = h_0;
// 		break;
// 	}
// 	case 2:
// 	{
// 		a_plus = 2. / 3.;
// 		a_zero = 4. / 3.;
// 		a_minus = -1. / 3.;
// 		hBDF = a_zero * h_0 + a_minus * state[tk - 1];
// 		break;
// 	}
// 	default:
// 		break;
// 	}

// 	// newton iteration jacobian - equation (3.47)
// 	J = fem.lumpedMassMatrix - dt * a_plus * (tmp + W * dPI);
// 	sol.compute(-J);
// 	sol.setTolerance(tolStateSolver);
// 	int iteration = 0;
// 	// nonlinear equation (3.44) that we need to solve
// 	F = fem.lumpedMassMatrix * (hNext - hBDF) - dt * a_plus * (tmp * hNext + W * Mesh::ehd(PI, hNext, f));
// 	F_res = F.cwiseAbs().maxCoeff();
// 	// newton iteration at each time step - equation (3.48)
// 	while (F_res > tolFRes && iteration < maxIteration)
// 	{
// 		iteration++;
// 		//hNext += sol.solveWithGuess(F,hNext);
// 		hNext += sol.solve(F);
// 		if (sol.info() != Eigen::Success)
// 		{
// 			printf("State solver failure ...\n");
// 			breakAlarm = true;
// 		}
// 		F = fem.lumpedMassMatrix * (hNext - hBDF) - dt * a_plus * (tmp * hNext + W * Mesh::ehd(PI, hNext, f));
// 		F_res = F.cwiseAbs().maxCoeff();
// 	}

// 	thinFilm::solverInfo solverInfoNow;
// 	solverInfoNow.deltaTime = dt;
// 	solverInfoNow.newtonIter = iteration;
// 	timeStamp(tk + 1) = timeStamp(tk) + dt;

// 	return solverInfoNow;
// };

// void thinFilm::backW(int tk, Eigen::VectorXd &lk, options OPTION)
// {
// 	int tk1 = clamp(tk + 1, 0, fem.info.tStep - 1);
// 	int tk2 = clamp(tk + 2, 0, fem.info.tStep - 1);
// 	int tk_1 = clamp(tk - 1, 0, fem.info.tStep - 1);

// 	const Eigen::VectorXd &lk1 = adjoint[tk1];
// 	const Eigen::VectorXd &lk2 = adjoint[tk2];
// 	const Eigen::VectorXd &f = control;
// 	const Eigen::VectorXd &hk1 = state[tk1];
// 	const Eigen::VectorXd &hk = state[tk];
// 	const Eigen::VectorXd &hk_1 = state[tk_1];
// 	const double &dtk = stepMonitor[tk].deltaTime;
// 	const double &dtk_1 = stepMonitor[tk_1].deltaTime;

// 	double a_plus_k_1 = 2. / 3.;
// 	double a_plus_k = 2. / 3.;
// 	double a_zero_k = 4. / 3.;
// 	double a_minus_k1 = -1. / 3.;

// 	if (tk == 1)
// 	{
// 		a_plus_k_1 = 1;
// 	}

// 	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> sol;
// 	Eigen::VectorXd lBDF(lk1.size());
// 	Eigen::SparseMatrix<double> dPI(fem.info.dof, fem.info.dof);
// 	dPI.setIdentity();
// 	W.setZero();
// 	dW.setZero();

// 	switch (fem.info.tStep - 1 - tk)
// 	{
// 	case 0:
// 	{
// 		lBDF = fem.lumpedMassMatrix * (hk - target);
// 		//lBDF =  fem.lumpedLaplaceMatrix* (hk - target);
// 		a_plus_k = 0.0;
// 		break;
// 	}
// 	case 1:
// 	{
// 		lBDF = a_zero_k * fem.lumpedMassMatrix * lk1;
// 		break;
// 	}
// 	default:
// 		lBDF = fem.lumpedMassMatrix * (a_zero_k * lk1 + a_minus_k1 * lk2);
// 		break;
// 	}

// 	Eigen::VectorXd pk1 = -fem.lumpedLaplaceMatrix * hk1 - Mesh::ehd(PI, hk1, f);
// 	fem.assembleWeightedStiff(W, Mesh::ehd(h3, hk_1), dW, Mesh::ehd(dh3dh, hk), pk1);

// 	dPI.diagonal() = Mesh::ehd(dPIdh, hk, f); // index k

// 	//Eigen::SparseMatrix<double> LHS = fem.lumpedMassMatrix + dtk_1 * a_plus_k_1 * Eigen::SparseMatrix<double>((W * (-fem.lumpedLaplaceMatrix - dPI)).transpose());
// 	Eigen::SparseMatrix<double> LHS = (W * (-fem.lumpedLaplaceMatrix - dPI)).transpose();
// 	LHS *= dtk_1 * a_plus_k_1;
// 	LHS += fem.lumpedMassMatrix;
// 	Eigen::VectorXd RHS = lBDF - dtk * a_plus_k * Eigen::SparseMatrix<double>(dW.transpose()) * lk1;

// 	sol.setTolerance(tolAdjointSolver);
// 	sol.compute(LHS);
// 	lk = sol.solve(RHS);

// 	if (OPTION == FORWARD_AND_BACKWARD_AND_CONSTRAINTFORCE)
// 	{
// 		Eigen::VectorXd tmp = -dtk_1 * a_plus_k_1 * dPI * W * lk;
// 		constraintforce[tk] = fem.inverseMassMatrix * tmp;
// 		rawGradient += tmp;
// 	}
// 	else
// 	{
// 		rawGradient += -dtk_1 * a_plus_k_1 * dPI * W * lk;
// 	}
// 	//lk = sol.solveWithGuess(RHS,lk1);
// 	if (sol.info() != Eigen::Success)
// 	{
// 		printf("Adjoint solver failure ...\n");
// 	}

// 	//printf("Time %5.5f --", omp_get_wtime() - start_time);
// }

// void thinFilm::forward(options OPTION, options SOLVERINFO)
// {
// 	int counter = 0;
// 	clearData();

// 	while (counter < fem.info.tStep - 1)
// 	{
// 		tic();
// 		if (counter > 0)
// 		{
// 			stepMonitor[counter] = BDF(2, counter, state[counter + 1]);
// 			if (breakAlarm)
// 				return;
// 		}
// 		else
// 		{
// 			stepMonitor[counter] = BDF(1, counter, state[counter + 1]);
// 			if (breakAlarm)
// 				return;
// 		}
// 		toc();
// 		if (SOLVERINFO == options::INFO_ON)
// 		{
// 			printf("\rStep %03d\tInnerIt %02d\tTime %5.5f", counter, stepMonitor[counter].newtonIter, tictoc());
// 			fflush(stdout);
// 		}
// 		counter++;
// 	}

// 	if (OPTION == FORWARD_AND_BACKWARD || OPTION == FORWARD_AND_BACKWARD_AND_CONSTRAINTFORCE)
// 	{
// 		rawGradient.setZero();
// 		counter = fem.info.tStep - 1;
// 		while (counter > 0)
// 		{
// 			backW(counter, adjoint[counter], OPTION);
// 			if (SOLVERINFO == options::INFO_ON)
// 			{
// 				printf("\rStep %03d", counter);
// 				fflush(stdout);
// 			}
// 			counter--;
// 		}
// 	}
// 	printf("\n");
// }

// // project raw gradient with mass matrix
// void thinFilm::processRawGradient()
// {
// 	rawGradient = fem.inverseMassMatrix * rawGradient;
// 	double favg = one.dot(fem.lumpedMassMatrix * rawGradient) / Mesh::oneMassOne();
// 	rawGradient += (0. - favg) * one;
// };

// //BFGS update of Hessian inverse
// void thinFilm::updateInvHessian(const Eigen::VectorXd &sk, const Eigen::VectorXd &yk, Eigen::MatrixXd &invHessian)
// {

// 	double skyk = sk.dot(fem.lumpedMassMatrix * yk);
// 	double ykInvHyk = yk.dot(fem.lumpedMassMatrix * invHessian * yk);
// 	Eigen::MatrixXd tmp = invHessian;
// 	invHessian.noalias() += -1. / skyk *
// 							(tmp * yk * sk.transpose() * fem.lumpedMassMatrix + sk * yk.transpose() * fem.lumpedMassMatrix * tmp);
// 	invHessian.noalias() += (skyk + ykInvHyk) / (skyk * skyk) * sk * sk.transpose() * fem.lumpedMassMatrix;
// };

// double thinFilm::objective(double c0)
// {
// 	const Eigen::VectorXd hDiff = state[state.size() - 1] - target;
// 	return 0.5 * hDiff.dot(fem.lumpedMassMatrix * hDiff) + c0 * 0.5 * control.dot(fem.stiffnessMatrix * control);
// };

// Eigen::MatrixXd thinFilm::computeFreeEnergy()
// {
// 	Eigen::MatrixXd freeEnergy(state.size(), 4);
// 	freeEnergy.setZero();
// 	Eigen::VectorXd ee;
// 	ee.setZero(fem.info.dof);
// 	for (int k = 0; k < state.size(); k++)
// 	{
// 		ee.setZero();
// 		const Eigen::VectorXd &hk = state[k];
// 		ee = Mesh::ehd(intPIdh, hk, control);
// 		freeEnergy(k, 3) = timeStamp[k];
// 		freeEnergy(k, 0) = 0.5 * hk.dot(fem.stiffnessMatrix * hk) - one.dot(fem.lumpedMassMatrix * ee);
// 		freeEnergy(k, 1) = hk.maxCoeff();
// 		freeEnergy(k, 2) = hk.minCoeff();
// 	}
// 	printf("Computing free energy for the current mask ...\n");
// 	return freeEnergy;
// };

} // namespace pde
} // namespace invEHL