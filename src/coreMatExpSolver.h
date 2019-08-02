/*
Author: Zack Taylor

Core solver for matrix expotentials. Computes y(t) = exp(A*t)v
the solution to:
	y' + Ay = 0

using the CRAM Method. The code was copied from pyne CRAM solver
*/

#include <Eigen/Core>
#include <complex>


class SolverType {

	private:
	
	Eigen::Matrix<std::complex<double>,8,1> theta;
	Eigen::Matrix<std::complex<double>,8,1> alpha;

	public:
	SolverType(){
		// Sets the values of theta
		theta(0,0) = (-10.843917078696988026, 19.277446167181652284);
		theta(1,0) = (5.2649713434426468895, 16.220221473167927305);
		theta(2,0) = (5.9481522689511774808, 3.5874573620183222829);
		theta(3,0) = (3.5091036084149180974, 8.4361989858843750826);
		theta(4,0) = (6.4161776990994341923, 1.1941223933701386874);
		theta(5,0) = (1.4193758971856659786, 10.925363484496722585);
		theta(6,0) = (4.9931747377179963991, 5.9968817136039422260);
		theta(7,0) = (-1.4139284624888862114, 13.497725698892745389);

		// Sets the values of alpha
		alpha(0,0) = (-.0000005090152186522491565,-.00002422001765285228797);
		alpha(1,0) = (.00021151742182466030907, .0043892969647380673918);
		alpha(2,0) = (113.39775178483930527, 101.9472170421585645);
		alpha(3,0) = (15.059585270023467528, -5.7514052776421819979);
		alpha(4,0) = (-64.500878025539646595, -224.59440762652096056);
		alpha(5,0) = (-1.4793007113557999718, 1.7686588323782937906);
		alpha(6,0) = (-62.518392463207918892, -11.19039109428322848);
		alpha(7,0) = (.041023136835410021273, -.15743466173455468191);
	}

	// Solver function
	Eigen::MatrixXd solve(Eigen::MatrixXd, Eigen::MatrixXd, double);
};
