/*
Author: Zack Taylor

Core solver for matrix expotentials. Computes y(t) = exp(A*t)v
the solution to:
	y' + Ay = 0

using the CRAM Method. Returns the solution to the system.
The code was copied from pyne CRAM solver
*/

#include <Eigen/Core>
#include <Eigen/Dense>
#include <complex>


class SolverType {

	private:
	
	Eigen::Matrix<std::complex<double>,8,1> theta;
	Eigen::Matrix<std::complex<double>,8,1> alpha;
	double alpha_0 = 2.1248537104952237488e-16;

	public:
	SolverType(){
		std::complex<double> t1(-10.843917078696988026, 19.277446167181652284);
		std::complex<double> t2(-5.2649713434426468895, 16.220221473167927305);
		std::complex<double> t3(5.9481522689511774808, 3.5874573620183222829);
		std::complex<double> t4(3.5091036084149180974, 8.4361989858843750826);
		std::complex<double> t5(6.4161776990994341923, 1.1941223933701386874);
		std::complex<double> t6(1.4193758971856659786, 10.925363484496722585);
		std::complex<double> t7(4.9931747377179963991, 5.9968817136039422260);
		std::complex<double> t8(-1.4139284624888862114, 13.497725698892745389);

		// Sets the values of theta
		theta(0,0) = t1; theta(1,0) = t2; theta(2,0) = t3; theta(3,0) = t4;
		theta(4,0) = t5; theta(5,0) = t6; theta(6,0) = t7; theta(7,0) = t8;

		std::complex<double> a1(-.0000005090152186522491565,-.00002422001765285228797);
		std::complex<double> a2(.00021151742182466030907, .0043892969647380673918);
		std::complex<double> a3(113.39775178483930527, 101.9472170421585645);
		std::complex<double> a4(15.059585270023467528, -5.7514052776421819979);
		std::complex<double> a5(-64.500878025539646595, -224.59440762652096056);
		std::complex<double> a6(-1.4793007113557999718, 1.7686588323782937906);
		std::complex<double> a7(-62.518392463207918892, -11.19039109428322848);
		std::complex<double> a8(.041023136835410021273, -.15743466173455468191);

		// Sets the values of alpha
		alpha(0,0) = a1; alpha(1,0) = a2; alpha(2,0) = a3; alpha(3,0) = a4; 
		alpha(4,0) = a5; alpha(5,0) = a6; alpha(6,0) = a7; alpha(7,0) = a8; 
	}

	// Solver function
	Eigen::MatrixXd solve(Eigen::MatrixXd, Eigen::MatrixXd, double);
};
