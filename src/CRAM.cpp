#include "CRAM.h"
#define MTAG1 1

//*************************************************************************
// Initialization of solver
//*************************************************************************
SolverType::SolverType(){
	// Method 0 == CRAM
	// Method 1 == Hyperbolic
	// Method 2 == Parabilic
	int order = 16;
	int method = 0;
	MatrixCLD theta_m16(8,1);
	MatrixCLD alpha_m16(8,1);
	MatrixCLD theta_m32(16,1);
	MatrixCLD alpha_m32(16,1);

	// Defines the complex values for CRAM of order 16
	std::complex<long double> tCRAM1(-10.843917078696988026, 19.277446167181652284);
	std::complex<long double> tCRAM2(-5.2649713434426468895, 16.220221473167927305);
	std::complex<long double> tCRAM3(5.9481522689511774808, 3.5874573620183222829);
	std::complex<long double> tCRAM4(3.5091036084149180974, 8.4361989858843750826);
	std::complex<long double> tCRAM5(6.4161776990994341923, 1.1941223933701386874);
	std::complex<long double> tCRAM6(1.4193758971856659786, 10.925363484496722585);
	std::complex<long double> tCRAM7(4.9931747377179963991, 5.9968817136039422260);
	std::complex<long double> tCRAM8(-1.4139284624888862114, 13.497725698892745389);

	// Defines the complex values for parabolic quadratures of order 16
	std::complex<long double> tP16QUAD1(2.020748077156870, 0.785398163397448);
	std::complex<long double> tP16QUAD2(1.431532694411836, 2.356194490192345);
	std::complex<long double> tP16QUAD3(0.253101928921766, 3.926990816987241);
	std::complex<long double> tP16QUAD4(-1.514544219313338, 5.497787143782138);
	std::complex<long double> tP16QUAD5(-3.871405750293477, 7.068583470577035);
	std::complex<long double> tP16QUAD6(-6.817482664018648, 8.639379797371930);
	std::complex<long double> tP16QUAD7(-10.352774960488860, 10.210176124166829);
	std::complex<long double> tP16QUAD8(-14.477282639704097, 11.780972450961723);

	// Defines the complex values for parabolic quadratures of order 32
	std::complex<long double> tP32QUAD1(4.151974038578435,0.785398163397448);
	std::complex<long double> tP32QUAD2(3.857366347205918,2.356194490192345);
	std::complex<long double> tP32QUAD3(3.268150964460883,3.926990816987241);
	std::complex<long double> tP32QUAD4(2.384327890343331,5.497787143782138);
	std::complex<long double> tP32QUAD5(1.205897124853261,7.068583470577035);
	std::complex<long double> tP32QUAD6(-0.267141332009325,8.639379797371930);
	std::complex<long double> tP32QUAD7(-2.034787480244431,10.210176124166829);
	std::complex<long double> tP32QUAD8(-4.097041319852049,11.780972450961723);
	std::complex<long double> tP32QUAD9(-6.453902850832191,13.351768777756622);
	std::complex<long double> tP32QUAD10(-9.105372073184846,14.922565104551516);
	std::complex<long double> tP32QUAD11(-12.051448986910021,16.493361431346415);
	std::complex<long double> tP32QUAD12(-15.292133592007705,18.064157758141310);
	std::complex<long double> tP32QUAD13(-18.827425888477920,19.634954084936208);
	std::complex<long double> tP32QUAD14(-22.657325876320640, 21.205750411731103);
	std::complex<long double> tP32QUAD15(-26.781833555535890,22.776546738526000);
	std::complex<long double> tP32QUAD16(-31.200948926123644,24.347343065320896);

	// Defines the complex values for hyperbolic quadratures of order 16
	std::complex<long double> tH16QUAD1(2.742831168621752,0.943848360919165);
	std::complex<long double> tH16QUAD2(2.135110463494564,2.848825632088136);
	std::complex<long double> tH16QUAD3(0.908542533283144,4.805960934046875);
	std::complex<long double> tH16QUAD4(-0.959329373015139,6.851086689176557);
	std::complex<long double> tH16QUAD5(-3.502703389291920,9.021646302543301);
	std::complex<long double> tH16QUAD6(-6.768145152543949,11.357379698510957);
	std::complex<long double> tH16QUAD7(-10.815440354851331,13.901050903424766);
	std::complex<long double> tH16QUAD8(-15.718689336157446,16.699230995391616);

	// Defines the complex values for hyperbolic quadratures of order 32
	std::complex<long double> tH32QUAD1(5.599231098854665,0.943309422231338);
	std::complex<long double> tH32QUAD2(5.296236884569082,2.834241007838773);
	std::complex<long double> tH32QUAD3(4.688863188784654,4.738130534415102);
	std::complex<long double> tH32QUAD4(3.774333143735294,6.663682444535401);
	std::complex<long double> tH32QUAD5(2.548465585467303,8.619700219595530);
	std::complex<long double> tH32QUAD6(1.005655937866020,10.615126628680125);
	std::complex<long double> tH32QUAD7(-0.861149411053042,12.659084614243690);
	std::complex<long double> tH32QUAD8(-3.060485357738284,14.760919001530999);
	std::complex<long double> tH32QUAD9(-5.602407104019437,16.930239222429300);
	std::complex<long double> tH32QUAD10(-8.498536128749128,19.176963249082540);
	std::complex<long double> tH32QUAD11(-11.762113320347574,21.511362938128986);
	std::complex<long double> tH32QUAD12(-15.408059513168350,23.944110992872538);
	std::complex<long double> tH32QUAD13(-19.453043704453300, 26.486329758095610);
	std::complex<long double> tH32QUAD14(-23.915559263759143, 29.149642070599510);
	std::complex<long double> tH32QUAD15(-28.816008483279890,31.946224397957653);
	std::complex<long double> tH32QUAD16(-34.176795855622940,34.888862508427690);

	// Sets the values of theta
	// CRAM
	if (method == 0){
		theta_m16(0,0) = tCRAM1; theta_m16(1,0) = tCRAM2; theta_m16(2,0) = tCRAM3; 
		theta_m16(3,0) = tCRAM4; theta_m16(4,0) = tCRAM5; theta_m16(5,0) = tCRAM6; 
		theta_m16(6,0) = tCRAM7; theta_m16(7,0) = tCRAM8;
	}
	// Quad hyperbolic
	else if (method == 1){
		theta_m16(0,0) = tH16QUAD1; theta_m16(1,0) = tH16QUAD2; theta_m16(2,0) = tH16QUAD3; 
		theta_m16(3,0) = tH16QUAD4; theta_m16(4,0) = tH16QUAD5; theta_m16(5,0) = tH16QUAD6; 
		theta_m16(6,0) = tH16QUAD7; theta_m16(7,0) = tH16QUAD8;
	}
	// Quad parabolic
	else if (method == 2){
		theta_m16(0,0) = tP16QUAD1; theta_m16(1,0) = tP16QUAD2; theta_m16(2,0) = tP16QUAD3; 
		theta_m16(3,0) = tP16QUAD4; theta_m16(4,0) = tP16QUAD5; theta_m16(5,0) = tP16QUAD6; 
		theta_m16(6,0) = tP16QUAD7; theta_m16(7,0) = tP16QUAD8;
	}
	// Quad hyperbolic
	if (method == 1){
		theta_m32(0,0) = tH32QUAD1; theta_m32(1,0) = tH32QUAD2; theta_m32(2,0) = tH32QUAD3; 
		theta_m32(3,0) = tH32QUAD4; theta_m32(4,0) = tH32QUAD5; theta_m32(5,0) = tH32QUAD6; 
		theta_m32(6,0) = tH32QUAD7; theta_m32(7,0) = tH32QUAD8; theta_m32(8,0) = tH32QUAD9; 
		theta_m32(9,0) = tH32QUAD10; theta_m32(10,0) = tH32QUAD11; theta_m32(11,0) = tH32QUAD12; 
		theta_m32(12,0) = tH32QUAD13; theta_m32(13,0) = tH32QUAD14; theta_m32(14,0) = tH32QUAD15; 
		theta_m32(15,0) = tH32QUAD16;
	}
	// Quad parabolic
	else if (method == 2){
		theta_m32(0,0) = tP32QUAD1; theta_m32(1,0) = tP32QUAD2; theta_m32(2,0) = tP32QUAD3; 
		theta_m32(3,0) = tP32QUAD4; theta_m32(4,0) = tP32QUAD5; theta_m32(5,0) = tP32QUAD6; 
		theta_m32(6,0) = tP32QUAD7; theta_m32(7,0) = tP32QUAD8; theta_m32(8,0) = tP32QUAD9; 
		theta_m32(9,0) = tP32QUAD10; theta_m32(10,0) = tP32QUAD11; theta_m32(11,0) = tP32QUAD12; 
		theta_m32(12,0) = tP32QUAD13; theta_m32(13,0) = tP32QUAD14; theta_m32(14,0) = tP32QUAD15; 
		theta_m32(15,0) = tP32QUAD16;
	}
	if (order == 16){
		theta = theta_m16;
	}
	else{
		theta = theta_m32;
	}

	// Defines the complex values for CRAM of order 16
	std::complex<long double> aCRAM1(-.0000005090152186522491565,-.00002422001765285228797);
	std::complex<long double> aCRAM2(.00021151742182466030907, .0043892969647380673918);
	std::complex<long double> aCRAM3(113.39775178483930527, 101.9472170421585645);
	std::complex<long double> aCRAM4(15.059585270023467528, -5.7514052776421819979);
	std::complex<long double> aCRAM5(-64.500878025539646595, -224.59440762652096056);
	std::complex<long double> aCRAM6(-1.4793007113557999718, 1.7686588323782937906);
	std::complex<long double> aCRAM7(-62.518392463207918892, -11.19039109428322848);
	std::complex<long double> aCRAM8(.041023136835410021273, -.15743466173455468191);

	// Defines the complex values for parabolic quadratures of order 16
	std::complex<long double> aP16QUAD1(-1.083477123084597, -1.583717738032942);
	std::complex<long double> aP16QUAD2(1.156101716894417, -0.323557679282414);
	std::complex<long double> aP16QUAD3(0.014170253999115, 0.441211659436057);
	std::complex<long double> aP16QUAD4(-0.089912139403344, -0.012162781789109);
	std::complex<long double> aP16QUAD5(0.002533198593825, -0.009897386337554);
	std::complex<long double> aP16QUAD6(5.926370534115512e-04, 2.056825502975817e-04);
	std::complex<long double> aP16QUAD7(-8.111258746289452e-06, 1.939107353658974e-05);
	std::complex<long double> aP16QUAD8(-3.477936036617725e-07, -1.653825957153982e-07);

	// Defines the complex values for parabolic quadratures of order 32
	std::complex<long double> aP32QUAD1(-10.182155502279953,-12.289466475490308);
	std::complex<long double> aP32QUAD2(10.723061492746975,-6.014335828841986);
	std::complex<long double> aP32QUAD3(2.465786018329411,6.819492524857861);
	std::complex<long double> aP32QUAD4(-3.177606102879742,0.659072785269153);
	std::complex<long double> aP32QUAD5(-0.092107516054386,-1.088673802832851);
	std::complex<long double> aP32QUAD6(0.274936752551191,0.004268698811926);
	std::complex<long double> aP32QUAD7(-0.005062447831339,0.051274795498324);
	std::complex<long double> aP32QUAD8(-0.007071558977973,-0.001194871857247);
	std::complex<long double> aP32QUAD9(1.653728209519868e-04,-7.219951339349447e-04);
	std::complex<long double> aP32QUAD10(5.461738811469894e-05,1.534908276125415e-05);
	std::complex<long double> aP32QUAD11(-1.000020090193583e-06,3.063390568371581e-06);
	std::complex<long double> aP32QUAD12(-1.274650608438357e-07,-4.671073963453931e-08);
	std::complex<long double> aP32QUAD13(1.582379724603487e-09,-3.936387183948117e-09);
	std::complex<long double> aP32QUAD14(9.025889160226451e-11,3.914918313349633e-11);
	std::complex<long double> aP32QUAD15(-7.106072486810770e-13,1.537125556911646e-12);
	std::complex<long double> aP32QUAD16(-1.944807234220184e-14,-9.492756162090556e-15);

	// Defines the complex values for hyperbolic quadratures of order 16
	std::complex<long double> aH16QUAD1(-2.135226899337539,-4.223725041937824);
	std::complex<long double> aH16QUAD2(2.836500326105343, 0.430540283727702);
	std::complex<long double> aH16QUAD3(-0.679943439267376,0.727328177347865);
	std::complex<long double> aH16QUAD4(-0.035908611305510,-0.181967489706694);
	std::complex<long double> aH16QUAD5(0.015348256628229,0.008542412336624);
	std::complex<long double> aH16QUAD6(-7.814476875623792e-04,1.801696997202231e-04);
	std::complex<long double> aH16QUAD7(1.187651730045994e-05,-1.161329890556223e-05);
	std::complex<long double> aH16QUAD8(-6.818445192754204e-08,1.280624286142434e-07);

	// Defines the complex values for hyperbolic quadratures of order 32
	std::complex<long double> aH32QUAD1(-42.384516271092950,-69.528306342879490);
	std::complex<long double> aH32QUAD2(61.773787127086870,-4.449430324218566);
	std::complex<long double> aH32QUAD3(-14.005943988790632,32.759074033995510);
	std::complex<long double> aH32QUAD4(-9.731768938063905,-11.879046341793380);
	std::complex<long double> aH32QUAD5(4.814098721661006,- 0.945326116980532);
	std::complex<long double> aH32QUAD6(-0.362115384751030,1.090501699188015);
	std::complex<long double> aH32QUAD7(-0.126088783670633,-0.148975837842977);
	std::complex<long double> aH32QUAD8(0.023638139561172,-0.002599527144529);
	std::complex<long double> aH32QUAD9(-0.001056888057661,0.001762268138188);
	std::complex<long double> aH32QUAD10(-3.818818068971985e-05,-1.183739799310812e-04);
	std::complex<long double> aH32QUAD11(4.597704920914048e-06,2.435079803696621e-06);
	std::complex<long double> aH32QUAD12(-1.454295038512129e-07,2.827065808637995e-08);
	std::complex<long double> aH32QUAD13(2.035093130313425e-09,-1.958860964202722e-09);
	std::complex<long double> aH32QUAD14(-1.205237419834894e-11,3.328813345697270e-11);
	std::complex<long double> aH32QUAD15(5.812062995639875e-15,-2.857853355407124e-13);
	std::complex<long double> aH32QUAD16(2.364377732512525e-16,1.434728826471849e-15);

	// Sets the values of alpha
	if (method == 0){
		alpha_m16(0,0) = aCRAM1; alpha_m16(1,0) = aCRAM2; alpha_m16(2,0) = aCRAM3; 
		alpha_m16(3,0) = aCRAM4; alpha_m16(4,0) = aCRAM5; alpha_m16(5,0) = aCRAM6; 
		alpha_m16(6,0) = aCRAM7; alpha_m16(7,0) = aCRAM8; 
		alpha_0 = 2.1248537104952237488e-16;
	}
	else if (method == 1){
		alpha_m16(0,0) = aH16QUAD1; alpha_m16(1,0) = aH16QUAD2; alpha_m16(2,0) = aH16QUAD3; 
		alpha_m16(3,0) = aH16QUAD4; alpha_m16(4,0) = aH16QUAD5; alpha_m16(5,0) = aH16QUAD6; 
		alpha_m16(6,0) = aH16QUAD7; alpha_m16(7,0) = aH16QUAD8; 
	}
	else if (method == 2){
		alpha_m16(0,0) = aP16QUAD1; alpha_m16(1,0) = aP16QUAD2; alpha_m16(2,0) = aP16QUAD3; 
		alpha_m16(3,0) = aP16QUAD4; alpha_m16(4,0) = aP16QUAD5; alpha_m16(5,0) = aP16QUAD6; 
		alpha_m16(6,0) = aP16QUAD7; alpha_m16(7,0) = aP16QUAD8; 
	}

	if (method == 1){
		alpha_m32(0,0) = aH32QUAD1; alpha_m32(1,0) = aH32QUAD2; alpha_m32(2,0) = aH32QUAD3; 
		alpha_m32(3,0) = aH32QUAD4; alpha_m32(4,0) = aH32QUAD5; alpha_m32(5,0) = aH32QUAD6; 
		alpha_m32(6,0) = aH32QUAD7; alpha_m32(7,0) = aH32QUAD8; alpha_m32(8,0) = aH32QUAD9; 
		alpha_m32(9,0) = aH32QUAD10; alpha_m32(10,0) = aH32QUAD11; alpha_m32(11,0) = aH32QUAD12; 
		alpha_m32(12,0) = aH32QUAD13; alpha_m32(13,0) = aH32QUAD14; alpha_m32(14,0) = aH32QUAD15; 
		alpha_m32(15,0) = aH32QUAD16; 
	}
	else if (method == 2){
		alpha_m32(0,0) = aP32QUAD1; alpha_m32(1,0) = aP32QUAD2; alpha_m32(2,0) = aP32QUAD3; 
		alpha_m32(3,0) = aP32QUAD4; alpha_m32(4,0) = aP32QUAD5; alpha_m32(5,0) = aP32QUAD6; 
		alpha_m32(6,0) = aP32QUAD7; alpha_m32(7,0) = aP32QUAD8; alpha_m32(8,0) = aP32QUAD9; 
		alpha_m32(9,0) = aP32QUAD10; alpha_m32(10,0) = aP32QUAD11; alpha_m32(11,0) = aP32QUAD12; 
		alpha_m32(12,0) = aP32QUAD13; alpha_m32(13,0) = aP32QUAD14; alpha_m32(14,0) = aP32QUAD15; 
		alpha_m32(15,0) = aP32QUAD16; 
	}
	if (order == 16){
		alpha = alpha_m16;
	}
	else{
		alpha = alpha_m32;
	}
}

//*****************************************************************************
// Matrix expotental solver
//	param A		The coefficient matrix for the system of ODEs
//	param w0		Initial condition 
//	param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
VectorD SolverType::solve(SparseMatrixD A, VectorD w0, double t){

	return solveBase(A, w0, t);
}
//*****************************************************************************
// Base Matrix expotental solver. No matrix scaling
//	param A		The coefficient matrix for the system of ODEs
//	param w0		Initial condition 
//	param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
VectorD SolverType::solveBase(SparseMatrixD A, VectorD w0, double t){

	// The sparse LU solver object
	Eigen::SparseLU<SparseMatrixCLD, COLAMDOrdering<int> > solver;
	
	// MPI stuff
	int myid = mpi.rank;
	int numprocs = mpi.size;
	int eleCount = A.rows();

	// Number of poles
	int s = theta.rows();
	SparseMatrixCLD At(A.rows(),A.cols());
	SparseMatrixCLD tempA(A.rows(),A.cols()); 
	VectorCLD w, tempB, w0cd, myW; 
	VectorD solutionVector;
	w0cd = w0.cast<std::complex<long double>>();
	SparseMatrixCLD ident(A.rows(),A.cols());

	myW = 0.*w0cd, w = 0*w0cd;
	ident.setIdentity();
	At = A.cast<std::complex<long double>>()*t;
	// analyze the sparsisty pattern
	solver.analyzePattern(tempA);

	// Loops over the imaginary poles. This is a linear solve over 8 lineary 
	// independent systems. The sum of all the independent solutions is w.
	for (int k = myid; k < s; k += numprocs){
		tempA = At - theta(k)*ident;
		tempB = alpha(k)*w0cd;
		// Compute the numerical factorization
		solver.factorize(tempA);

		myW = myW + solver.solve(tempB);
	}
	if (myid != 0){
		// Sends solution data to the master node 
		mpi.send(myW, eleCount, 0, MTAG1);
	}
	else {
		w = myW;
		// Receives data from the slave nodes
		for (int islave = 1; islave < numprocs; islave++) {
			myW = mpi.recv(myW, eleCount, islave, MTAG1);
			w = w + myW;
		}
	}
	w = 2.*w.real();
	w = w + alpha_0*w0cd;
	solutionVector = w.real().cast<double>();

	return solutionVector;
}

//*****************************************************************************
// Scaled Matrix expotental solver.
//	param A		The coefficient matrix for the system of ODEs
//	param w0		Initial condition 
//	param t		Time of the solve
//	
//	return w	Solution vector
//*****************************************************************************
VectorD SolverType::solveScale(SparseMatrixD A, VectorD w0, double t){

	// MPI stuff
	int myid = mpi.rank;
	int numprocs = mpi.size;
	int eleCount = A.rows()*A.rows();
	int beta = 10;
	double matrixReduction = pow(2.,beta);
	int matrixPower = pow(2,beta);

	// Number of poles
	int s = theta.rows();
	SparseMatrixCLD At(A.rows(),A.cols());
	SparseMatrixCLD tempA(A.rows(),A.cols()); 
	SparseMatrixCLD myExpA;
	SparseMatrixCLD expA;
	SparseMatrixCLD ident(A.rows(),A.cols());
	SparseMatrixLD expAReal;
	SparseMatrixLD expASquaredReal;
	VectorD solutionVector;

	myExpA = 0*A.cast<std::complex<long double>>();
	ident.setIdentity();
	expA = 0*A.cast<std::complex<long double>>();
	At = A.cast<std::complex<long double>>()*t/matrixReduction;

	// Loops over the imaginary poles to compute the matrix exponential
	for (int k = myid; k < s; k += numprocs){
		tempA = At - theta(k)*ident;

		// Multipy A-1*alpha(K) because the matrix is 1-1 and only left invertible 
		// because it has linearly independent columns but not rows. 
		myExpA = myExpA + MoorePenroseInv(tempA)*alpha(k);
	}
	if (myid != 0){
		// Sends solution data to the master node 
		mpi.send(myExpA, eleCount, 0, MTAG1);
	}
	else {
		expA = myExpA;
		// Receives data from the slave nodes
		for (int islave = 1; islave < numprocs; islave++) {
			myExpA = mpi.recv(myExpA, eleCount, islave, MTAG1);
			expA = expA + myExpA;
		}
	}
	expAReal = 2.*expA.real();
	expAReal = expAReal + alpha_0*ident.real();
	expASquaredReal = MatrixSquare(expAReal, matrixPower);
	solutionVector = expASquaredReal.cast<double>()*w0.cast<double>();

	return solutionVector;
}
