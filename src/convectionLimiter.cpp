//*****************************************************************************
// Author: Zack Taylor
//
// Reference for functions: An introduction to Computational Fluid Dynamics: 
// The finite Volume Method Second Edition
//*****************************************************************************
#include "convectionLimiter.h"
#include <iostream>

//*****************************************************************************
// Constructor
//
// @param limitType	Type of flux limiter
//    0 = SUPERBEE
//    1 = VanLeer
//    2 = Van Albada
//    3 = Min-Mod
//    4 = Sweby
//    5 = QUICK
//    6 = UMIST
//*****************************************************************************
fluxLimiter::fluxLimiter(int limitType){
	assert(limitType <= 6 and limitType >= 0);
	limiterType = limitType;

	switch (limitType) {
		case 0: fluxLimiterPtr = &fluxLimiter::superbee; break;
		case 1: fluxLimiterPtr = &fluxLimiter::vanLeer; break;
		case 2: fluxLimiterPtr = &fluxLimiter::vanAlbada; break;
		case 3: fluxLimiterPtr = &fluxLimiter::minMod; break;
		case 4: fluxLimiterPtr = &fluxLimiter::sweby; break;

	}
}


//*****************************************************************************
// Gets the Psi(r) function value
//
// @param r		Species concentration slope
//*****************************************************************************
double fluxLimiter::getPsi(const double r){
	double psi = (this->*fluxLimiterPtr)(r);
	//std::cout << r << " " << psi << std::endl;
	return psi;
}

//*****************************************************************************
// Van Leer limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::vanLeer(const double r){
	double psi = (r + std::abs(r))/(1. + r);
	return psi;
}

//*****************************************************************************
// Van Albada limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::vanAlbada(const double r){
	double psi = 1.;
	return psi;
}

//*****************************************************************************
// Min-Mod limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::minMod(const double r){
	double psi, val1, val2;
	val1 = 2./(1.+r);
	val2 = (2.*r)/(1.+r);
	psi = std::min(val1, val2);
	return psi;
}

//*****************************************************************************
// SUPERBEE limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::superbee(const double r){
	double psi, val1, val2;
	val1 = std::min(4.*r/(1.+r), 2./(1.+r));
	val2 = std::min(2.*r/(1.+r), 4./(1.+r));
	psi = std::max(val1, val2);
	return psi;
}

//*****************************************************************************
// Sweby limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::sweby(const double r){
	double beta = 1.5;
	double psi = 1.;
	return psi;
}

//*****************************************************************************
// QUICK limiter
//
// @param r		Flux slope
//*****************************************************************************
//double fluxLimiter::quick(double r){
//}

//*****************************************************************************
// UMIST limiter
//
// @param r		Flux slope
//*****************************************************************************
//double fluxLimiter::umist(double r){
//}
