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
double fluxLimiter::getPsi(double r){
	double psi = (this->*fluxLimiterPtr)(r);
	return psi;
}

//*****************************************************************************
// Van Leer limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::vanLeer(double r){
	double psi = (r + std::abs(r))/(1. + r);
	return psi;
}

//*****************************************************************************
// Van Albada limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::vanAlbada(double r){
	double psi = (r + std::pow(r,2.0))/(1. + std::pow(r,2.0));
	return psi;
}

//*****************************************************************************
// Min-Mod limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::minMod(double r){
	double psi;
	if (r > 0.0){
		psi = std::min(r,1.0);
	}
	else{
		psi = 0.0;
	}
	return psi;
}

//*****************************************************************************
// SUPERBEE limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::superbee(double r){
	double psi = std::max({0.0, std::min(2.*r,1.), std::min(r, 2.)});
	return psi;
}

//*****************************************************************************
// Sweby limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::sweby(double r){
	double beta = 1.5;
	double psi = std::max({0.0, std::min(beta*r,1.), std::min(r, beta)});
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
