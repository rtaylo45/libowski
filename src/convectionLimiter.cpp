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
//    5 = First order upwind
//    6 = QUICK
//    7 = UMIST
//*****************************************************************************
fluxLimiter::fluxLimiter(int limitType){
	setLimiterFunction(limitType);
}

//*****************************************************************************
// Changes the flux limiter function
//
// @param limitType	Type of flux limiter
//    0 = SUPERBEE
//    1 = VanLeer
//    2 = Van Albada
//    3 = Min-Mod
//    4 = Sweby
//    5 = First order upwind
//    6 = QUICK
//    7 = UMIST
//*****************************************************************************
void fluxLimiter::setLimiterFunction(int limitType){
	assert(limitType <= 6 and limitType >= 0);
	limiterType = limitType;

	switch (limitType) {
		case 0: fluxLimiterPtr = &fluxLimiter::superbee; break;
		case 1: fluxLimiterPtr = &fluxLimiter::vanLeer; break;
		case 2: fluxLimiterPtr = &fluxLimiter::vanAlbada; break;
		case 3: fluxLimiterPtr = &fluxLimiter::minMod; break;
		case 4: fluxLimiterPtr = &fluxLimiter::sweby; break;
		case 5: fluxLimiterPtr = &fluxLimiter::firstOrder; break;
		case 6: fluxLimiterPtr = &fluxLimiter::muscl; break;

	}
}

//*****************************************************************************
// Gets the Psi(r) function value
//
// @param r		Species concentration slope
//*****************************************************************************
double fluxLimiter::getPsi(const double r){
	double psi = (this->*fluxLimiterPtr)(r);
	return psi;
}

//*****************************************************************************
// First order upwind difference
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::firstOrder(const double r){
	double psi = 0.0;
	return psi;
}

//*****************************************************************************
// Van Leer limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::vanLeer(const double r){
	double psi = (r + std::abs(r))/(1. + r);
	if (isnan(psi)){psi = 0.0;};
	return psi;
}

//*****************************************************************************
// Van Albada limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::vanAlbada(const double r){
	double psi = (r + std::pow(r, 2.))/(1. + std::pow(r,2.));
	if (isnan(psi)){psi = 0.0;};
	return psi;
}

//*****************************************************************************
// Min-Mod limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::minMod(const double r){
	double psi, val1, val2;
	if (r > 0){
		psi = std::min(r,1.);
	}
	else{
		psi = 0.0;
	}
	if (isnan(psi)){psi = 0.0;};
	return psi;
}

//*****************************************************************************
// SUPERBEE limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::superbee(const double r){
	double psi, val1, val2;
	val1 = std::min(2.*r, 1.);
	val2 = std::min(r, 2.);
	psi = std::max(val1, val2);
	psi = std::max(0., psi);
	if (isnan(psi)){psi = 0.0;};
	return psi;
}

//*****************************************************************************
// Sweby limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::sweby(const double r){
	double beta = 1.5, psi, val1, val2;
	val1 = std::min(beta*r, 1.);
	val2 = std::min(r, beta);
	psi = std::max(val1, val2);
	psi = std::max(0., psi);
	if (isnan(psi)){psi = 0.0;};
	return psi;
}
//*****************************************************************************
// MUSCL limiter
//
// @param r		Flux slope
//*****************************************************************************
double fluxLimiter::muscl(const double r){
	double psi, val1, val2;
	val1 = std::min(2.*r, (r+1.)/2.);
	val2 = std::min(val1, 2.);
	psi = std::max(0., val2);
	if (isnan(psi)){psi = 0.0;};
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
