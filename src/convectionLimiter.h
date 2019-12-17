//*****************************************************************************
// Author: Zack Taylor
//
// flux limiter class  used for the convective species flux. Flux limiters
// give the second order upwind approximation the TVD property.
//*****************************************************************************
#ifndef FLUXLIMITER_H
#define FLUXLIMITER_H 
#include <cmath>
#include <algorithm>
#include <assert.h>

class fluxLimiter {
	// Class attributes
	public:
	// Type of flux limiter:
	//		0 = SUPERBEE
	//		1 = VanLeer
	//		2 = Van Albada
	//		3 = Min-Mod
	//		4 = Sweby
	//		5 = QUICK
	//		6 = UMIST	
	// REF: An introduction to Computational Fluid Dynamics: The finite Volume 
	//		  Method Second Edition
	int limiterType = -1;

	// Class Methods
	public:
	// Constructor
	fluxLimiter(int);
	// Function pointer that applies the limiter function
	double getPsi(double);
	private:
	// Pointer to specific limiter function
	double (fluxLimiter::*fluxLimiterPtr)(double) = nullptr;
	// Van Leer limiter
	double vanLeer(double);
	//	Van Albada limiter
	double vanAlbada(double);
	// Min-Mod limiter
	double minMod(double);
	// SUPERBEE limiter
	double superbee(double);
	// Sweby limiter
	double sweby(double);	
	// QUICK limiter
	//double quick(double);
	// UMIST limiter
	//double umist(double);	

};
#endif
