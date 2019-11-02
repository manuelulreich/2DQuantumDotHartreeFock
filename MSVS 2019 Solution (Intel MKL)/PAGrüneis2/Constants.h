#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

namespace PA
{
#define DOUBLEPR
#ifdef DOUBLEPR
	using constType = double;
	using coeffType = double;
#else
	using constType = float;
	using coeffType = float;
#endif // DOUBLEPR


	//simulation parameters (compile-time cost for now)
	constexpr constType dR = 0.0625; //grid spacing in the radial direction
	constexpr constType omega = 1.; //harmonic oscillator frequency
	constexpr constType dT = 1.e-4; //has to have twice the order of of magnitude of dR, ie if dR = 1e-2 then dT = 1e-4

	//preprocessed constants
	constexpr constType radial_start = dR / 2.; //innermost ring has this radius
	constexpr constType dArea = M_PI * dR*dR;
	constexpr constType invRDeltaSquared = 1. / (dR * dR);
	//constexpr constType invThetaDeltaSquared = 1. / (dThetaRad * dThetaRad);
	constexpr constType invRDeltaSquaredHalf = 0.5 / (dR * dR);
	constexpr constType invRDelta = 1. / dR;
	constexpr constType omegaSquaredHalf = omega * omega*0.5;

}