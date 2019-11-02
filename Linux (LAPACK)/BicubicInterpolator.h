#pragma once
#include "Constants.h"
#include <stdint.h>

namespace PA
{
	class BicubicInterpolator //Based on implementation suggested by Numerical Recipes
	{
	public:
		BicubicInterpolator() {};
		~BicubicInterpolator() {};
	private:
		//Always the same
		const int32_t coeffMatrix[256] =
		{	 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
			 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
			-3 , 0 , 0 , 3 , 0 , 0 , 0 , 0 ,-2 , 0 , 0 ,-1 , 0 , 0 , 0 , 0 ,
			 2 , 0 , 0 ,-2 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ,
			 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
			 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ,
			 0 , 0 , 0 , 0 ,-3 , 0 , 0 , 3 , 0 , 0 , 0 , 0 ,-2 , 0 , 0 ,-1 ,
			 0 , 0 , 0 , 0 , 2 , 0 , 0 ,-2 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 1 ,
			-3 , 3 , 0 , 0 ,-2 ,-1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
			 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,-3 , 3 , 0 , 0 ,-2 ,-1 , 0 , 0 ,
			 9 ,-9 , 9 ,-9 , 6 , 3 ,-3 ,-6 , 6 ,-6 ,-3 , 3 , 4 , 2 , 1 , 2 ,
			-6 , 6 ,-6 , 6 ,-4 ,-2 , 2 , 4 ,-3 , 3 , 3 ,-3 ,-2 ,-1 ,-1 ,-2 ,
			 2 ,-2 , 0 , 0 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
			 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 ,-2 , 0 , 0 , 1 , 1 , 0 , 0 ,
			-6 , 6 ,-6 , 6 ,-3 ,-3 , 3 , 3 ,-4 , 4 , 2 ,-2 ,-2 ,-2 ,-1 ,-1 ,
			 4 ,-4 , 4 ,-4 , 2 , 2 ,-2 ,-2 , 2 ,-2 ,-2 , 2 , 1 , 1 , 1 , 1 
			};


		coeffType firstDerivative(const coeffType pIdxMinus, const coeffType pIdxPlus,
			const coeffType pDistance)
		{
			return (pIdxPlus - pIdxMinus) / pDistance;
		}

		coeffType crossDerivative(const coeffType pIdxXPlusYPlus, const coeffType pIdxXPlusYMinus,
			const coeffType pIdxXMinusYPlus, const coeffType pIdxXMinusYMinus,
			const coeffType pXRDistance, const coeffType pYDistance)
		{
			return (pIdxXPlusYPlus - pIdxXPlusYMinus - pIdxXMinusYPlus + pIdxXMinusYMinus) /
				(pXRDistance*pYDistance);
		}

		//Input values are 16 functions values -> gradients and crossderivatives are calculated
		//We input a 4 x 4 grid, and will interpolate the central square.
		//(x0,y0) (x1,y0) (x2,y0) (x3, y0),
		//(x0,y1) (x1,y1) (x2,y1) (x3, y1),
		//(x0,y2) (x1,y2) (x2,y2) (x3, y2),
		//(x0,y3) (x1,y3) (x2,y3) (x3, y3)
		void setupInputFromPoints(coeffType(&inputVector)[16], coeffType(&pFunctionValues)[16],
			const coeffType pXDistance, const coeffType pYDistance)
		{
			//First, assign coordinates for target square starting in bottom left in counterclockwise order.
			inputVector[0] = pFunctionValues[5]; //(x1,y1)
			inputVector[1] = pFunctionValues[6]; //(x2,y1)
			inputVector[2] = pFunctionValues[10]; //(x2,y2)
			inputVector[3] = pFunctionValues[9]; //(x1,y2)

			//get first derivatives in first coordinates direction
			inputVector[4] = firstDerivative(pFunctionValues[4], pFunctionValues[6], 2*pXDistance);
			inputVector[5] = firstDerivative(pFunctionValues[5], pFunctionValues[7], 2*pXDistance);
			inputVector[6] = firstDerivative(pFunctionValues[9], pFunctionValues[11], 2*pXDistance);
			inputVector[7] = firstDerivative(pFunctionValues[8], pFunctionValues[10], 2*pXDistance);

			//get first deriviatives in the second coordinates direction
			inputVector[8] = firstDerivative(pFunctionValues[1], pFunctionValues[9], 2*pYDistance);
			inputVector[9] = firstDerivative(pFunctionValues[2], pFunctionValues[10], 2*pYDistance);
			inputVector[10] = firstDerivative(pFunctionValues[6], pFunctionValues[14], 2*pYDistance);
			inputVector[11] = firstDerivative(pFunctionValues[5], pFunctionValues[13], 2*pYDistance);

			//get crossderivatives
			inputVector[12] = crossDerivative(pFunctionValues[10], pFunctionValues[2],
				pFunctionValues[8], pFunctionValues[0], 2*pXDistance, 2*pYDistance);

			inputVector[13] = crossDerivative(pFunctionValues[11], pFunctionValues[3],
				pFunctionValues[9], pFunctionValues[1], 2*pXDistance, 2*pYDistance);

			inputVector[14] = crossDerivative(pFunctionValues[15], pFunctionValues[7],
				pFunctionValues[13], pFunctionValues[5], 2*pXDistance, 2*pYDistance);

			inputVector[15] = crossDerivative(pFunctionValues[14], pFunctionValues[6],
				pFunctionValues[12], pFunctionValues[4], 2*pXDistance, 2*pYDistance);
		}

	public:

		//Input values are 4 function values, 4 radial gradients, 4 theta gradients, 4 crossderivatives
		//These are always for the respective corner in the order 0,1,2,3 (counterclockwise),
		//starting in the lower left of a rectangular grid cell.
		void calcCoeffs(coeffType(&pInputValues)[16], coeffType(&pOutputCoeffs)[16])
		{
			for (size_t i = 0; i < 16; i++)
			{
				coeffType coeff = 0;
				for (size_t k = 0; k < 16; k++)
				{
					coeff += coeffMatrix[k + 16 * i] * pInputValues[k];
				}
				pOutputCoeffs[i] = coeff;
			}
		}

		//Input values are 16 functions values -> gradients and crossderivatives are calculated
		//We input a 4 x 4 grid. See method setupInputFromPoints above.
		void calcDerivativesAndCoeffs(coeffType(&pFunctionValues)[16], coeffType(&pOutputCoeffs)[16],
			const coeffType pXDistance, const coeffType pYDistance)
		{
			coeffType inputVector[16];
			setupInputFromPoints(inputVector, pFunctionValues, pXDistance, pYDistance);
			calcCoeffs(inputVector, pOutputCoeffs); //now calc bicubic interpolation coefficients and save them in pOutputCoeffs.
		}


		coeffType interpolateBicubic(const coeffType pXCoord, const coeffType pYCoord,
			coeffType(&pInputCoeffs)[16])
		{
			coeffType result = 0;
			for (int i = 3; i >= 0; i--)
			{
				result = pXCoord * result + ((pInputCoeffs[i * 4 + 3] * pYCoord +
					pInputCoeffs[i * 4 + 2])*pYCoord + pInputCoeffs[i * 4 + 1])*pYCoord +
					pInputCoeffs[i * 4];
				//could also calculate (interpolated) first derivatives here, if desired
			}
			return result;
		}
	};
}

