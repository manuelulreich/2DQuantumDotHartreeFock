#pragma once
#include "Constants.h"
#include <cstring>
#include <array>
#include <algorithm>
#include <iostream>
//#include <mkl.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <functional>
#include <complex>

namespace PA
{

	struct Hamiltonian
	{
	public:
		Hamiltonian(long pGridSizeR, long pGridSizeTheta)
			: gridSizeR(pGridSizeR), gridSizeTheta(pGridSizeTheta), matrixDim(pGridSizeR*pGridSizeTheta),
			matrixSize(pGridSizeR*pGridSizeR*pGridSizeTheta*pGridSizeTheta), dThetaRad(2 * M_PI / gridSizeTheta),
			invThetaDeltaSquared(gridSizeTheta*gridSizeTheta/(4 * M_PI * M_PI))
		{
			mData = nullptr;
		}

		~Hamiltonian()
		{
			if (mData != nullptr)
			{
				delete[] mData;
				mData = nullptr;
			}
		}

		coeffType *mData;
		//constants but dependent on startup
		const long gridSizeR;
		const long gridSizeTheta;
		const long matrixSize;
		const long matrixDim;
		const constType dThetaRad;
		const constType invThetaDeltaSquared;

		///Theta-major row setup
		void AsymmetricSetup()
		{
			if (mData != nullptr) delete[] mData;
			mData = new coeffType[matrixSize]();
			memset(mData, 0, matrixSize * sizeof(coeffType));

			//Executed with ifs for code-brevity/-legibility. Could be optimized by creating code without ifs (becomes much longer).
			for (long Ri = 0; Ri < gridSizeR; Ri++)
			{
				coeffType current_r = radial_start + Ri * dR;
				coeffType rSquared = current_r * current_r;
				coeffType temp = -invThetaDeltaSquared * coeffType(0.5) / rSquared;
				//long radialOffset = gridSizeR * Ri;
				for (long Tj = 0; Tj < gridSizeTheta; Tj++)
				{
					long currentIdx = gridSizeTheta * Ri + Tj; //offset in our row
					long rowOffset = matrixDim*currentIdx; //offset from start of matrix to start of current row

					//Main diagonal: 1/DeltaR^2 + 1/(DeltaTheta^2*R^2) + omega^2*R^2/2
					mData[rowOffset + currentIdx] = invRDeltaSquared + invThetaDeltaSquared / rSquared + omegaSquaredHalf * rSquared;

					//Right and left offdiagonal: -1/(2*DeltaTheta^2*R^2)
					if (Tj == 0)
					{
						mData[rowOffset + currentIdx + 1] = temp;
						mData[rowOffset + currentIdx - 1 + gridSizeTheta] = temp;
					}
					else if (Tj == gridSizeTheta - 1)
					{
						mData[rowOffset + currentIdx + 1 - gridSizeTheta] = temp;
						mData[rowOffset + currentIdx - 1] = temp;
					}
					else
					{
						mData[rowOffset + currentIdx + 1] = temp;
						mData[rowOffset + currentIdx - 1] = temp;
					}

					//Far right offdiagonal: -1/(2*DeltaR^2) -1/(DeltaR*R)
					if (gridSizeTheta + currentIdx < matrixDim) //else index out of bounds
					{
						mData[rowOffset + gridSizeTheta + currentIdx] = -invRDeltaSquaredHalf - coeffType(0.25)*invRDelta / current_r;
					}

					//Far left offdiagonal -1/(2*DeltaR^2) +1/(DeltaR*R)
					if (-gridSizeTheta + currentIdx >= 0)
					{
						mData[rowOffset - gridSizeTheta + currentIdx] = -invRDeltaSquaredHalf + coeffType(0.25)*invRDelta / current_r;
					}
				}
			}

			//if (ApplyTransitionCondition) //we don't need this because with our current setup the contribution from oppositeTheta is 0!
		}

		void freeResources()
		{
			if (mData != nullptr)
			{
				delete[] mData;
				mData = nullptr;
			}
		}

		///We transform the problem to u(r,theta) = psi(r,theta)*r^(1/2) and multiply the entire problem by 2*r^(1/2)
		void SymmetricSetup(bool useBandedStorage = false)
		{
			if (mData != nullptr) delete[] mData;
			mData = new coeffType[matrixSize]();
			memset(mData, 0, matrixSize * sizeof(coeffType));

			for (long Ri = 0; Ri < gridSizeR; Ri++)
			{
				//long radialOffset = gridSizeR * Ri;
				coeffType current_r = radial_start + Ri * dR;
				coeffType rSquared = current_r * current_r;
				coeffType temp = -invThetaDeltaSquared / rSquared*coeffType(0.5);
				for (long Tj = 0; Tj < gridSizeTheta; Tj++)
				{
					long currentIdx = gridSizeTheta * Ri + Tj; //offset in our row
					long rowOffset = matrixDim*currentIdx; //offset from start of matrix to start of current row

					//Main diagonal: 1/DeltaR^2 + 1/(DeltaTheta^2*R^2) + omega^2*R^2/2
					mData[rowOffset + currentIdx] = invRDeltaSquared + invThetaDeltaSquared / rSquared
						+ omega*omega * rSquared * coeffType(0.5) - coeffType(0.125)/ rSquared;

					//Right and left offdiagonal: -1/(2*DeltaTheta^2*R^2)
					if (Tj == 0)
					{
						mData[rowOffset + currentIdx + 1] = temp;
						mData[rowOffset + currentIdx - 1 + gridSizeTheta] = temp;
					}
					else if (Tj == gridSizeTheta - 1)
					{
						mData[rowOffset + currentIdx + 1 - gridSizeTheta] = temp;
						mData[rowOffset + currentIdx - 1] = temp;
					}
					else
					{
						mData[rowOffset + currentIdx + 1] = temp;
						mData[rowOffset + currentIdx - 1] = temp;
					}

					//Far right offdiagonal: -1/(2*DeltaR^2) -1/(DeltaR*R)
					if (gridSizeTheta + currentIdx < matrixDim) //else index out of bounds
					{
						//current_r = radial_start + (Ri + 1)*dR;
						mData[rowOffset + gridSizeTheta + currentIdx] = -invRDeltaSquared * coeffType(0.5);
					}

					//Far left offdiagonal -1/(2*DeltaR^2) +1/(DeltaR*R)
					if (-gridSizeTheta + currentIdx >= 0)
					{
						//current_r = radial_start + (Ri + 1)*dR;
						mData[rowOffset - gridSizeTheta + currentIdx] = -invRDeltaSquared * coeffType(0.5);
					}
				}
			}


			for (long Tj = 0; Tj < gridSizeTheta; Tj++) //attempting to find approximation for ghost point at r = -dR/2
														//(does not make necessarily sense analytically, since we are in polar coords)
														//necessary for 2nd derivative of point r= +dR/2
			{
				long rowOffset = matrixDim * Tj; //offset from start of matrix to start of current row

				//Main diagonal: 1/DeltaR^2 + 1/(DeltaTheta^2*R^2) + omega^2*R^2/2
				mData[rowOffset + Tj] -= -invRDeltaSquared * coeffType(0.5); //dual negative, once to switch the sign of current psi value, once for sign of coefficient
			}

			//Print();

			if (useBandedStorage) //basically deprecated/superfluous, but kept if desirable in future
			{
				//PA::coeffType * buffer = new PA:coeffType[(2 * gridSizeTheta + 1)]
				const long maxRow = 2 * gridSizeTheta + 1;
				for (long col = 0; col < gridSizeTheta; col++) //first gridSizeTheta rows
				{
					for (long row = gridSizeTheta+col; row >=0; row--)
					{
						long sourceRowOffset = matrixDim *row;
						long moveRows = gridSizeTheta - col;
						long targetRowOffset = matrixDim *(row+moveRows);
						mData[targetRowOffset + col] = mData[sourceRowOffset + col];
						mData[sourceRowOffset + col] = 0;
					}
				}

				for (long col = gridSizeTheta + 1; col < matrixDim; col++) //rest of the columns
				{
					for (long row = 0; row < 2*gridSizeTheta+1; row++)
					{
						long moveRows = col - gridSizeTheta;
						long targetRowOffset = matrixDim * row;
						long sourceRowOffset = matrixDim * (row + moveRows);
						if (sourceRowOffset >= matrixSize) continue;
						mData[targetRowOffset + col] = mData[sourceRowOffset + col];
						mData[sourceRowOffset + col] = 0;
					}
				}
			}

			//Print();
		}

		enum CoeffApplicationType {AsymNegCmplxTimeStep = 0, SymNegCmplxTimeStep = 1, AsymHamil = 2};

		inline void SetupCoeffs(coeffType &cRPlusOne, const coeffType &current_r,
			coeffType &cRMinusOne, coeffType &cRThetaSame, const coeffType &rSquared, coeffType &cThetaPMOne, CoeffApplicationType pApplicationType)
		{
			switch (pApplicationType)
			{
			case PA::Hamiltonian::AsymNegCmplxTimeStep:
				//Asymmetric setup - Negative Complex Timestep
				cRPlusOne = -(-invRDeltaSquaredHalf - coeffType(0.25) * invRDelta / current_r) * dT;
				cRMinusOne = -(-invRDeltaSquaredHalf + coeffType(0.25) * invRDelta / current_r) * dT;
				cRThetaSame = -(invRDeltaSquared + invThetaDeltaSquared / rSquared + omegaSquaredHalf * rSquared) * dT + coeffType(1.);
				cThetaPMOne = -(-invThetaDeltaSquared * coeffType(0.5) / rSquared) * dT;
				break;
			case PA::Hamiltonian::SymNegCmplxTimeStep:
				//Symmetric setup - Negative Complex Timestep
				cRPlusOne = (-invRDeltaSquared)*dT;
				cRMinusOne = (-invRDeltaSquared)*dT;
				cRThetaSame = (invRDeltaSquared*coeffType(2.) + invThetaDeltaSquared * coeffType(2.) / rSquared
					+ omega * omega * rSquared - coeffType(0.25) / rSquared)*dT + coeffType(1.);
				cThetaPMOne = (-invThetaDeltaSquared / rSquared)*dT;
				break;
			case PA::Hamiltonian::AsymHamil:
				//Asymmetric setup - Apply Hamiltonian to Psi
				cRPlusOne = (-invRDeltaSquaredHalf - coeffType(0.25)*invRDelta / current_r);
				cRMinusOne = (-invRDeltaSquaredHalf + coeffType(0.25)*invRDelta / current_r);
				cRThetaSame = (invRDeltaSquared + invThetaDeltaSquared / rSquared + omegaSquaredHalf * rSquared);
				cThetaPMOne = (-invThetaDeltaSquared * coeffType(0.5) / rSquared);
				break;
			default:
				break;
			}
		}

		void ApplyCoeffsToPsi(coeffType* PsiOld, coeffType* PsiNew, CoeffApplicationType pApplicationType)
		{
			//ri == 0
			{
				long ri = 0;
				coeffType current_r = radial_start + ri * dR, rSquared = current_r * current_r;
				coeffType cRPlusOne = 0, cRMinusOne = 0, cRThetaSame = 0, cThetaPMOne = 0;
				SetupCoeffs(cRPlusOne, current_r, cRMinusOne, cRThetaSame, rSquared, cThetaPMOne, pApplicationType);

				long tj = 0;
				long oppositeThetaIdx = (gridSizeTheta / 2 + tj) % gridSizeTheta;
				long currentIdx = gridSizeTheta*ri+tj; //target index for new psi
				//TODO: PROPERLY EMPLOY CURRENTIDX VARIABLE, AND POSSIBLY EXTRACT FUNCTION FOR GENERATING PSINEW!
				PsiNew[currentIdx] = PsiOld[currentIdx + gridSizeTheta] * cRPlusOne + /*PsiOld[oppositeThetaIdx] * cRMinusOne +*/
					PsiOld[currentIdx] * cRThetaSame + (PsiOld[currentIdx + 1] + PsiOld[currentIdx + gridSizeTheta - 1])*cThetaPMOne;
				for (tj = 1; tj < gridSizeTheta-1; tj++)
				{
					currentIdx = gridSizeTheta * ri + tj; //target index for new psi
					oppositeThetaIdx = (gridSizeTheta / 2 + tj) % gridSizeTheta;
					PsiNew[currentIdx] = PsiOld[currentIdx + gridSizeTheta] * cRPlusOne + /*PsiOld[oppositeThetaIdx] * cRMinusOne +*/
						PsiOld[currentIdx] * cRThetaSame + (PsiOld[currentIdx + 1] + PsiOld[currentIdx - 1])*cThetaPMOne;
				}
				tj = gridSizeTheta - 1;
				currentIdx = gridSizeTheta * ri + tj; //target index for new psi
				oppositeThetaIdx = (gridSizeTheta / 2 + tj) % gridSizeTheta;
				PsiNew[currentIdx] = PsiOld[currentIdx + gridSizeTheta] * cRPlusOne + /*PsiOld[oppositeThetaIdx] * cRMinusOne +*/
					PsiOld[currentIdx] * cRThetaSame + (PsiOld[currentIdx - gridSizeTheta + 1] + PsiOld[currentIdx - 1])*cThetaPMOne;
			}

			for (long ri = 1; ri < gridSizeR-1; ri++) //not for ri == 0 and ri == gridSizeR-1
			{
				coeffType current_r = radial_start + ri * dR, rSquared = current_r * current_r;
				coeffType cRPlusOne = 0, cRMinusOne = 0, cRThetaSame = 0, cThetaPMOne = 0;
				SetupCoeffs(cRPlusOne, current_r, cRMinusOne, cRThetaSame, rSquared, cThetaPMOne, pApplicationType);

				long tj = 0;
				long currentIdx = gridSizeTheta * ri + tj; //target index for new psi
				PsiNew[currentIdx] = PsiOld[currentIdx + gridSizeTheta] * cRPlusOne + PsiOld[currentIdx - gridSizeTheta] * cRMinusOne +
					PsiOld[currentIdx] * cRThetaSame + (PsiOld[currentIdx + 1] + PsiOld[currentIdx + gridSizeTheta - 1])*cThetaPMOne;
				for (tj = 1; tj < gridSizeTheta-1; tj++)
				{
					currentIdx = gridSizeTheta * ri + tj; //target index for new psi
					PsiNew[currentIdx] = PsiOld[currentIdx + gridSizeTheta] * cRPlusOne + PsiOld[currentIdx - gridSizeTheta] * cRMinusOne +
						PsiOld[currentIdx] * cRThetaSame + (PsiOld[currentIdx + 1] + PsiOld[currentIdx - 1])*cThetaPMOne;
				}
				tj = gridSizeTheta - 1;
				currentIdx = gridSizeTheta * ri + tj; //target index for new psi
				PsiNew[currentIdx] = PsiOld[currentIdx + gridSizeTheta] * cRPlusOne + PsiOld[currentIdx - gridSizeTheta] * cRMinusOne +
					PsiOld[currentIdx] * cRThetaSame + (PsiOld[currentIdx - gridSizeTheta  + 1] + PsiOld[currentIdx - 1])*cThetaPMOne;
			}

			//ri == gridSizeR-1, but tj < gridSizeTheta-1
			{
				long ri = gridSizeR - 1;
				coeffType current_r = radial_start + ri * dR, rSquared = current_r * current_r;
				coeffType cRPlusOne = 0, cRMinusOne = 0, cRThetaSame = 0, cThetaPMOne = 0;
				SetupCoeffs(cRPlusOne, current_r, cRMinusOne, cRThetaSame, rSquared, cThetaPMOne, pApplicationType);

				long tj = 0;
				long currentIdx = gridSizeTheta * ri + tj; //target index for new psi
				PsiNew[currentIdx] = PsiOld[currentIdx - gridSizeTheta] * cRMinusOne +
					PsiOld[currentIdx] * cRThetaSame + (PsiOld[currentIdx + 1] + PsiOld[currentIdx + gridSizeTheta - 1])*cThetaPMOne;
				for (tj = 1; tj < gridSizeTheta-1; tj++)
				{
					currentIdx = gridSizeTheta * ri + tj; //target index for new psi
					PsiNew[currentIdx] = PsiOld[currentIdx - gridSizeTheta] * cRMinusOne +
						PsiOld[currentIdx] * cRThetaSame + (PsiOld[currentIdx + 1] + PsiOld[currentIdx - 1])*cThetaPMOne;
				}
				tj = gridSizeTheta - 1;
				currentIdx = gridSizeTheta * ri + tj; //target index for new psi
				PsiNew[currentIdx] = PsiOld[currentIdx - gridSizeTheta] * cRMinusOne +
					PsiOld[currentIdx] * cRThetaSame + (PsiOld[currentIdx - 1] + PsiOld[currentIdx - gridSizeTheta + 1])*cThetaPMOne;
			}
		}

		bool Print(std::ostream& pTarget = std::cout, std::function<coeffType(coeffType)> pFunc = [](coeffType val) { return val; })
		{
			if (matrixDim > 40 && &pTarget == &std::cout)
			{
				std::cout << "Matrix too large to print!\n";
				for (size_t i = 0; i < matrixDim; i++)
				{
					pTarget << mData[matrixDim*i + i] << " ";
				}
				return false;
			}
			else
			{
				pTarget << std::fixed << std::setprecision(9);
				std::cout << std::fixed << std::setprecision(3) << "Printing Hamiltonian.\n";
				for (size_t row = 0; row < matrixDim; row++)
				{
					for (size_t i = 0; i < matrixDim; i++)
					{
						pTarget << pFunc(mData[matrixDim*row + i]) << " ";
					}
					pTarget << "\n";
				}
				return true;
			}
		}

		coeffType Reference1DHarmonicCart(long n, coeffType x)
		{
			return exp(-x*x*coeffType(0.5))*std::hermite(n, x);
		}

		coeffType Reference2DHarmonicCart(long nx, long ny, coeffType x, coeffType y)
		{
			return Reference1DHarmonicCart(nx, x)*Reference1DHarmonicCart(ny, y);
		}

		void Print2DHarmonicCart(std::ofstream& pFileStream, long nx, long ny)
		{
			//normalizing
			coeffType norm = 0;
			for (long Ri = 0; Ri < gridSizeR; Ri++)
			{
				for (long Tj = 0; Tj < gridSizeTheta; Tj++)
				{
					coeffType radius = Ri * dR + radial_start, theta = dThetaRad * Tj;
					coeffType x = radius * cos(theta), y = radius * sin(theta);
					auto val = Reference2DHarmonicCart(nx, ny, x, y);
					norm += (2 * Ri + 1)*val * val;
				}
			}

			norm = 1 / norm;

			for (long Ri = 0; Ri < gridSizeR; Ri++)
			{
				coeffType radius = Ri * dR + radial_start;
				for (long Tj = 0; Tj < gridSizeTheta; Tj++)
				{
					coeffType theta = dThetaRad * Tj;
					coeffType x = radius * cos(theta), y = radius * sin(theta);
					auto val = Reference2DHarmonicCart(nx, ny, x, y);
					val *= val*norm;
					pFileStream << radius << " " << theta << " " << val << "\n";
				}
				{
					coeffType theta = 0; //repeat first point for gnuplot
					coeffType x = radius * cos(theta), y = radius * sin(theta);
					auto val = Reference2DHarmonicCart(nx, ny, x, y);
					val *= val * norm;
					pFileStream << radius << " " << theta << " " << val << "\n";
				}
				pFileStream << "\n";
			}
			
		}

		std::complex<coeffType> Reference2DHarmonicPolar(long k, long m, coeffType r, coeffType theta)
		{
			coeffType radial = std::exp(-r * r*0.5)*std::pow(r, abs(m*0.5))*std::assoc_laguerre(k, abs(m), r*r);
			std::complex<coeffType> angular = std::exp(std::complex<coeffType>(0, 1)*coeffType(m)*theta);
			return radial * angular;
		}

		void Print2DHarmonicPolar(std::ofstream& pFileStream, long k, long m)
		{
			//normalizing

			pFileStream << std::scientific << std::setprecision(6);
			coeffType norm = 0;
			for (long Ri = 0; Ri < gridSizeR; Ri++)
			{
				for (long Tj = 0; Tj < gridSizeTheta; Tj++)
				{
					coeffType radius = Ri * dR + radial_start, theta = dThetaRad * Tj;
					//coeffType x = radius * cos(theta), y = radius * sin(theta);
					auto val = Reference2DHarmonicPolar(k, m, radius, theta);
					norm += (2 * Ri + 1)*std::abs(val);
				}
			}

			norm = 1 / norm;

			for (long Ri = 0; Ri < gridSizeR; Ri++)
			{
				coeffType radius = Ri * dR + radial_start;
				for (long Tj = 0; Tj < gridSizeTheta; Tj++)
				{
					coeffType theta = dThetaRad * Tj;
					//coeffType x = radius * cos(theta), y = radius * sin(theta);
					auto val = Reference2DHarmonicPolar(k, m, radius, theta);
					float out = std::abs(val) * norm;
					pFileStream << radius << " " << theta << " " << val.real() << "\n";
				}
				{
					coeffType theta = 0; //repeat first point for gnuplot
					//coeffType x = radius * cos(theta), y = radius * sin(theta);
					auto val = Reference2DHarmonicPolar(k, m, radius, theta);
					float out = std::abs(val) * norm;
					pFileStream << radius << " " << theta << " " << val.real() << "\n";
				}
				pFileStream << "\n";
			}

		}
	};
}
