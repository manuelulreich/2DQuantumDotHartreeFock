#pragma once
#include "Hamiltonian.h"
#include "BicubicInterpolator.h"
#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <future>

namespace PA
{
	using valIdx = std::pair<coeffType, size_t>;

	class Solver
	{

	public:
		Solver(bool pSymmetric, long pGridSizeR, long pGridSizeTheta, bool pBinaryOutput);
		~Solver();

		//Initializes Psi with gauss test peak.
		void InitPsi();

		//Executes symmetric or asymmetric setup of hamiltonian.
		void InitHamiltonian();

		//Initializes work ptr with N*matrixDim values corresponding to N eigenvectors;
		void initWork(const size_t pEigenvectorCount);

		//Applies generalized eigenvalue + eigenvector algorithm to hamiltonian.
		void Solve(size_t pEigenVecCount, bool pAlternative = false);

		//Loads eigenvectors from file to memory using work ptr.
		void LoadEigenVecs(long pCount, bool fromBinary, std::string pFileName = std::string());

		//Loads eigenvalues from file to memory using real ptr.
		void loadEigenValues(size_t pCount, std::string pFileName = std::string());

		//Converts .bdata to .txt file
		void convertData(long pCount, bool fromBinary, std::string pFileName = std::string());

	private:
		//Calculates inner product between two row-major vectors.
		PA::coeffType innerProductRowMajor(PA::coeffType * work, size_t vec1, size_t vec2);

		//Calculates inner product between two column-major vectors.
		PA::coeffType innerProductColMajor(PA::coeffType* work, size_t vec1, size_t vec2);

		//Projects vector 1 onto vector 2, saving result in vector 2. (For Gram-Schmidt).
		void ProjectRowMajor(PA::coeffType * work, size_t vec1, size_t vec2);

		//Projects vector 1 onto vector 2, saving result in vector 2. (For Gram-Schmidt).
		void ProjectColMajor(PA::coeffType* work, size_t vec1, size_t vec2);

		//Employs modified Gram-Schmidt algorithm to orthonormalize all desired eigenvectors.
		void Orthonormalize(PA::coeffType * work, std::vector<valIdx> & desiredVecs, const size_t & pEigenVecCount);

		//Only normalizes eigenvectors.
		void NormalizeEigenVecs(const size_t &pEigenVecCount, PA::coeffType * work);

		//Searches resulting eigenvalues for the lowest N eigenvalues, and remembers their index.
		void FindLowestEigenVecs(const size_t &pEigenVecCount, std::vector<valIdx> &desiredEigenVecs, PA::coeffType * real);

		//Prints desired eigenvalues to console.
		void PrintEV(std::vector<PA::valIdx> &desiredEigenVecs);

		//Applies a set of lapack routines for eigenvecs + eigenvalues.
		void GeneralSolve(const size_t pEigenVecCount);

		//Returns stem of the filename.
		std::string getFilenameStem() {
			return std::to_string(gridSizeR) + "x" + std::to_string(gridSizeTheta) + "_dr_" + std::to_string(dR);
		}

		//Returns path to file in which eigenvectors should be saved.
		std::filesystem::path getFilePath(const size_t pEigenVecCount, const bool isBinary = false);

		//Sub-method for PrintEigenVecsAsText. Eigenvectors are stored in columns.
		void PrintEVecRow_EVecsColMajor(const size_t &i, const size_t &j, const size_t & pEigenVecCount, std::ofstream & myfile,
			std::vector<PA::valIdx> * targetColumns, PA::coeffType * work, const size_t &count);

		//Sub-method for PrintEigenVecsAsText. Eigenvectors are stored in rows.
		void PrintEVecRow_EVecsRowMajor(const size_t& i, const size_t& j, const size_t& pEigenVecCount, std::ofstream& myfile,
			std::vector<PA::valIdx>* targetRows, PA::coeffType* work, const size_t& count);

		//Prints desired eigenvectors to console, or optionally to specified ofstream.
		void PrintEigenVecsAsText(const size_t &pEigenVecCount, PA::coeffType * work, std::vector<valIdx> *targetColumns,
			bool pEVecsColMajor, std::ofstream * myfile = nullptr);

		//Prints desired eigenvectors in binary format to file.
		void printEigenVecsBinary(const size_t &pEigenVecCount, PA::coeffType * work, std::vector<valIdx> *targetColumns,
			std::ofstream * myfile = nullptr);

		//Writes eigenvector (representad as column in row-major matrix) to psiOld.
		void ExtractPsiMinor(long EVTarget);

		//Writes eigenvector (representad as row in row-major matrix) to psiOld.
		void ExtractPsiMajor(long EVTarget);

		//Inner product specifically between psiOld and psiNew.
		PA::coeffType PsiOldPsiNew();

		//Applies general lapack eigenvector + eigenvalue algorithm.
		void DriverRoutineSolve(const size_t pEigenVecCount);

		//Uses driver routine from MKL to solve symmetric matrix for eigenvalues and -vectors.
		void SymmetricDriverRoutineSolve(const size_t pEigenVecCount);

		void ApplyFofR(const size_t& pEigenVecCount);

		//Applies general lapack eigenvector and eigenvalue algorithms.
		void SymmetricSolve(const size_t pEigenVecCount);

		//Outputs the orthogonality relations between <Psi_N|H|Psi_M> for all N and M.
		void PrintOrthoRelations(const size_t &pEigenVecCount, std::vector<PA::valIdx> &desiredEigenVecs, const bool isPsiMajor);

		//Calculates [NM|PQ] for use with Slater-Condon rules -> two body propagators.
		coeffType integralNMPQ(coeffType* psiN, coeffType* psiM, coeffType* psiP, coeffType* psiQ, 
			size_t pGridSizeR, size_t pGridSizeTheta);

		//Calculates coulomb energy correction for matrix element with same wavefunctions.
		coeffType slaterCondonPsiCoulombPsi(std::vector<PA::valIdx> &desiredEigenVecs);

	public:

		//Returns filename for eigenvalue file
		std::string getEValFilename() {
#ifdef DOUBLEPR
			std::string prefix = "d";
#else
			std::string prefix = "s";
#endif // DOUBLEP

			std::string symmetry = std::string();
			if (symmetric) symmetry = "sym";
			else symmetry = "asym";

			return prefix + "_" + symmetry + "_" + getFilenameStem() + "_eigenvals.txt";
		}

		//Returns filename for eigenvector file
		std::string getEVecFilename(const size_t pEigenVecCount, const bool pIsBinary) {
			return getFilePath(pEigenVecCount, pIsBinary).filename().string();
		}


		void extendPsi(coeffType * &pSourcePsi, coeffType * &pTargetPsi);

		//interpolates single psi 
		void interpolatePsi(const size_t &intervalsTheta, const size_t &intervalsR,
			coeffType * &interpolatedPsi, coeffType * &startPsi);

		//Interpolates single psi for test purposes
		void interpolatePsiTest();

		//Specialized negative imaginary time propagation to find groundstate.
		void SpecializedSolve(bool pShouldPrint);

		//Only calculates eigenvalues with lapack routines.
		void CalcEigenValues(bool pSaveData);

		//Prints psiOld to file.
		void customPrintToFile(std::string pFileName, bool polarCoords, bool pWeighted,
			bool pGNUplotFormat, const size_t pGridSizeR, const size_t pGridSizeTheta,
			const coeffType pDeltaR, const coeffType pDeltaTheta);
		//void PrintEigenVector(coeffType* pSource, size_t pCol, std::ostream& pTarget = std::cout);

		//Allows for early release of hamiltonian matrix.

		//Calculates energy of wavefunction with n orbitals based on Hartree-Fock method using Slater-Condon rules.
		coeffType calcHFEnergy(size_t pOrbitalCount);

	private:
		Hamiltonian mHamil;
		const bool symmetric, binaryOutput;
		const long gridSizeR = 0;
		const long gridSizeTheta = 0;
		const long matrixSize;
		const long matrixDim;
		const constType dThetaRad;

		coeffType *real, *imag, *work, *tauOut;

		//inline coeffType fastMagSquared(const std::complex<double> in) // use std::norm()?
		//{
		//	return in.imag()*in.imag() + in.real()*in.real();
		//}

		//wavefunction storage (may be able to use only one psi -> gauss-seidel?)
		coeffType *psiOld, *psiNew;

		//USING I AS INDEX FOR R AND J AS INDEX FOR THETA!!
		inline coeffType GetOld(unsigned radial_i, unsigned polar_j)
		{
			return psiOld[polar_j + radial_i * gridSizeTheta];
		}

		//public: //TODO: remove after debugging

		//Returns area for a radial index without prefactor of pi*dr^2/gridSizeTheta.
		inline coeffType getReducedArea(int pRadialIdx, int pInterpolatedIntervals)
		{
			coeffType upperBound = 0, lowerBound = 0;
			int interpolatedGridSize = (gridSizeR - 1)*pInterpolatedIntervals + 1;
			
			if (pRadialIdx == 0) //lower boundary is zero and does not have to be calculated
			{
				upperBound = (0.5 + 0.5 / pInterpolatedIntervals);
			}
#ifndef NDEBUG
			else if (pRadialIdx < 0 || pRadialIdx >= interpolatedGridSize) //invalid index
			{
				throw std::invalid_argument("Radial index out of range.");
			}
#endif // DEBUG
			else if (pRadialIdx == interpolatedGridSize - 1) //upper boundary is always the same (we have to preserve total area)
			{
				upperBound = gridSizeR;
				lowerBound = gridSizeR - 0.5 - 0.5 / pInterpolatedIntervals;
			}
			else
			{
				upperBound = 0.5 + (0.5 + pRadialIdx) / pInterpolatedIntervals;
				lowerBound = 0.5 + (0.5 + pRadialIdx - 1) / pInterpolatedIntervals ;
			}
			return upperBound * upperBound - lowerBound * lowerBound;
		}

	public:
		//Renozmalizes specified target Vector, defaults to psiOld.
		void Renormalize(coeffType* pTarget = nullptr);

		//Prints reference 2D polar harmonic functions up to the max main quantum number.
		void PrintReferenceHarmonics(long maxQuantumNumber);
	};
}

