#include "stdafx.h"
#include "Solver.h"
#include <mkl.h>
#include <algorithm>
//#include <iostream>
//#include <fstream>
namespace PA
{

	Solver::Solver(bool pSymmetric, long pGridSizeR, long pGridSizeTheta, bool pBinaryOutput)
		: symmetric(pSymmetric), binaryOutput(pBinaryOutput), gridSizeR(pGridSizeR), gridSizeTheta(pGridSizeTheta),
		mHamil(pGridSizeR,pGridSizeTheta), matrixDim(pGridSizeR*pGridSizeTheta),
		matrixSize(pGridSizeR*pGridSizeR*pGridSizeTheta*pGridSizeTheta), dThetaRad(2 * M_PI / gridSizeTheta)
	{
		psiNew = new coeffType[matrixDim]();
		memset(psiNew, 0, matrixDim * sizeof(coeffType));
		psiOld = new coeffType[matrixDim]();
		memset(psiOld, 0, matrixDim * sizeof(coeffType));
		real = new coeffType[matrixDim]();
		imag = new coeffType[matrixDim]();
		tauOut = new coeffType[matrixDim]();
		work = nullptr;
	}


	Solver::~Solver()
	{
		if (real != nullptr)
			delete[] real;
		if (imag != nullptr)
			delete[] imag;
		if (work != nullptr)
			delete[] work;
		if (tauOut != nullptr)
			delete[] tauOut;
		if (psiNew != nullptr)
			delete[] psiNew;
		if (psiOld != nullptr)
			delete[] psiOld;
	}

	void Solver::InitPsi()
	{
		for (size_t i = 0; i < gridSizeR; i++)
		{
			coeffType temp = std::exp(-(dR * i + radial_start)*(dR * i + radial_start)*coeffType(.5));
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				psiOld[j + i * gridSizeTheta] = temp; //testing gauss peak
			}
		}
		Renormalize();
		customPrintToFile("InitialPsi.dat",false,true,true,
			gridSizeR, gridSizeTheta, dR, dThetaRad);

	}

	void Solver::InitHamiltonian()
	{
		if (symmetric)
			mHamil.SymmetricSetup();
		else
			mHamil.AsymmetricSetup();
#ifndef NDEBUG
		mHamil.Print();
/*		std::cout << "\nNow in banded Storage!\n";
		mHamil.SymmetricSetup(true);
		mHamil.Print()*/;
#endif // !NDEBUG

	}

	void Solver::initWork(const size_t pEigenvectorCount)
	{
		work = new coeffType[pEigenvectorCount*matrixDim]();
	}

	void Solver::Solve(size_t pEigenVecCount, bool pAlternative)
	{
		if (pEigenVecCount > matrixDim) pEigenVecCount = matrixDim;

		InitHamiltonian();

		//initializing Work for all eigenvectors
		//it is necessary for the lapack algorithms
		initWork(matrixDim);

		if (pAlternative)
			if (symmetric)
				SymmetricDriverRoutineSolve(pEigenVecCount);
			else
				DriverRoutineSolve(pEigenVecCount);
		else
			if(symmetric)
				SymmetricSolve(pEigenVecCount);
			else
				GeneralSolve(pEigenVecCount);

		std::vector<valIdx> desiredEigenVecs = std::vector<valIdx>();
		desiredEigenVecs.reserve(pEigenVecCount);

		FindLowestEigenVecs(pEigenVecCount, desiredEigenVecs, real);

		PrintEV(desiredEigenVecs);

		if (pAlternative)
		{
			std::swap(work, mHamil.mData); 
		}

		if (symmetric)
		{
			//Symmetric solver uses relatively robust representations algorithm, no gram-schmidt necessary.
			NormalizeEigenVecs(pEigenVecCount, work);
			//PrintOrthoRelations(pEigenVecCount, desiredEigenVecs, false);
		}
		else
			//orthogonalization also necessary in asymmetric case
			Orthonormalize(work, desiredEigenVecs, pEigenVecCount);

		//For debugging
		//PrintOrthoRelations(pEigenVecCount, desiredEigenVecs, false);

		//printing eigenvectors to file
		if (binaryOutput)
			printEigenVecsBinary(matrixDim, work, &desiredEigenVecs);
		else
			PrintEigenVecsAsText(matrixDim, work, &desiredEigenVecs,true);
		//mHamil.Print(myfile, [](coeffType val) {return val * val; });
		if (pAlternative)
			std::swap(work, mHamil.mData);
	}

	void Solver::LoadEigenVecs(long pCount, bool fromBinary, std::string pFileName)
	{
		if (pFileName.empty()) pFileName = getEVecFilename(pCount, fromBinary);

		using namespace std::filesystem;
		path targetPath = absolute("Results");
		targetPath.append(pFileName);
		if (fromBinary && !(targetPath.extension() == ".bdata"))
			throw std::invalid_argument("Invalid file");
		if (!fromBinary && !(targetPath.extension() == ".txt"))
			throw std::invalid_argument("Invalid file");
		
		std::ifstream targetFile;
		if(fromBinary) targetFile = std::ifstream(targetPath.native(), std::ios::binary | std::ios::in);
		else targetFile = std::ifstream(targetPath.native(), std::ios::in);

		if (targetFile.is_open())
		{
			initWork(pCount);
			if (!fromBinary)
			{
				std::string line = std::string();
				long elements = 0;
				long eigenVecIdx = 0;
				for (long rows = 0; rows < (gridSizeR + 1)*gridSizeTheta; rows++)
				{
					coeffType trash = 0;
					targetFile >> trash >> trash;
					if (rows % gridSizeR == 0 && rows != 0) //we have a duplicate row for gnuplot plotting
					{
						for (long elem = 0; elem < pCount; elem++)
						{
							targetFile >> trash;
						}
					}
					else //we have unique eigenvector data
					{
						for (long column = 0; column < pCount; column++)
						{
							targetFile >> work[column*matrixDim + eigenVecIdx];
							elements++;
						}
						eigenVecIdx++;
					}
				}

				std::cout << "Number of elements read: " << elements << "\n";
				std::cout << "Number of eigenvectors read: " << (elements + 1) / matrixDim << "\n";
			}
			else
			{
				targetFile.read((char*)work, sizeof(coeffType)*matrixDim*pCount);
			}
		}
		else
			std::cout << "Could not open file: \"" << pFileName << "\"\n";

		targetFile.close();
	}

	void Solver::loadEigenValues(size_t pCount, std::string pFileName)
	{
		if (pFileName.empty()) pFileName = getEValFilename();

		std::string line = std::string();
		long idx = 0;
		using namespace std::filesystem;
		path targetPath = absolute("Results");
		targetPath.append(pFileName);
		std::ifstream targetFile(targetPath.native(), std::ios::in);
		if (targetFile.is_open())
		{
			long idx = 0;
			for (long rows = 0; rows < pCount; rows++)
			{
				std::string trash = std::string();
				targetFile >> trash;
				targetFile >> trash;
				targetFile >> trash;
				targetFile >> real[idx++];
			}

			std::cout << "Number of eigenvalues read: " << idx << "\n";
		}
		else
			std::cout << "Could not open file: \"" << pFileName << "\"\n";

		targetFile.close();
	}

	void Solver::convertData(long pCount, bool fromBinary, std::string pFileName)
	{
		LoadEigenVecs(pCount, fromBinary, pFileName);

		std::vector<valIdx> targetColumns = std::vector<valIdx>();
		for (size_t i = 0; i < pCount && i < matrixDim; i++)
		{
			targetColumns.push_back({ real[i],i }); //initializing desiredEigenVecs
		}

		//transposing from row major eigenvectors this time, so last argument == false
		PrintEigenVecsAsText(pCount, work, &targetColumns, false);
	}

	PA::coeffType Solver::innerProductRowMajor(PA::coeffType * work, size_t vec1, size_t vec2)
	{
		PA::coeffType result = 0;
		size_t offset1 = matrixDim * vec1, offset2 = matrixDim * vec2;
		for (size_t i = 0; i < gridSizeR; i++)
		{
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				size_t currentRowOffset = (gridSizeTheta*i + j);
				result += (2 * i + 1)*work[currentRowOffset+offset1] * work[currentRowOffset+offset2];// norm of psi squared i.e. |psi|� = 1
			}
		}
		return result *dArea / gridSizeTheta;
	}

	//projection of vec1 onto vec2, modifies vec2
	void Solver::ProjectRowMajor(PA::coeffType * work, size_t vec1, size_t vec2)
	{
		//projection <u,v>/<u,u>*u where <a,b> is the dot product between a and b
		PA::coeffType udotv = innerProductRowMajor(work, vec1, vec2);//Should be normalized already / innerProductRowMajor(work, vec2, vec2);
		//PA::coeffType norm = 0;
		for (size_t i = 0; i < gridSizeR; i++)
		{
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				size_t currentRowOffset = (gridSizeTheta*i + j);
				work[currentRowOffset + matrixDim * vec2] -= udotv * work[currentRowOffset + matrixDim * vec1];// u2 = v1 - projection of u1 onto v1
			}
		}
	}

	PA::coeffType Solver::innerProductColMajor(PA::coeffType* work, size_t vec1, size_t vec2)
	{
		PA::coeffType result = 0;
		size_t offset1 = vec1, offset2 = vec2;
		for (size_t i = 0; i < gridSizeR; i++)
		{
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				size_t currentOffset = (gridSizeTheta * i + j) * matrixDim;
				result += (2 * i + 1) * work[currentOffset + offset1] * work[currentOffset + offset2];// norm of psi squared i.e. |psi|� = 1
			}
		}
		return result * dArea / gridSizeTheta;
	}

	//projection of vec1 onto vec2, modifies vec2
	void Solver::ProjectColMajor(PA::coeffType* work, size_t vec1, size_t vec2)
	{
		//projection <u,v>/<u,u>*u where <a,b> is the dot product between a and b
		PA::coeffType udotv = innerProductColMajor(work, vec1, vec2);//Should be normalized already / innerProductRowMajor(work, vec2, vec2);
		//PA::coeffType norm = 0;
		for (size_t i = 0; i < gridSizeR; i++)
		{
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				size_t currentRowOffset = (gridSizeTheta * i + j) * matrixDim;
				work[currentRowOffset + vec2] -= udotv * work[currentRowOffset + vec1];// u2 = v1 - projection of u1 onto v1
			}
		}
	}

	void Solver::Orthonormalize(PA::coeffType * work, std::vector<valIdx> & desiredVecs, const size_t & pEigenVecCount)
	{
		coeffType norm = 0;
		size_t vecCount = desiredVecs.size();
		for (size_t i = 0; i < vecCount; i++) //start with second vector for gram-schmidt
		{
			for (size_t j = 0; j < i; j++) //previous vectors
			{
				ProjectColMajor(work, desiredVecs[j].second, desiredVecs[i].second);
			}
			PA::coeffType norm = sqrt(innerProductColMajor(work, desiredVecs[i].second, desiredVecs[i].second));
			for (size_t k = 0; k < matrixDim; k++)
			{
				work[k*matrixDim + desiredVecs[i].second] /= norm;// norm of psi squared i.e. |psi|� = 1
			}
		}
	}

	void Solver::NormalizeEigenVecs(const size_t &pEigenVecCount, PA::coeffType * work)
	{
		std::vector<coeffType> norms = std::vector<coeffType>(pEigenVecCount);
		for (size_t i = 0; i < gridSizeR; i++)
		{
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				size_t currentRowOffset = (gridSizeTheta*i + j)*matrixDim;
				for (size_t k = 0; k < pEigenVecCount; k++)
				{
					norms[k] += (2 * i + 1)*work[currentRowOffset + k] * work[currentRowOffset + k];// norm of psi squared i.e. |psi|� = 1
				}
			}
		}

		for (size_t k = 0; k < pEigenVecCount; k++)
		{
			norms[k] *= dArea / gridSizeTheta; //coefficients of norm (pi*dr*dr/gridSizeTheta)
		}
		

		for (size_t i = 0; i < gridSizeR; i++)
		{
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				size_t currentRowOffset = (gridSizeTheta*i + j)* matrixDim;
				for (size_t k = 0; k < pEigenVecCount; k++)
				{
					work[currentRowOffset + k] /= sqrt(norms[k]);// norm of psi squared i.e. |psi|� = 1
				}
			}
		}
	}

	void Solver::FindLowestEigenVecs(const size_t &pEigenVecCount, std::vector<valIdx> &desiredEigenVecs, PA::coeffType * real)
	{
		for (size_t i = 0; i < pEigenVecCount && i < matrixDim; i++)
		{
			desiredEigenVecs.push_back({ real[i],i }); //initializing desiredEigenVecs
		}

		struct {
			bool operator()(valIdx left, valIdx right) const
			{
				return left.first < right.first;
			}
		} customLess;

		std::sort(desiredEigenVecs.begin(), desiredEigenVecs.end(), customLess);

		for (size_t i = pEigenVecCount; i < matrixDim; i++)
		{
			coeffType val = real[i];
			if (val <= 0) continue;
			for (size_t j = 0; j < pEigenVecCount; j++)
			{
				if (val < desiredEigenVecs[j].first)
				{
					desiredEigenVecs.insert(desiredEigenVecs.begin() + j, { val, i });
					desiredEigenVecs.resize(pEigenVecCount);//should "forget" the last 
					break;
				}
			}
		}
	}

	void Solver::PrintEV(std::vector<PA::valIdx> &desiredEigenVecs)
	{
		using namespace std::filesystem;

		std::cout << "Printing eigenvalues:\n";
		path filePath = absolute("Results");//current_path();
		create_directory(filePath);

		filePath.append(getEValFilename());
		std::ofstream myfile = std::ofstream();
		myfile.open(filePath.native());
		if (myfile.is_open())
		{
			std::memcpy(imag, real, sizeof(coeffType)*matrixDim);

			if (symmetric)
				std::sort(imag, imag + desiredEigenVecs.size());
			else
				std::sort(imag, imag + matrixDim);
			
			for (int i = 0; i < matrixDim; i++)
			{
				myfile << "Eigenvalue # " << i << ": " << *(imag + i) << "\n";
			}
			myfile.close();
		}
		for (int i = 0; i < desiredEigenVecs.size(); i++)
		{
			std::cout << "Eigenvalue # " << i << ": " << desiredEigenVecs[i].first << "\n";
		}
	}

	void Solver::GeneralSolve(const size_t pEigenVecCount)
	{
		std::cout << "Pursuing general solving scheme for nonsymmetric matrix.\n";

		int ilo = 1, ihi = matrixDim;
		lapack_int info = 0;
		//InitPsi();

		//Hessenberg form
#ifdef DOUBLEPR
		info = LAPACKE_dgehrd(
#else
		info = LAPACKE_sgehrd(
#endif // DOUBLEPR
		LAPACK_ROW_MAJOR, matrixDim, ilo, ihi,
			mHamil.mData, matrixDim, tauOut);
		std::cout << "Reduction to upper Hessenberg form info: " << info << ". (0 means success.)\n";

		memcpy(work, mHamil.mData, matrixSize * sizeof(coeffType));//store hessenberg form in *work for now

		//explicit construction of orthogonal matrix Q from hessenberg reduction
#ifdef DOUBLEPR
		info = LAPACKE_dorghr(
#else
		info = LAPACKE_sorghr(
#endif // DOUBLEPR
		LAPACK_ROW_MAJOR, matrixDim, ilo, ihi, work, matrixDim, tauOut);
		std::cout << "Construction of orthogonal matrix Q info: " << info << ". (0 means success.)\n";

		//eigenvalue calculation -> modifies mData => we need to restore hessenberg form afterwards
#ifdef DOUBLEPR
		info = LAPACKE_dhseqr(
#else
		info = LAPACKE_shseqr(
#endif // DOUBLEPR
		LAPACK_ROW_MAJOR, 'S', 'V', matrixDim, 1, matrixDim,
			mHamil.mData, matrixDim, real, imag, work, matrixDim);
		std::cout << "Eigenvalue calculation info: " << info << ". (0 means success.)\n";

		//std::swap(mHamil.mData, work);

		int numEigenVecs = 0;
		//eigenvectors of matrix in shur form
#ifdef DOUBLEPR
		info = LAPACKE_dtrevc(
#else
		info = LAPACKE_strevc(
#endif // DOUBLEPR
		LAPACK_ROW_MAJOR, 'R', 'B', nullptr, matrixDim, mHamil.mData, matrixDim,
			nullptr, matrixDim, work, matrixDim, matrixDim, &numEigenVecs);
		std::cout << "Eigenvector calculation info: " << info << ". (0 means success.)\n";

	}

	std::filesystem::path Solver::getFilePath(const size_t pEigenVecCount, const bool isBinary)
	{
		using namespace std::filesystem;
		path workingDir = absolute("Results");//current_path();
		//workingDir.append("Results/");
		create_directory(workingDir);
#ifdef DOUBLEPR
		std::string prefix = "d";
#else
		std::string prefix = "s";
#endif // DOUBLEP

		std::string suffix = std::string();
		if (isBinary) suffix = ".bdata";
		else suffix = ".txt";

		std::string symmetry = std::string();
		if (symmetric) symmetry = "sym";
		else symmetry = "asym";

		std::string name = prefix + "_" + symmetry + "_" + getFilenameStem() +
			"_eigenvecs_" + std::to_string(pEigenVecCount) + suffix;
		workingDir.append(name);
		return workingDir;
	}

	void Solver::PrintEVecRow_EVecsColMajor(const size_t &i, const size_t &j, const size_t & pEigenVecCount, std::ofstream & myfile, std::vector<PA::valIdx> * targetColumns, PA::coeffType * work, const size_t &count)
	{
		size_t currentOffset = (gridSizeTheta*i + j)*pEigenVecCount;
		coeffType coord1 = 0, coord2 = 0;
		coord1 = i * dR + radial_start; //radial coord
		coord2 = dThetaRad * j; //in radians

		myfile << coord1 << " " << coord2;
		if (targetColumns == nullptr)
		{
			for (size_t k = 0; k < pEigenVecCount; k++)
			{
				myfile << " " << work[currentOffset + k];
			}
		}
		else
		{
			for (size_t i = 0; i < count; i++)
			{
				size_t idx = (*targetColumns)[i].second;
				myfile << " " << work[currentOffset + idx];
			}
		}
	}

	void Solver::PrintEVecRow_EVecsRowMajor(const size_t& i, const size_t& j, const size_t& pEigenVecCount, std::ofstream& myfile, std::vector<PA::valIdx>* targetRows, PA::coeffType* work, const size_t& count)
	{
		size_t currentOffset = gridSizeTheta * i + j;
		coeffType coord1 = 0, coord2 = 0;
		coord1 = i * dR + radial_start; //radial coord
		coord2 = dThetaRad * j; //in radians

		myfile << coord1 << " " << coord2;
		if (targetRows == nullptr)
		{
			for (size_t k = 0; k < pEigenVecCount; k++)
			{
				myfile << " " << work[currentOffset + k* matrixDim];
			}
		}
		else
		{
			for (size_t i = 0; i < count; i++)
			{
				size_t idx = (*targetRows)[i].second;
				myfile << " " << work[currentOffset + idx* matrixDim];
			}
		}
	}

	void Solver::PrintEigenVecsAsText(const size_t &pEigenVecCount, PA::coeffType * work,
		std::vector<valIdx> *targetColumns, bool pEVecsColMajor, std::ofstream * myfile)
	{
		bool filePtrOwned = false;
		if (myfile == nullptr)
		{
			filePtrOwned = true;
			myfile = new std::ofstream;
			std::filesystem::path filePath = getFilePath(targetColumns->size());
			std::cout << "Writing results to file: " << filePath.filename() << "\n";
			//printing eigenvectors to file
			myfile->open(filePath);
		}
		size_t count = 0;
		if(targetColumns != nullptr)
			count = targetColumns->size();
		
#ifdef DOUBLEPR
		//*myfile << std::scientific << std::setprecision(15);
		*myfile << std::scientific << std::setprecision(std::numeric_limits<double>::digits10);
#else
		*myfile << std::scientific << std::setprecision(6);
#endif // DOUBLEPR
		for (size_t i = 0; i < gridSizeR; i++)
		{
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				if (pEVecsColMajor == true)
					PrintEVecRow_EVecsColMajor(i, j, pEigenVecCount, *myfile, targetColumns, work, count);
				else
					PrintEVecRow_EVecsRowMajor(i, j, pEigenVecCount, *myfile, targetColumns, work, count);
				*myfile << "\n";
			}
			//j = 0, this is for gnuplot only
			if (pEVecsColMajor == true)
				PrintEVecRow_EVecsColMajor(i, 0, pEigenVecCount, *myfile, targetColumns, work, count);//gnuplot extra dataset done
			else
				PrintEVecRow_EVecsRowMajor(i, 0, pEigenVecCount, *myfile, targetColumns, work, count);//gnuplot extra dataset done
			*myfile << "\n\n";
		}

		if (filePtrOwned == true && myfile != nullptr)
		{
			myfile->close();
			std::cout << "Finished writing\n";
			delete myfile;
		}
	}

	void Solver::printEigenVecsBinary(const size_t & pEigenVecCount, PA::coeffType * work, std::vector<valIdx>* targetColumns, std::ofstream * myfile)
	{
		bool filePtrOwned = false;
		if (myfile == nullptr)
		{
			filePtrOwned = true;
			std::filesystem::path filePath = getFilePath(targetColumns->size(),true);
			std::cout << "Writing results to file: " << filePath.filename() << "\n";
			myfile = new std::ofstream(filePath,std::ios::binary);
			//printing eigenvectors to file
			//myfile->open(filePath);
		}
		size_t count = 0;
		if (targetColumns != nullptr)
			count = targetColumns->size();
		coeffType *buffer = new coeffType[count*matrixDim];
		coeffType *psiOldPlaceHolder = nullptr;
		std::swap(psiOld, psiOldPlaceHolder);//save psiOld in placeholder var
		for (size_t i = 0; i < count; i++)
		{
			psiOld = buffer + matrixDim * i;//point psiOld to position in buffer
			ExtractPsiMinor((*targetColumns)[i].second);//extract to psiOld pointer
		}
		myfile->write((char*)(buffer), sizeof(coeffType)*matrixDim*count);
		myfile->close();
		std::swap(psiOld, psiOldPlaceHolder);//recall psiOld from placeholder var
		if (buffer != nullptr) delete[] buffer;
	}

	void Solver::ExtractPsiMinor(long EVTarget)
	{
		for (long i = 0; i < matrixDim; i++)
		{
			psiOld[i] = work[matrixDim*i + EVTarget];
			//psiNew[i] = work[matrixDim*i + EVTarget];
		}
	}

	void Solver::ExtractPsiMajor(long EVTarget)
	{
		for (long i = 0; i < matrixDim; i++)
		{
			psiOld[i] = work[i + matrixDim * EVTarget];
		}
	}

	PA::coeffType Solver::PsiOldPsiNew()
	{
		coeffType result = 0;
		for (long Ri = 0; Ri < gridSizeR; Ri++)
		{
			for (long Tj = 0; Tj < gridSizeTheta; Tj++)
			{
				long currentIdx = gridSizeTheta * Ri + Tj;
				result += psiOld[currentIdx] * psiNew[currentIdx] *(2 * Ri + 1);
			}
		}
		return result * dArea/gridSizeTheta;
	}

	void Solver::DriverRoutineSolve(const size_t pEigenVecCount) //currently working (in comparison to general solve)
	{
		std::cout << "Using driver routine for nonsymmetric matrix (good for error checking).\n";
		InitHamiltonian();

		//initializing Work for all eigenvectors
		//it is necessary for the lapack algorithms
		initWork(matrixDim);
		//InitPsi();

		lapack_int info = 0;
#ifdef DOUBLEPR
		info = LAPACKE_dgeev(
#else
		info = LAPACKE_sgeev(
#endif // DOUBLEPR

		LAPACK_ROW_MAJOR, 'N', 'V', matrixDim, mHamil.mData,
			matrixDim, real, imag, nullptr, matrixDim, work, matrixDim);
		std::cout << "Right eigenvector calculation info: " << info << "\n";

	}

	void Solver::ApplyFofR(const size_t& pEigenVecCount)
	{
		for (size_t ri = 0; ri < gridSizeR; ri++)
		{
			PA::coeffType FofR = 1.0 / sqrt(double(ri) * dR + radial_start);
			for (size_t tj = 0; tj < gridSizeTheta; tj++)
			{
				for (size_t vecCount = 0; vecCount < pEigenVecCount; vecCount++)
				{
					//transform vectors back in symmetric case since we are using psi(r,theta) = f(r)*u(r,theta)
					work[(ri * gridSizeTheta + tj) * matrixDim + vecCount] *= FofR;
				}
			}

		}
	}

	void Solver::SymmetricSolve(const size_t pEigenVecCount)
	{
		std::cout << "Using custom routine for symmetric matrix.\n";
		InitHamiltonian();
		initWork(matrixDim);

		lapack_int info = 0;

		//Reduction to tridiagonal form.
#ifdef DOUBLEPR
		info = LAPACKE_dsytrd( //full storage
#else
		info = LAPACKE_ssytrd(
#endif // DOUBLEPR
			LAPACK_ROW_MAJOR, 'U', matrixDim, mHamil.mData, matrixDim,
			real, imag, tauOut);

		PA::coeffType* evals = new PA::coeffType[matrixDim]();
		int* support = new lapack_int[2*matrixDim]();

		int ev_found = 0;
		lapack_logical high_accuracy = true;
		//Calculation of eigenvalues and eigenvectors.
#ifdef DOUBLEPR
		info = LAPACKE_dstemr(
#else
		info = LAPACKE_sstemr(
#endif // DOUBLEPR
			LAPACK_ROW_MAJOR, 'V', 'I', matrixDim, real, imag, 0, 0, 1, pEigenVecCount,
			&ev_found, evals, work, matrixDim, pEigenVecCount, support, & high_accuracy);

#ifdef DOUBLEPR
		info = LAPACKE_dormtr( //full storage
#else
		info = LAPACKE_sormtr(
#endif // DOUBLEPR
			LAPACK_ROW_MAJOR, 'L', 'U', 'N', matrixDim, matrixDim, mHamil.mData,
			matrixDim, tauOut, work, matrixDim);

		ApplyFofR(pEigenVecCount);

		std::swap(evals, real);

		if (evals != nullptr)
			delete[] evals;
		if (support != nullptr)
			delete[] support;

	}

	void Solver::SymmetricDriverRoutineSolve(const size_t pEigenVecCount)
	{
		std::cout << "Using driver routine for symmetric matrix (good for error checking).\n";
		InitHamiltonian();
		initWork(matrixDim);

		lapack_int info = 0;
#ifdef DOUBLEPR
		info = LAPACKE_dsyevd(
#else
		info = LAPACKE_ssyevd(
#endif // DOUBLEPR
			LAPACK_ROW_MAJOR, 'V', 'U', matrixDim, mHamil.mData, matrixDim, real);

		ApplyFofR(matrixDim);
	}

	void Solver::PrintOrthoRelations(const size_t &pEigenVecCount, std::vector<PA::valIdx> &desiredEigenVecs, const bool isPsiMajor)
	{
		auto results = std::vector<coeffType>(pEigenVecCount*pEigenVecCount);

		bool pWithHamil = true;

		for (long i = 0; i < pEigenVecCount; i++)
		{
			if (isPsiMajor) ExtractPsiMajor(desiredEigenVecs[i].second);
			else ExtractPsiMinor(desiredEigenVecs[i].second); //load into PsiOld

			if (pWithHamil)
				mHamil.ApplyCoeffsToPsi(psiOld, psiNew, PA::Hamiltonian::CoeffApplicationType::AsymHamil); //write into PsiNew
			else
				std::swap(psiNew, psiOld);

			
			for (long j = 0; j < pEigenVecCount; j++)
			{
				if (isPsiMajor) ExtractPsiMajor(desiredEigenVecs[j].second);
				else ExtractPsiMinor(desiredEigenVecs[j].second); //load into PsiOld
				auto buffer = PsiOldPsiNew(); //dot product PsiOld*PsiNew with area weighting
				if (buffer < coeffType(1.e-13))
					results[i + pEigenVecCount * j] = 0;
				else results[i + pEigenVecCount * j] = buffer;

			}
		}

		std::cout << "Printing orthogonality relations:\n";
		std::cout << std::fixed << std::setprecision(5);
		for (long i = 0; i < pEigenVecCount; i++)
		{
			if (pWithHamil)
				std::cout << results[i + pEigenVecCount * i] << " ";
			else
			{
				for (long j = 0; j < pEigenVecCount; j++)
				{
					std::cout << results[i + pEigenVecCount * j] << " ";
				}
			}
			std::cout << "\n";
		}
	}

	void Solver::extendPsi(coeffType * &pSourcePsi, coeffType * &pTargetPsi)
	{
		if (pTargetPsi == nullptr) pTargetPsi = new coeffType[matrixDim + gridSizeTheta * 2];
		std::memcpy(pTargetPsi + gridSizeTheta, pSourcePsi, sizeof(coeffType)*matrixDim);
		for (size_t i = 0; i < gridSizeTheta; i++)
		{
			pTargetPsi[matrixDim + gridSizeTheta + i] = 0; //set outer values for rMax to 0
		}
		for (size_t i = 0; i < gridSizeTheta; i++)
		{
			//pTargetPsi[i] = 2 * pTargetPsi[gridSizeTheta + i] - pTargetPsi[2 * gridSizeTheta + i]; //approximate derivative by extrapolation
			pTargetPsi[i] = pSourcePsi[(gridSizeTheta/2 + i) % gridSizeTheta]; //use values on other side of r=0
		}
	}

	void Solver::interpolatePsi(const size_t &intervalsTheta, const size_t &intervalsR,
		coeffType * &interpolatedPsi, coeffType * &startPsi)
	{
		BicubicInterpolator interpolator = BicubicInterpolator();


		const size_t newGridSizeTheta(gridSizeTheta*intervalsTheta);
		const size_t newGridSizeR((gridSizeR - 1)* intervalsR + 1);
		const size_t newMatrixDim(newGridSizeR*newGridSizeTheta);
		const coeffType newDR(dR / intervalsR), newDTheta(dThetaRad / intervalsTheta);
		//coeffType* coeffStorage = new coeffType[matrixDim*16];
		//if (psiOld != nullptr) delete[] psiOld;
		//psiOld = new coeffType[newMatrixDim];
		//customPrintToFile("oldTest.txt", false, false, true,
		//gridSizeR, gridSizeTheta, dR, dThetaRad);
		//coeffType coeffTest[16] = { 0 }, inputTest[16] = {	//1,3,7,13,
		//													-1,1,5,11,
		//													-1,1,5,11,
		//													1,3,7,13,
		//													5,7,11,17};


		interpolatedPsi = new coeffType[newMatrixDim];
		coeffType* extendedPsi = nullptr;

		extendPsi(startPsi, extendedPsi);

		//for each gridpoint in the eigenvector
		for (long oldRi = 0; oldRi < gridSizeR; oldRi++)
		{
			for (long oldTj = 0; oldTj < gridSizeTheta; oldTj++)
			{
				coeffType inputVals[16] = { 0 };
				//coeffType interpolated[intervalsR*intervalsTheta] = { 0 };

				//determine all surrounding values necessary for interpolation
				for (long y = -1; y < 3; y++)
				{
					for (long x = -1; x < 3; x++)
					{
						long targetRi = oldRi + y, targetTj = oldTj + x;
						if (targetTj >= gridSizeTheta) targetTj -= gridSizeTheta;
						else if (targetTj < 0) targetTj += gridSizeTheta;
						inputVals[(x + 1) + (y + 1) * 4] = extendedPsi[targetRi*gridSizeTheta + targetTj + gridSizeTheta];
					}
				}

				//for (size_t Ri = 0; Ri < 4; Ri++)
				//{
				//	for (size_t Tj = 0; Tj < 4; Tj++)
				//	{
				//		inputVals[Ri * 4 + Tj] = GetOld(Ri, Tj);
				//	}
				//}

				coeffType outputCoeffs[16] = { 0 };
				interpolator.calcDerivativesAndCoeffs(inputVals, outputCoeffs, 1, 1);

				for (size_t Ri = 0; Ri < intervalsR; Ri++)
				{
					for (size_t Tj = 0; Tj < intervalsTheta; Tj++)
					{
						//interpolated[Ri*intervalsTheta + Tj] =
						long newRi = oldRi * intervalsR + Ri, newTj = oldTj * intervalsTheta + Tj;
						if (newRi >= newGridSizeR) continue; //out of bounds, fix later
						interpolatedPsi[newRi*newGridSizeTheta + newTj] =
							interpolator.interpolateBicubic(Tj / coeffType(intervalsTheta), Ri / coeffType(intervalsR), outputCoeffs);
					}
				}
			}
		}
		if (extendedPsi != nullptr) delete[] extendedPsi;
	}

	void Solver::interpolatePsiTest()
	{
		//amount of points we want to add by interpolation
		const size_t intervalsR(3), intervalsTheta(3);
		coeffType* startPsi = psiOld;
		coeffType *interpolatedPsi = nullptr;
		ExtractPsiMajor(0);

		interpolatePsi(intervalsTheta, intervalsR, interpolatedPsi, startPsi);

		//memcpy(interpolatedPsi, psiOld, sizeof(coeffType)*matrixDim);

		std::swap(psiOld, interpolatedPsi);//swap ptrs
		//customPrintToFile("newTest.txt", false, false, true,
			//newGridSizeR, newGridSizeTheta, newDR, newDTheta);
		std::swap(psiOld, interpolatedPsi);//swap back

		if (interpolatedPsi != nullptr) delete[] interpolatedPsi;

		std::cout << "Interpolation done." << std::endl;
		//if (psiOld != nullptr) delete[] psiOld;
		//if (coeffStorage != nullptr) delete[] coeffStorage;
	}

	void Solver::SpecializedSolve(bool pShouldPrint)
	{
		//mHamil.GivensTest();
		//mHamil.SpecializedEigenvalueSolve();
		InitPsi();

		long iteration = 0;
		while(std::abs(psiNew[0]-psiOld[0]) > coeffType(1.e-12) && iteration++ < 10000)
		{
			mHamil.ApplyCoeffsToPsi(psiOld, psiNew, PA::Hamiltonian::CoeffApplicationType::AsymNegCmplxTimeStep);
			std::swap(psiNew, psiOld); //no memcpy, just swapping pointers
			Renormalize(psiOld);
			if(pShouldPrint) customPrintToFile("Iteration" + std::to_string(iteration) + ".dat",false,false,true,
				gridSizeR, gridSizeTheta, dR, dThetaRad);
		}
		std::cout << "Groundstate converged after " << iteration << " iteration(s). Delta is: "
			<< std::abs(psiNew[0] - psiOld[0]) << "\n";

		Renormalize(psiOld);

		customPrintToFile("GroundState.dat", false, false, true,
			gridSizeR,gridSizeTheta,dR,dThetaRad);

	}

	void Solver::CalcEigenValues(bool pSaveData)
	{
		InitHamiltonian();

		//Hessenberg form
		lapack_int info = 0;

#ifdef DOUBLEPR
		info = LAPACKE_dgehrd(
#else
		info = LAPACKE_sgehrd(
#endif // DOUBLEPR
		LAPACK_ROW_MAJOR, matrixDim, 1, matrixDim,
			mHamil.mData, matrixDim, tauOut);
		std::cout << "Reduction to upper Hessenberg form info: " << info << ". (0 means success.)\n";

		//eigenvalue calculation
#ifdef DOUBLEPR
		info = LAPACKE_dhseqr(
#else
		info = LAPACKE_shseqr(
#endif // DOUBLEPR
		LAPACK_ROW_MAJOR, 'E', 'N', matrixDim, 1, matrixDim,
			mHamil.mData, matrixDim, real,imag,nullptr, matrixDim);

		std::cout << "Eigenvalue calculation info: " << info << ". (0 means success.)\n";

		std::sort(real, real + matrixDim);
		//for (size_t i = 0; i < gridSizeR*gridSizeTheta; i++)
		for (size_t i = 0; i < 50 && i < matrixDim; i++)
		{
			coeffType out = 0;
			out = real[i];
			std::cout << "Eigenvalue #" << i << " : (" << out << "," << imag[i] << ")\n";
		}

		if (pSaveData)
		{
			std::ofstream myfile;
			std::string fileName = std::to_string(gridSizeR) + "x" + std::to_string(gridSizeTheta) + "_dr_equals_" + std::to_string(dR) + ".txt";
			std::cout << "Writing results to file: \"" << fileName << "\"\n";
			myfile.open(fileName);
			for (size_t i = 0; i < matrixDim; i++)
			{
				myfile << "Eigenvalue #" << i << " : (" << real[i] << "," << imag[i] << ")\n";
			}
			myfile.close();
			std::cout << "Finished writing\n";
		}
	}

	coeffType Solver::integralNMPQ(coeffType* psiN, coeffType* psiM, coeffType* psiP, coeffType* psiQ,
		size_t pGridSizeR, size_t pGridSizeTheta)
	{
		//we want to calculate the slater-condon integral for a two body operator -> see Wikipedia for details
		//calculation essentially corresponds to element [nm|pq] -> thus the name integralNMPQ

		//TODO: consider interpolation
		const long interpolationIntervals = (pGridSizeR - 1) / (gridSizeR - 1);
		const coeffType newDR(dR*gridSizeR / coeffType(pGridSizeR));
		const coeffType newDTheta(dThetaRad*gridSizeTheta / coeffType(pGridSizeTheta));
		const coeffType newDArea(M_PI * newDR*newDR);

		coeffType result = 0;
		for (size_t ri = 0; ri < pGridSizeR; ri++) //radial coordinate for idx I
		{
			for (size_t ti = 0; ti < pGridSizeTheta; ti++) //theta coordinate for idx I
			{
				for (size_t rj = 0; rj < pGridSizeR; rj++) //radial coordinate for idx J
				{
					for (size_t tj = 0; tj < pGridSizeTheta; tj++) //theta coordinate for idx J
					{
						size_t idxI = ri * pGridSizeTheta + ti , idxJ = rj * pGridSizeTheta + tj; //true index into the wavefunction-array

						//first get all psi values
						coeffType pN = psiN[idxI], pM = psiM[idxJ], pP = psiP[idxI], pQ = psiQ[idxJ];
						//calculate area elements from indices (coefficients multipled to result at end)
						//coeffType areaI = 2 * ri + 1, areaJ = 2 * rj + 1;
						coeffType areaI = getReducedArea(ri, interpolationIntervals), areaJ = getReducedArea(rj, interpolationIntervals);
						if (idxI == idxJ)
						{
							//step 1: generate coefficients for interpolation
							//coeffType localPsiN[25] = { 0 }, localPsiM[25] = { 0 }, localPsiP[25] = { 0 }, localPsiQ[25] = { 0 }; //gridpoints for interpolation
							//for (size_t m = 0; m < 5; m++)
							//{

							//}
							//step 2: interpolate with an odd number of intervals, so there is no overlap between gridpoints (and their respective areas)
							//step 2.5: write new gridpoints with psi values into local array ready for calculation
							//step 3: calculate contribution to integral with correct (!) weights
							result += pN * pM * 2 * std::pow(areaI, coeffType(1.5)) * std::pow(M_PI, coeffType(0.5)); //approximation of contribution around singularity

							//result += 0;	//integral diverges for ri == rj since integrand is proportional to 1/(ri-rj)
											//could potentially interpolate only in this scenario
						}
						else
						{
							//calculate true radius and theta from indices
							coeffType radiusI = ri*newDR + 0.5*dR, radiusJ = rj*newDR + 0.5*dR, thetaI = ti* newDTheta, thetaJ = tj* newDTheta;
							//calculate distance between two points
							coeffType rIJ = sqrt(radiusI*radiusI + radiusJ*radiusJ - 2*radiusI*radiusJ*cos(thetaI - thetaJ));
							//finally calculate contribution to full integral
							result += (pN * pM * pP * pQ * areaI * areaJ)/rIJ;
						}
					}
				}
			}
		}
		//return result * newDArea / coeffType(pGridSizeTheta) *newDArea / coeffType(pGridSizeTheta);
		return result * dArea / coeffType(pGridSizeTheta) *dArea / coeffType(pGridSizeTheta);
	}

	coeffType Solver::slaterCondonPsiCoulombPsi(std::vector<PA::valIdx> &desiredEigenVecs)
	{
		coeffType result = 0;
		size_t orbitals = desiredEigenVecs.size();

		auto hartree_handles = std::vector<std::future<coeffType>>();
		auto fock_handles = std::vector<std::future<coeffType>>();
		auto func = [this](coeffType* psi1, coeffType* psi2, coeffType* psi3, coeffType* psi4, const long pGridSizeR, const long pGridSizeTheta)
		{
			return integralNMPQ(psi1, psi2, psi3, psi4, pGridSizeR, pGridSizeTheta);
		};

		for (size_t i = 0; i < orbitals; i++)
		{
			for (size_t j = 0; j < orbitals; j++)
			{
				//if (j == i) continue;
				//else
				//{
					coeffType *psiI = work + desiredEigenVecs[i].second*matrixDim, *psiJ = work + desiredEigenVecs[j].second*matrixDim;
					//coeffType *interpolatedPsiI = nullptr, *interpolatedPsiJ = nullptr;
					
					//size_t intervalsX = 2, intervalsY = 2;
					//size_t intervalsX = 1, intervalsY = 1;

					//interpolatePsi(intervalsX, intervalsX, interpolatedPsiI, psiI);
					//interpolatePsi(intervalsX, intervalsX, interpolatedPsiJ, psiJ);
					//const size_t newGridSizeTheta(gridSizeTheta*intervalsX);
					//const size_t newGridSizeR((gridSizeR - 1)* intervalsY + 1);
					//auto result1 = integralNMPQ(psiI, psiJ, psiI, psiJ, gridSizeR, gridSizeTheta);
					//auto result2 = integralNMPQ(psiI, psiJ, psiJ, psiI, gridSizeR, gridSizeTheta);

					//For asynchronous work. An "easy" way to multithread independent calculations.
					hartree_handles.push_back(std::async(func, psiI, psiJ, psiI, psiJ, gridSizeR, gridSizeTheta));
					fock_handles.push_back(std::async(func, psiI, psiJ, psiJ, psiI, gridSizeR, gridSizeTheta));

					//auto result1 = integralNMPQ(interpolatedPsiI, interpolatedPsiJ, interpolatedPsiI, interpolatedPsiJ,
					//	newGridSizeR, newGridSizeTheta);
					//auto result2 = integralNMPQ(interpolatedPsiI, interpolatedPsiJ, interpolatedPsiJ, interpolatedPsiI,
					//	newGridSizeR, newGridSizeTheta);
					//std::cout << result1 << " " << result2 << std::endl;
					//result += 2*result1 - result2;
					//result += integralNMPQ(interpolatedPsiI, interpolatedPsiJ, interpolatedPsiI, interpolatedPsiJ,
					//	newGridSizeR, newGridSizeTheta) -
					//	integralNMPQ(interpolatedPsiI, interpolatedPsiJ, interpolatedPsiJ, interpolatedPsiI,
					//		newGridSizeR, newGridSizeTheta);

					//if (interpolatedPsiI != nullptr) delete[] interpolatedPsiI;
					//if (interpolatedPsiJ != nullptr) delete[] interpolatedPsiJ;
				//}
			}
		}

		for (size_t i = 0; i < fock_handles.size(); i++)
		{
			//resynchronization point for asynchronously executed tasks/calculations from before.
			auto res1 = hartree_handles[i].get(), res2 = fock_handles[i].get();
			std::cout << res1 << " " << res2 << std::endl;
			result += 2 * res1 - res2;
		}
		//result *= coeffType(0.5);	//consider removing this -> is (i,j) == (j,i) in the loop?
									//then it may be possible to reduce the work done in the previous for statement
									//have to choose indices correctly though

		return result;
	}

	coeffType Solver::calcHFEnergy(size_t pOrbitalCount)
	{
		coeffType result = 0;
		if (pOrbitalCount < 1) pOrbitalCount = 1; //no less than 1 orbital allowed

		std::vector<valIdx> eigenVectors = std::vector<valIdx>();
		if (real != nullptr)
			FindLowestEigenVecs(pOrbitalCount, eigenVectors, real);
		else
			throw std::logic_error("Branch not implemented.");

		//PrintOrthoRelations(pOrbitalCount, eigenVectors,true);

		for (size_t i = 0; i < pOrbitalCount; i++)
		{
			result += 2* eigenVectors[i].first;//energies of individual wave functions
		}

		coeffType temp = 0;
		if (pOrbitalCount > 1) //need at least two orbitals for coulomb correction
		{
			temp += slaterCondonPsiCoulombPsi(eigenVectors); //coulomb correction
		}

		return result + temp;
	}

	///Prints in format r,theta (in degrees), psi�
	void Solver::customPrintToFile(std::string pFileName, bool pPolarCoords, bool pWeighted,
		bool pGNUplotFormat, const size_t pGridSizeR, const size_t pGridSizeTheta,
		const coeffType pDeltaR, const coeffType pDeltaTheta)
	{
		coeffType areaElement = M_PI * pDeltaR*pDeltaR;
		auto out = std::ofstream(pFileName);
#ifdef DOUBLEPR
		out << std::scientific << std::setprecision(15);
#else
		out << std::scientific << std::setprecision(6);
#endif // DOUBLEPR
		if (out.is_open())
		{
			for (size_t i = 0; i < pGridSizeR; i++)
			{
				for (size_t j = 0; j < pGridSizeTheta; j++)
				{
					coeffType psi = 0;
					if (pWeighted)
					{
						psi = (2 * i + 1)*GetOld(i, j)*GetOld(i, j)*areaElement;
						//psi /= M_PI * gridSizeR*gridSizeR;
					}
					else psi = psiOld[i*pGridSizeTheta + j];//GetOld(i, j);// *GetOld(i, j);
					coeffType coord1 = 0, coord2 = 0;
					if (pPolarCoords || pGNUplotFormat)
					{
						coord1 = i* pDeltaR + pDeltaR / 2.0; //radial coord
						coord2 = pDeltaTheta * j; //in radians
					}
					else
					{
						coeffType radius = i* pDeltaR + pDeltaR / 2.0, theta = pDeltaTheta * j;
						coord1 = radius * cos(theta);
						coord2 = radius * sin(theta);
					}
					out /*<< std::fixed << std::setprecision(9)*/ << coord1 << " " << coord2 << " " << psi << "\n";
				}
				if (pGNUplotFormat)
				{
					size_t j = 0;
					coeffType psi = 0;
					if (pWeighted)
					{
						psi = (2 * i + 1)*GetOld(i, j)*GetOld(i, j)*areaElement;
						//psi /= M_PI * gridSizeR*gridSizeR;
					}
					else psi = psiOld[i*pGridSizeTheta + j];//GetOld(i, j);// *GetOld(i, j);
					coeffType coord1 = 0, coord2 = 0;
					coord1 = i * pDeltaR + pDeltaR/2.0; //radial coord
					coord2 = pDeltaTheta * j; //in radians
					out /*<< std::fixed << std::setprecision(9)*/ << coord1 << " " << coord2 << " " << psi << "\n";

					out << "\n";//spacing for pm3d map
				}
			}
			out.close();
		}

	}

	void Solver::Renormalize(coeffType* pTarget)
	{
		if (pTarget == nullptr)
			pTarget = psiOld;

		double sum = 0.f;
		//simple integral (for now)
		for (size_t i = 0; i < gridSizeR; i++)
		{
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				//Weighting by area -> pi*(rmax**2 - rmin**2)/thetasegments = area per point
				sum += (2 * i + 1)*pTarget[j + i * gridSizeTheta] * pTarget[j + i * gridSizeTheta];
			}
		}
		sum *= dArea/gridSizeTheta;
		//sum /= M_PI * gridSizeR*gridSizeR; //normalzing weights to 1;
#ifndef NDEBUG
		std::cout << "Sum is: " << sum << "\n"; //for debugging
#endif // !NDEBUG
		sum = std::sqrt(sum);
		double debugSum = 0.f;
		for (size_t i = 0; i < gridSizeR; i++)
		{
			for (size_t j = 0; j < gridSizeTheta; j++)
			{
				psiOld[j + i * gridSizeTheta] /= sum;
#ifndef NDEBUG
				debugSum += (2 * i + 1)*pTarget[j + i * gridSizeTheta] * pTarget[j + i * gridSizeTheta];
#endif // !NDEBUG
			}
		}
#ifndef NDEBUG
		debugSum *= dArea/ gridSizeTheta;
		//debugSum /= M_PI * gridSizeR*gridSizeR;
		std::cout << "Integration over entire space of <psi|psi> should be 1. Result: " << debugSum << "\n";
#endif // !NDEBUG
	}

	void Solver::PrintReferenceHarmonics(long maxQuantumNumber)
	{
		std::cout << "Writing out first " << maxQuantumNumber * (maxQuantumNumber + 1) / 2 << " eigenfunctions.\n";
		//for (long i = 0; i < maxQuantumNumber; i++)
		//{
		//	for (long j = 0; j <= i; j++)
		//	{
		//		long NX = i-j;
		//		long NY = j;
		//		std::ofstream myfile;
		//		std::string fileName = "RefEV" + std::to_string(i) + "NX" + std::to_string(NX) + "NY" + std::to_string(NY) + ".dat";
		//		myfile.open(fileName);
		//		mHamil.Print2DHarmonicCart(myfile, NX, NY);
		//		myfile.close();
		//	}
		//}

		for (long E = 1; E <= maxQuantumNumber; E++)
		{
			for (long m = 0; m < E; m++)
			{
				if ((E - m - 1) >= 0 && fmod((E - m - 1) / 2.f, 1.f) == 0.0)
				{ //we have valid quantum numbers
					long k = (E - m - 1) / 2;
					if (m == 0)
					{
						std::ofstream myfile;
						std::string fileName = "RefEV" + std::to_string(E) + "K"
							+ std::to_string(k) + "M" + std::to_string(0) + ".dat";
						myfile.open(fileName);
						mHamil.Print2DHarmonicPolar(myfile, k, 0);
						myfile.close();
					}
					else
					{
						std::ofstream myfile;
						std::string fileName = "RefEV" + std::to_string(E) + "K"
							+ std::to_string(k) + "M" + std::to_string(m) + ".dat";
						myfile.open(fileName);
						mHamil.Print2DHarmonicPolar(myfile, k, m);
						myfile.close();

						fileName = "RefEV" + std::to_string(E) + "K"
							+ std::to_string(k) + "M-" + std::to_string(m) + ".dat";
						myfile.open(fileName);
						mHamil.Print2DHarmonicPolar(myfile, k, -m);
						myfile.close();
					}
				}
			}
		}



		std::cout << "Finished writing\n";
	}
}
