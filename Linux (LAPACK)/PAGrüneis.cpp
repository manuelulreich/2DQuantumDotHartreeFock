#include "Solver.h"
#include "SolverInput.h"
#include <chrono>
#include <vector>


enum SolverOptions
{
	EVOnly = 1,
	EVandHF = 2,
	HFOnly = 3,
	ConvertData = 4
};

float routine(PA::SolverInput pInput)
{
	
	using namespace std::chrono;
	high_resolution_clock clock = high_resolution_clock();
	auto t1 = clock.now();
	 
	if (pInput.solverOptionInput == ConvertData) //convert binary file to text
	{
		PA::Solver mySolver(pInput.symmetric, pInput.gridSizeR, pInput.gridSizeTheta, pInput.binaryFile);
		mySolver.convertData(pInput.evCount, pInput.binaryFile);
	}
	if (pInput.solverOptionInput == EVOnly || pInput.solverOptionInput == EVandHF) //only ev & ev with hf
	{
		PA::Solver mySolver(pInput.symmetric, pInput.gridSizeR, pInput.gridSizeTheta, pInput.binaryFile);
		mySolver.Solve(pInput.evCount); //calculation of eigenvals + eigenvecs via Intel MKL
	}
	if (pInput.solverOptionInput == EVandHF || pInput.solverOptionInput == HFOnly) //ev with hf & only hf
	{
		PA::Solver mySolver(pInput.symmetric, pInput.gridSizeR, pInput.gridSizeTheta, pInput.binaryFile);
		mySolver.LoadEigenVecs(pInput.evCount, pInput.binaryFile);
		mySolver.loadEigenValues(pInput.gridSizeR*pInput.gridSizeTheta);
		std::cout << mySolver.calcHFEnergy(pInput.orbitalCount) << std::endl;
	}
	
	auto t2 = clock.now();
	auto s = duration_cast<microseconds>(t2 - t1).count() / float(1000000);
	std::cout << "Calculation took: " << s << " seconds.\n";
	return s;
}


int main()
{
	std::cout << "Solver for the 2D quantum harmonic oscillator, with additional hartree fock calcution of multiparticle energies.\n";

	PA::SolverInput input = PA::SolverInput();
	input.readInput();
	routine(input);
	
	///THE FOLLOWING COMMENTS DISPLAY A FEW EXAMPLE FOR SETTING UP MULTIPLE CALCULATIONS VIA CODE
	///Naturally reading the input, and starting the routine above needs to be removed in that case.

	//for (size_t i = 50; i <= 190; i += 10) //symmetric loop
	//{
	//	input.hardcodeInput(i, i, 50, true, true, 3, EVandHF);
	//	routine(input);//calculate values, vectors, hf energy
	//	input.hardcodeInput(i, i, 50, true, true, 3, ConvertData);
	//	routine(input);//convert data from binary to text for plotting
	//}

	//size_t j = 190, k = 160, l = 240;

	//input.hardcodeInput(j, j, 50, false, true, 3, EVandHF);
	//routine(input);//calculate values, vectors, hf energy
	//input.hardcodeInput(j, j, 50, false, true, 3, ConvertData);
	//routine(input);//convert data from binary to text for plotting

	//input.hardcodeInput(j, j, 50, false, true, 3, EVOnly);
	//routine(input);//calculate values, vectors, hf energy
	//input.hardcodeInput(j, j, 50, false, true, 3, ConvertData);
	//routine(input);//convert data from binary to text for plotting
	//input.hardcodeInput(j, j, 50, false, false, 3, HFOnly);
	//routine(input);//calculate values, vectors, hf energy

	//input.hardcodeInput(k, l, 50, false, true, 3, EVOnly);
	//routine(input);//calculate values, vectors, hf energy
	//input.hardcodeInput(k, l, 50, false, true, 3, ConvertData);
	//routine(input);//convert data from binary to text for plotting

	//input.hardcodeInput(k, l, 50, true, true, 3, EVOnly);
	//routine(input);//calculate values, vectors, hf energy
	//input.hardcodeInput(k, l, 50, true, true, 3, ConvertData);
	//routine(input);//convert data from binary to text for plotting
	//std::cout << "Starting symmetric 6 orbital HF calculations:\n";

	//for (size_t i = 180; i <= 180; i += 10) //asymmetric loop
	//{
	//	input.hardcodeInput(i, i, 50, false, false, 6, HFOnly);
	//	routine(input);//calculate values, vectors, hf energy
	//}

	//std::cout << "Starting asymmetric 6 orbital HF calculations:\n";

	//for (size_t i = 180; i <= 180; i += 10) //asymmetric loop
	//{
	//	input.hardcodeInput(i, i, 50, true, false, 6, HFOnly);
	//	routine(input);//calculate values, vectors, hf energy
	//}


#ifndef NDEBUG
	std::cin.get();//wait for input
#endif // NDEBUG

	return 0;
}
