#include "stdafx.h"
#include "SolverInput.h"
namespace PA
{
	void SolverInput::readInput()
	{
		using namespace std;
		bool confirmed = false;
		while (!confirmed)
		{
			determineOption();
			cout << "There are some parameters to enter. You will be asked to confirm entries in the end.\n";
			cout << "1) Number of gridpoints in the radial direction?\n";
			cin >> gridSizeR;
			while (!validateGridSizeR())
			{
				std::cout << "Bad input. Value has to be greater than 0 and an integer. Please Try again.\n";
				RetryInput(gridSizeR);
			}
			cout << "2) Number of gridpoints in the azimuthal direction?\n";
			cin >> gridSizeTheta;
			while (!validateGridSizeTheta())
			{
				std::cout << "Bad input. Value has to be greater than 0 and an integer. Please Try again.\n";
				RetryInput(gridSizeTheta);
			}
			cout << "3) Number of eigenvectors and values to solve for? Max is " << gridSizeR * gridSizeTheta << ".\n";
			cin >> evCount;
			while (!validateEvCount())
			{
				cout << "Bad input. Value has to be greater than 0, less than " << gridSizeR * gridSizeTheta << ", and an integer. Please Try again.\n";
				RetryInput(evCount);
			}
			cout << "4) Symmetric or asymmetric approach? (enter 's' or 'a')\n";
			cin >> approachInput;
			while (!validateApproachInput())
			{
				cout << "Bad input. Value has to be \"s\" or \"a\".\n";
				RetryInput(approachInput);
			}
			if (approachInput == "s")
				symmetric = true;
			else if (approachInput == "a")
				symmetric = false;

			//cout << "5) Save/load eigenvectors from/to binary file or from .txt? Text format is less precise than binary format. (b/t)...\n";
			//cin >> binaryInput;
			//while (!validateBinaryInput())
			//{
			//	cout << "Bad input. Value has to be \"b\" or \"t\".\n";
			//	RetryInput(binaryInput);
			//}
			//if (binaryInput == "b")
			//	binaryFile = true;
			//else if (binaryInput == "t")
			//	binaryFile = false;

			if (solverOptionInput == 2 || solverOptionInput == 3)
			{
				cout << "6) Number of orbitals to calculate HF energy for? Max is " << evCount << ".\n";
				cin >> orbitalCount;
				while (!validateOrbitalCount())
				{
					cout << "Bad input. Value has to be greater than 0, less than " << evCount << ", and an integer. Please Try again.\n";
					RetryInput(orbitalCount);
				}
			}

			cout << "Radial grid points: " << gridSizeR << endl;
			cout << "Azimuthal grid points: " << gridSizeTheta << endl;
			cout << "First n eigenvectors to be saved as output: " << evCount << endl;
			cout << "Symmetric approach: " << boolalpha << symmetric << noboolalpha << endl;
			
			if (solverOptionInput == 2 || solverOptionInput == 3)
			{
				cout << "Number of orbitals to calculate HF energy for: " << orbitalCount << endl;
			}

			cout << "Confirm parameters: (y/n)...\n";
			cin >> confirmInput;
			while (!validateConfirmInput())
			{
				cout << "Bad input. Value has to be \"y\" or \"n\".\n";
				RetryInput(confirmInput);
			}
			if (confirmInput == "y")
				confirmed = true;
		}
	}
	bool SolverInput::validateInput()
	{
		return validateApproachInput() && validateConfirmInput() && validateEvCount() &&
			validateGridSizeR() && validateGridSizeTheta() && validateOption();
	}

	void SolverInput::hardcodeInput(long pGridSizeR, long pGridSizeTheta, long pEVCount, bool pSymmetric, bool pBinaryData, long pOrbitalCount, long pSolverOption)
	{
		gridSizeR = pGridSizeR;
		gridSizeTheta = pGridSizeTheta;
		evCount = pEVCount;
		symmetric = pSymmetric;
		binaryFile = pBinaryData;
		orbitalCount = pOrbitalCount;
		solverOptionInput = pSolverOption;
	}

	void SolverInput::determineOption()
	{
		using namespace std;
		cout << "Options...\n";
		cout << "1 ... Just calculate eigenvalues and -vectors\n";
		cout << "2 ... Calculate eigenvalues and vectors, and HF energy (orbital count entered later on)\n";
		cout << "3 ... Load eigenvalues and vectors from file, then calculate HF energy\n";
		cout << "4 ... Convert binary data to ascii data for plotting/reading\n";
		cout << "Please enter option number:\n";
		cin >> solverOptionInput;
		while (!validateOption())
		{
			cout << "Bad input. Value has to be a listed option.\n";
			RetryInput(solverOptionInput);
		}
	}
}