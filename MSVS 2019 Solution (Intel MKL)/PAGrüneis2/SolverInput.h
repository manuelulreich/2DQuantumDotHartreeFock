#pragma once
#include <string>
#include <iostream>
namespace PA
{
	template<typename T>
	void RetryInput(T& inputTarget)
	{
		std::cin.clear();
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin >> inputTarget;
	}

	class SolverInput
	{
	public:
		long gridSizeR = 0, gridSizeTheta = 0, evCount = 0, orbitalCount = 0, solverOptionInput = 0;
		std::string approachInput = "",confirmInput = "", binaryInput = "";
		bool symmetric = false , binaryFile = true;
		void hardcodeInput(long pGridSizeR, long pGridSizeTheta, long pEVCount, bool pSymmetric,
			bool pBinaryData, long pOrbitalCount, long pSolverOption);
	private:
		void determineOption();
		bool validateBinaryInput() {
			return (binaryInput == "b") || (binaryInput == "t"); }
		bool validateGridSizeR() {
			return gridSizeR > 0; }
		bool validateGridSizeTheta() {
			return gridSizeTheta > 0; }
		bool validateEvCount() {
			return (evCount > 0) && (evCount <= gridSizeR * gridSizeTheta); }
		bool validateApproachInput() {
			return (approachInput == "s") || (approachInput == "a"); }
		bool validateConfirmInput() {
			return (confirmInput == "y") || (confirmInput == "n"); }
		bool validateOption() {
			return (solverOptionInput > 0) && (solverOptionInput <= 4); }
		bool validateOrbitalCount() { return (orbitalCount > 0) && (orbitalCount <= evCount); }

	public:

		void readInput();
		bool validateInput();

	};
}
