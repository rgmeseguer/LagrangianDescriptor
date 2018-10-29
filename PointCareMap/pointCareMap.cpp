#include <cstdio>
#include <vector>			// std::vector
#include <sstream>			// std::stringstream
#include <iomanip>			// std::setprecision
#include <fstream>			// std::ofstream
#include <ios>

#include "..\LagrangianDescriptor\utilities.h"
#include "..\LagrangianDescriptor\Oscillator.h"
#include "..\LagrangianDescriptor\Dynamics.h"

utilities ut;

/* Potential and Gradient of the system */
double DOS_V(std::vector<double> coeff, std::vector<double> r)
{
	double ep = 0.;
	for (int i = 0; i < 5; ++i)
	{
		ep += coeff[i] * ut.powerd(r[0], i);
	}
	ep += coeff[5] * ut.powerd(r[1] - coeff[6], 2);
	ep += coeff[7] / ut.powerd(r[1] - r[0], 12);

	return ep;
}
std::vector<double> DOS_G(std::vector<double> coeff, std::vector<double> r)
{
	std::vector<double> grad(2, 0.);
	grad[0] = 0;
	for (int j = 1; j < 5; ++j)
	{
		grad[0] += double(j)*coeff[j] * ut.powerd(r[0], j - 1);
	}
	grad[0] += 12. * coeff[7] / ut.powerd(r[1] - r[0], 13);
	grad[1] = 0;
	grad[1] = 2. * coeff[5] * (r[1] - coeff[6]) - 12. * coeff[7] / ut.powerd(r[1] - r[0], 13);
	return grad;
}

int main(int argc, char *argv[])
{
	/* Check the correct Parameters */
#pragma region Program Usage
	if ((argc != 4))
	{
		std::cout << "Program Usage" << std::endl;
		std::cout << "PointCareMap.out [bathMass] [TrjTime] [Trj Number] [Crossing Point]" << std::endl;
		return 1;
	}
#pragma endregion

	/* Variables introduced by the User */
#pragma region Input variables
	double bathMass = strtof(argv[1], NULL);						// Mass of the Bath
	double TrjTime = strtof(argv[2], NULL);							// Time of the Trj
	int I = int(strtof(argv[3], NULL));								// Number of the trajectory
	double crossPoint = strtof(argv[4], NULL);						// Point at wich the crossing trj is saved

#pragma endregion

#pragma region System Initial Conditions

	/* Initiate the Oscillator */
	Oscillator Osc({ 321.904484,-995.713452,1118.689573,-537.856726,92.976121,1.0,1.0,0.01 },	/* Coefficients from the oscillator */\
	{ 1., bathMass },																			/* Mass Values */\
	DOS_V, DOS_G);																				/* Potential and Gradient Functions */

	std::vector<double> R0 = { 1.36561 ,2.161769 };					// Initial Position
	double Energy = 3.691966889;									// Energy of the system
	double timeStep = 1.e-3;										// Set the Time Step/Precision of the Dynamic

#pragma endregion

	return 0;
}