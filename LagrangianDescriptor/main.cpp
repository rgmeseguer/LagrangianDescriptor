#include <cstdio>
#include <vector>			// std::vector
#include <sstream>			// std::stringstream
#include <iomanip>			// std::setprecision
#include <fstream>			// std::ofstream

#include "utilities.h"
#include "LDSurfaceCreator.h"
#include "Oscillator.h"
#include "LagDesc.h"
#include "Dynamics.h"

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
	if ((argc != 5))
	{
		std::cout << "Program Usage" << std::endl;
		std::cout << "LagD_Trj.out [tau] [bathMass] [TrjTime > (2xtau+1)] [Trj Number]" << std::endl;
		return 1;
	}
#pragma endregion

/* Variables introduced by the User */
#pragma region Input variables
	double tau = strtof(argv[1], NULL);								// Tau value
	double bathMass = strtof(argv[2], NULL);						// Mass of the Bath
	double TrjTime = strtof(argv[3], NULL);							// Time of the Trj
	int I = int(strtof(argv[4], NULL));								// Number of the trajectory
#pragma endregion

/* Set Variables for writing the output files */
#pragma region Variables for printing
	std::string BthM, stau;																						//String version of tau and Mass
	std::stringstream stream;																					//Temporary string
	std::ofstream outputFile;																					//Saving file


	/* Print the bath mass in the correct format  */
	if (bathMass < 1) { stream << std::fixed << std::setprecision(1) << bathMass;	BthM = stream.str(); }
	else { stream << std::fixed << std::setprecision(0) << bathMass; BthM = stream.str(); }
	stream.str(std::string());																					//clear the stringstream

	/* Print the Tau in the correct format */
	if (tau < 0) { stream << std::fixed << std::setprecision(1) << tau; stau = stream.str(); }
	else { stream << std::fixed << std::setprecision(0) << tau; stau = stream.str(); }
	stream.str(std::string());																					//clear the stringstream

#if KEY_DETERM==1
	std::string calc = "LDT" + stau + "m" + BthM;																//Name of the calculation
#else
	std::string calc = "sLDT" + stau + "m" + BthM;																//Name of the calculation
#endif


#pragma endregion

/* Set the Initial Conditions */
#pragma region System Initial Conditions

	/* Initiate the Oscillator */
	Oscillator Oscf({ 321.904484,-995.713452,1118.689573,-537.856726,92.976121,1.0,1.0,0.01 },	/* Coefficients from the oscillator */\
	{ 1., bathMass },																/* Mass Values */\
		DOS_V, DOS_G);																/* Potential and Gradient Functions */
	Oscillator Oscb({ 321.904484,-995.713452,1118.689573,-537.856726,92.976121,1.0,1.0,0.01 },	/* Coefficients from the oscillator */\
	{ 1., bathMass },																/* Mass Values */\
		DOS_V, DOS_G);																/* Potential and Gradient Functions */

	srand(time(NULL));												// Initialize the random seed

	std::vector<double> R0 = { 1.36561 ,2.161769 };					// Initial Position
	double Energy = 3.691966889;									// Energy of the system
	double timeStep = 1.e-3;										// Set the Time Step/Precision of the Dynamic

#pragma endregion

/* Initialize the Variables */
#pragma region Process Required Variables
	/* We are going to perform two Dynamics at the same time,
	one in Foward and the other in Backward time.
	To do so, we need to initialize 2 Oscillators and 2 Dynamics*/

	/* Initialize the Oscillators */
	Oscf.setInitP(R0, true);										// Initiate at the shadle point
	Oscf.setInitialVel_NVE(Energy, I);								// Set initial Random Velocities to keep energy
	Oscb.setInitP(R0, true);										// Initiate at the shadle point
	Oscb.setInitialVel_NVE(Energy, I);								// Set initial Random Velocities to keep energy

	/* Initialize the dynamics */
	Dynamics Dynf;													// Set the forward Dynamic							
	Dynf.setTimeStep(timeStep);										// Set the Time step
	Dynf.setTime(TrjTime);											// and the total time (nsteps = totalTime/timeTtep)
	Dynamics Dynb;													// Set the backward Dynamic
	Dynb.setTimeStep(-1 * timeStep);								// Set the Time step
	Dynb.setTime(TrjTime);											// and the total time (nsteps = totalTime/timeTtep)
	bool growingState = true;										// state of the trj growing or updating

	/* Variables fot the Stochastic Dynamics */
#if KEY_DETERM==0													
	double beta = 4.;												// Effective temperature
	double gamma = 5.;												// Disipation factor
	Dyn.setLangevin(Osc.mu.back(), beta, gamma);

	long T = time(0);												// Seed for the random 
	std::vector<double> gasdev = ut.box_muller(Dyn.nstep, &T);		// Gas deviation
	double rtherm;
#endif

	/* Initialize the Lagrangian Descriptor Object*/
	LagDesc LD;														// Lagrangian Descriptor variable

	/* Create the surface of the where the LD values will be saved */
	LDSurfaceCreator Surface;
	std::vector<std::vector<std::string>> openPoints_f;				// Tracker of the points at which the LD forward is beig calculated	
	std::vector<std::vector<std::string>> openPoints_b;				// Tracker of the points at which the LD backward is beig calculated	
	std::vector<std::string> zeroPoint(Oscf._size * 2, " ");		// Initial point

	/* Variables to create the Keys of the surface's points */
	std::vector<std::string> keyf(Oscf._size * 2, " ");				// Key of the position for the forward Trj 
	std::vector<std::string> keyb(Oscf._size * 2, " ");				// Key of the position for the backward Trj
	std::vector<std::string> key(Oscf._size * 2, " ");				// Temporary key


	/* Surface Variables Required to Udate the Points */
	std::vector<std::vector<double>> LDList;					// Buffer of the LD value at each step
	std::vector<double> LDTot = { 0.,0. };						// Buffer of the LD value
	double totLength = 0.;										// Lenght of the Buffers
	int bufferCounter = 0;										// Counter to follow the values we need to remove from the buffer
#pragma endregion

/* Open the inital Point */
#pragma region Open first Point
	key = Surface.keyWrite(Oscf._position, Oscf.calcMomenta());
	Surface.addPoint(key);														// Create the new point
	zeroPoint = key;															// Save the Initial point

#pragma endregion

/* Initate the trayectory */
	for (size_t j = 0; j < Dynf._numberStep; j++)
	{
		/* Once we reach time tau points begin to close opened Points */ 
		if (growingState && j*Dynf._timeStep >= tau) { growingState = false; } 

/* Perform a Dinamic Step in forward and backward direction */
#pragma region Dynamic Step
#if KEY_DETERM==1
		Dynf.DynamicStep(Oscf);										//Dynamic forward step
		Dynb.DynamicStep(Oscb);										//Dynamic backward step
#else
		rtherm = Dynf.sigma*gasdev[j];								//Bath Random
		Dynf.DynamicStep(Oscf, rtherm);								//Dynamic foward step
		rtherm = Dynb.sigma*gasdev[j];								//Bath Random
		Dynb.DynamicStep(Oscb, rtherm);								//Dynamic backward step

#endif
#pragma endregion

/* With the new positions write the two new keys obtained */
#pragma region Key Write
		/* Forward Key Write*/
		keyf = Surface.keyWrite(Oscf._position, Oscf.calcMomenta());

		/* Backward Key Write*/
		keyb = Surface.keyWrite(Oscb._position, Oscb.calcMomenta());
#pragma endregion

/* And Calculate the value of LD at those points */
#pragma region Lagrangian Descriptor Calculation

#if KEY_MF==1
		/* We can use LDTot as a temporal Variable to save the value of the LD at that point */

		/* Foward LD Calculation */
		LDTot = LD.MAction(Oscf.calcMomenta(), Oscf._velocity);
		LDTot[0] *= timeStep;
		LDTot[1] *= timeStep;

		LDList.push_back(LDTot);											// Save it at the end

		/* Backward LD Calculation */
		LDTot = LD.MAction(Oscb.calcMomenta(), Oscb._velocity);
		LDTot[0] *= timeStep;
		LDTot[1] *= timeStep;

		LDList.insert(LDList.begin(), LDTot);							// Save it at the begining

#endif
#pragma endregion

/* At every step Open those points that require to do so */
#pragma region Open Points

	/* Open new point only until reach Time-tau otherwise we will have unfinished points */
		if (j*Dynf._timeStep < Dynf._numberStep*Dynf._timeStep - tau)
		{
#pragma region Forward Opening Points

			/* If the key does not exist create a new one */
			if (Surface.doesKeyExist(keyf) == false)
			{

				Surface.addPoint(keyf);									// Create the new point
				Surface.openPoint(keyf);								// Open the point
				openPoints_f.push_back(keyf);							// Include the point in the list of open points
			}
			/* And if exist but is not open, open the point again */
			else if (Surface._pointStatComplete[keyf] == true)
			{
				Surface.openPoint(keyf);								// Open the point
				openPoints_f.push_back(keyf);							// Include the point in the list of open points
			}

#pragma endregion

#pragma region Backward Opening Points

			/* If the key does not exist create a new one */
			if (Surface.doesKeyExist(keyb) == false)
			{
				Surface.addPoint(keyb);									// Create the new point
				Surface.openPoint(keyb);								// Open the point
				openPoints_b.push_back(keyb);							// Include the point in the list of open points

			}
			/* And if exist but is not open, open the point again */
			else if (Surface._pointStatComplete[keyb] == true)
			{
				Surface.openPoint(keyb);								// Open the point
				openPoints_b.push_back(keyb);							// Include the point in the list of open points
			}

#pragma endregion

		}


#pragma endregion

/* Once we reach the end of the growing state we can start closing points*/
#pragma region Closing Points

		if (growingState == false)
		{
#pragma region Close the Initial Point
			/* If the initial point is not closed close it as soon as we leave the growing state */
			if (Surface._pointStatComplete[zeroPoint] == false)
			{

				/* Add all the values of the LDList as
				at this point this is the value of the initial position */
				for (size_t i = 0; i < LDList.size(); i++)
				{
					for (size_t k = 0; k < 2; k++) { LDTot[k] += LDList[i][k]; }
				}
				Surface.closePoint(zeroPoint);
			}
#pragma endregion
	/* If the initial point is closed then we have to close the firt point in each list */
#pragma region Close the Points
			else
			{
#pragma region Forward Closing
				/* Close points only if there are opened Points */
				if (openPoints_f.size() >= 1)
				{

					/* Get the first point */
					key = openPoints_f[0];

					bufferCounter = 0;
					totLength = 0.;

					/* Add the values of the LDList from backwards
					until we reach 2*tau */
					while (totLength < (2 * tau))
					{
						LDTot[0] += LDList[LDList.size() - (1 + bufferCounter)][0];
						LDTot[1] += LDList[LDList.size() - (1 + bufferCounter)][1];

						bufferCounter += 1;					// Increase the buffer to select the correct value of the list
						totLength += timeStep;				// Increase the length of the Point
					}

					/* Save the Value of the point */
					Surface._LDallPoints[key].push_back(LDTot);

					/* Close the point */
					Surface.closePoint(key);
					openPoints_f.erase(openPoints_f.begin());	//Remove the point from the list

				}
#pragma endregion

#pragma region Backward Closing
				/* Close points only if there are opened Points */
				if (openPoints_b.size() >= 1)
				{

					/* Get the first point */
					key = openPoints_b[0];

					bufferCounter = 0;
					totLength = 0.;

					/* Add the values of the LDList from 0 to end
					until we reach 2*tau */
					while (totLength < (2 * tau))
					{
						LDTot[0] += LDList[bufferCounter][0];
						LDTot[1] += LDList[bufferCounter][1];

						bufferCounter += 1;					// Increase the buffer to select the correct value of the list
						totLength += timeStep;				// Increase the length of the Point
					}

					/* Save the Value of the point */
					Surface._LDallPoints[key].push_back(LDTot);

					/* Close the point */
					Surface.closePoint(key);
					openPoints_b.erase(openPoints_b.begin());	//Remove the point from the list

				}


#pragma endregion

			}
#pragma endregion
		}

#pragma endregion
	}

/* Once the trj has finished Save the values */
#pragma region Saving
	outputFile.open(calc + "_" + std::to_string(I) + "_" + "LD.txt", std::ios::out | std::ios::trunc);
	Surface.SavePointAver(outputFile);											//Calc the average and print it
	outputFile.close();

#pragma endregion

	return 0;
}
