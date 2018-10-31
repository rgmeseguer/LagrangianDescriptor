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

void programUsage()
{
	std::cout << "Program Usage" << std::endl;
	std::cout << "LagD_Trj.out [tau] [bathMass] [TrjTime > (2xtau+1)] [Trj Number] [Direction {-1,1}] [[Crossing Point]]" << std::endl;
}

int main(int argc, char *argv[])
{

/* Check the correct Parameters */
#pragma region Program Usage
	bool selecSaving;
	if ((argc == 6))
	{
		selecSaving = false;
	}
	else if ((argc != 7))
	{
		programUsage();
		return 1;
	}
	else
	{
		selecSaving = true;
	}
#pragma endregion

/* Variables introduced by the User */
#pragma region Input variables
	double tau = strtof(argv[1], NULL);								// Tau value
	double bathMass = strtof(argv[2], NULL);						// Mass of the Bath
	double TrjTime = strtof(argv[3], NULL);							// Time of the Trj
	int I = int(strtof(argv[4], NULL));								// Number of the trajectory
	int direction = int(strtof(argv[5], NULL));						// Direction of the time Step
	double crossPoint;												// Point at wich the crossing trj is saved
	if (selecSaving)
	{
		crossPoint = strtof(argv[6], NULL);					
	}
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
	if (tau < 1) { stream << std::fixed << std::setprecision(1) << tau; stau = stream.str(); }
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

	std::vector<double> R0 = { 1.4 ,2.18 };							// Initial Position
	double Energy = 3.691966889;									// Energy of the system
	double timeStep = 1.e-3 ;										// Set the Time Step/Precision of the Dynamic multiply by the direction

#pragma endregion

/* Initialize the Variables */
#pragma region Process Required Variables
	/* We are going to perform two Dynamics at the same time,
	one in Foward and the other in Backward time.
	To do so, we need to initialize 2 Oscillators and 2 Dynamics*/

	/* Initialize the Oscillators */
#pragma region Oscillator Variables
	Oscf.setInitP(R0, true);										// Initiate at the shadle point
	Oscf.setInitialVel_NVE(Energy, I, 4);							// Set initial Random Velocities to keep energy
	Oscb.setInitP(R0, true);										// Initiate at the shadle point
	Oscb.setInitialVel_NVE(Energy, I, 4);							// Set initial Random Velocities to keep energy

#pragma endregion

	/* Initialize the dynamics */
#pragma region Dynamics Variables
	Dynamics Dynf;													// Set the forward Dynamic							
	Dynf.setTimeStep(timeStep * direction);							// Set the Time step
	Dynf.setTime(TrjTime);											// and the total time (nsteps = totalTime/timeTtep)
	Dynamics Dynb;													// Set the backward Dynamic
	Dynb.setTimeStep(-1 * timeStep * direction);					// Set the Time step
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

#pragma endregion
	
	/* Create the surface of the where the LD values will be saved */
#pragma region Surface Variables
	LDSurfaceCreator Surface(selecSaving);							// Surface creator with selective saving
	std::vector<std::vector<std::string>> openPoints;				// Tracker of the points at which the LD forward is beig calculated	
	std::vector<std::string> previousPoint(Oscf._size * 2, " ");	// Variable to keep the two points where the foward trajectory crosses zero

	std::vector<std::string> zeroPoint(Oscf._size * 2, " ");		// Initial point

	/* Variables to create the Keys of the surface's points */
	std::vector<std::string> key(Oscf._size * 2, " ");				// Key of the position for the Trj 


	/* Surface Variables Required to Udate the Points */
	std::vector<std::vector<double>> LDList;						// Buffer of the LD value at each step
	std::vector<double> LDTot = { 0.,0. };							// Buffer of the LD value
	double totLength = 0.;											// Lenght of the Buffers
	int bufferCounter = 0;											// Counter to follow the values we need to remove from the buffer
#pragma endregion


	/* Initialize the Lagrangian Descriptor Object*/
	LagDesc LD;														// Lagrangian Descriptor variable

#pragma endregion

/* Open the inital Point */
#pragma region Open first Point
	key = Surface.keyWrite(Oscf._position, Oscf.calcMomenta());
	Surface.addPoint(key);											// Create the new point
	Surface.openPoint(key);											// Open the point
	zeroPoint = key;												// Save the Initial point
	previousPoint = key;											// Set the first previous point
	
	/* If we have activated the selective saving and
	 the initial condition has a zero momenta on the
	 bath then save that point too */
	if (Surface.selectiveSaving)
	{
		Dynf.doesItCross(Oscf.calcMomenta()[1], crossPoint);		// Call the function once to get the previous value saved
	}

#pragma endregion

/* Initate the trayectory */
	for (size_t j = 0; j < Dynf._numberStep; j++)
	{
		/* Once we reach time tau points begin to close opened Points and the Backwar Trj Stops */ 
		if (growingState && j*timeStep >= tau) { growingState = false; } 

/* Perform a Dinamic Step in forward and backward direction */
#pragma region Dynamic Step
#if KEY_DETERM==1
		Dynf.DynamicStep(Oscf);										//Dynamic forward step
		if (growingState == true) { Dynb.DynamicStep(Oscb); }		//Dynamic backward step
#else
		rtherm = Dynf.sigma*gasdev[j];								//Bath Random
		Dynf.DynamicStep(Oscf, rtherm);								//Dynamic foward step
		if (growingState == true)
		{
			rtherm = Dynb.sigma*gasdev[j];								//Bath Random
			Dynb.DynamicStep(Oscb, rtherm);								//Dynamic backward step
		}

#endif
#pragma endregion

/* With the new positions write the two new keys obtained */
#pragma region Key Write
		/* Forward Key Write*/
		key = Surface.keyWrite(Oscf._position, Oscf.calcMomenta());

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

		/* Only save the points in the backward direction during the growing state */
		if (growingState == true)
		{
			/* Backward LD Calculation */
			LDTot = LD.MAction(Oscb.calcMomenta(), Oscb._velocity);
			LDTot[0] *= timeStep;
			LDTot[1] *= timeStep;

			LDList.insert(LDList.begin(), LDTot);							// Save it at the begining
		}


#endif
#pragma endregion

/* At every step Open those points that require to do so */
#pragma region Open Points

		/* Open new point only until reach Time-tau otherwise we will have unfinished points */
		if (j*timeStep < Dynf._numberStep*timeStep - tau)
		{
			/* If the key does not exist create a new one */
			if (Surface.doesKeyExist(key) == false)
			{
				Surface.addPoint(key);								// Create the new point
				Surface.openPoint(key);								// Open the point
				openPoints.push_back(key);							// Include the point in the list of open points

				/* If we have activated the selective saving and
				the trajectory crossed the crossPoint in the bath momenta 
				from the previous step save those two points */
				if (Surface.selectiveSaving && Dynf.doesItCross(stof(key[3]), crossPoint))
				{
					//std::cout << previousPoint[1] << " " << previousPoint[3] << " " << key[1] << " " << key[3] << " " << j << std::endl;
					Surface.savingPointAdd(previousPoint);			// Because both points are toghether in the list we will be able to know between which two interpolate
					Surface.savingPointAdd(key);				

				} // we don't need to save this again because same point will always cross zero and will have it values saved.

			}
			/* And if exist but is not open, open the point again */
			else if (Surface._pointStatComplete[key] == true)
			{
				Surface.openPoint(key);								// Open the point
				openPoints.push_back(key);							// Include the point in the list of open points
			}
			
			previousPoint = key;									// Actualize the previous point
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
				
				/* Save the Value of the point */
				Surface._LDallPoints[zeroPoint].push_back(LDTot);
				Surface.closePoint(zeroPoint);
			}
#pragma endregion
	/* If the initial point is closed then we have to close the firt point in each list */
#pragma region Close the Points
			else
			{
				/* Close points only if there are opened Points */
				if (openPoints.size() >= 1)
				{
					/* Get the first point */
					key = openPoints[0];

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
					openPoints.erase(openPoints.begin());	//Remove the point from the list

				}

			}
#pragma endregion
		}

#pragma endregion
	}

/* Once the trj has finished Save the values */
#pragma region Saving
	//outputFile.open(calc + "_" + std::to_string(I) + "_" + "LD.txt", std::ios::out | std::ios::trunc);
	//Surface.SavePointAver(outputFile);											//Calc the average and print it

	outputFile.open(calc + "_" + std::to_string(I) + "_" + std::to_string(direction) + "_" + "LD.txt", std::ios::out | std::ios::trunc);
	Surface.SavePointAver(outputFile);
	outputFile.close();

#pragma endregion

	return 0;
}
