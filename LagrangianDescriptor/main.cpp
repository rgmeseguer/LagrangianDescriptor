#include <cstdio>
#include <vector>				// Vector
#include <sstream>				// stringstream
#include <iomanip>				// setprecision

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
	if ((argc != 5))
	{
		std::cout << "Program Usage" << std::endl;
		std::cout << "LagD_Trj.out [tau] [bathMass] [TrjTime > (2xtau+1)] [Trj Number]" << std::endl;
		return 1;
	}


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
#pragma region System Initial Conditions and Variables
	srand(time(NULL));																			// Initialize the random seed

	/* Initiate the Oscillator */
	Oscillator Oscf({ 321.904484,-995.713452,1118.689573,-537.856726,92.976121,1.0,1.0,0.01 },	/* Coefficients from the oscillator */\
	{ 1., bathMass },																/* Mass Values */\
		DOS_V, DOS_G);																/* Potential and Gradient Functions */
	Oscillator Oscb({ 321.904484,-995.713452,1118.689573,-537.856726,92.976121,1.0,1.0,0.01 },	/* Coefficients from the oscillator */\
	{ 1., bathMass },																/* Mass Values */\
		DOS_V, DOS_G);																/* Potential and Gradient Functions */

	std::vector<double> R0 = { 1.36561 ,2.161769 };					// Initial Position
	double Energy = 3.691966889;									// Energy of the system
	Oscf.setInitP(R0, true);										// Initiate at the shadle point
	Oscb.setInitialVel_NVE(Energy, I);								// Set initial Random Velocities to keep energy
	Oscf.setInitP(R0, true);										// Initiate at the shadle point
	Oscb.setInitialVel_NVE(Energy, I);								// Set initial Random Velocities to keep energy

	/* Initialize the dynamics */
	double timeStep = 1.e-3;										// Set the Time Step/Precision of the Dynamic
	Dynamics Dynf;													// Set the forward Dynamic							
	Dynf.setTimeStep(timeStep);										// Set the Time step
	Dynf.setTime(TrjTime);											// and the total time (nsteps = totalTime/timeTtep)
	Dynamics Dynb;													// Set the backward Dynamic
	Dynb.setTimeStep(-1 * timeStep);									// Set the Time step
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
	LagDesc LD;													// Lagrangian Descriptor variable

	/* Create the surface of the where the LD values will be saved */
	LDSurfaceCreator Surface;
	std::vector<std::vector<std::string>> openPoints_f;			// Tracker of the points at which the LD forward is beig calculated	
	std::vector<std::vector<std::string>> openPoints_b;			// Tracker of the points at which the LD backward is beig calculated	
	std::vector<std::string> zeroPoint(Oscf._size * 2, " ");	// Initial point

	/* Variables to create the Keys of the surface's points */
	//std::vector<std::string> stringCoord(Oscf._size, " "), stringMoment(Oscf._size, " ");		// Strings needed for the key of the point
	std::vector<std::string> keyf(Oscf._size * 2, " ");			// Key of the position for the forward Trj 
	std::vector<std::string> keyb(Oscf._size * 2, " ");			// Key of the position for the backward Trj
	std::vector<std::string> key(Oscf._size * 2, " ");			// Temporary key


	/* Surface Variables Required to Udate the Points */
	//std::vector<std::vector<std::string>> OpenPoints;			// Tracker of the points at which the LD is beig calculated							
	std::vector<std::vector<double>> LDList;					// Buffer of the LD value at each step
	std::vector<double> LDTot = { 0.,0. };						// Buffer of the LD value
	//std::vector<double> LDTot_b = { 0.,0. };					// Buffer of the LD value
	double totLength = 0.;										// Lenght of the Buffers
	std::vector<double> tmp = { 0.,0. };						// Temporary Buffer
	int bufferCounter = 0;										// Counter to follow the values we need to remove from the buffer
#pragma endregion


	/* Open the inital Point */
	key = Surface.keyWrite(Oscf._position, Oscf.calcMomenta());
	Surface.addPoint(key);										// Create the new point
	zeroPoint = key;											// Save the Initial point


	for (size_t j = 0; j < Dynf._numberStep; j++)
	{
		if (growingState && j*Dynf._timeStep >= tau) { growingState = false; } // Once we reach time tau points begin to close

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
		/* Foward LD Calculation */
		tmp = LD.MAction(Oscf.calcMomenta(), Oscf._velocity);
		tmp[0] *= timeStep;
		tmp[1] *= timeStep;

		LDList.push_back(tmp);											// Save it at the end

		/* Backward LD Calculation */
		tmp = LD.MAction(Oscb.calcMomenta(), Oscb._velocity);
		tmp[0] *= timeStep;
		tmp[1] *= timeStep;

		LDList.insert(LDList.begin(), tmp);								// Save it at the begining

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
				//Surface._LDallPoints[keyf].push_back(LDTot_f);		// Add the value to the end of the values it has (if any)
				//Surface._pointCompletion[keyf] = totLength_f;			// Sets the completion status of the point
				openPoints_f.push_back(keyf);							// Include the point in the list of open points
			}
			/* And if exist but is not open, open the point again */
			else if (Surface._pointStatComplete[keyf] == true)
			{
				Surface.openPoint(keyf);								// Open the point
				//Surface._LDallPoints[keyf].push_back(LDTot_f);		// Add the value to the end of the values it has (if any)
				//Surface._pointCompletion[keyf] = totLength_f;			// Sets the completion status of the point
				openPoints_f.push_back(keyf);							// Include the point in the list of open points
			}

#pragma endregion
#pragma region Backward Opening Points

			/* If the key does not exist create a new one */
			if (Surface.doesKeyExist(keyb) == false)
			{
				Surface.addPoint(keyb);									// Create the new point
				Surface.openPoint(keyb);								// Open the point
				//Surface._LDallPoints[keyb].push_back(LDTot_b);		// Add the value to the end of the values it has (if any)
				//Surface._pointCompletion[keyb] = totLength_f;			// Sets the completion status of the point
				openPoints_b.push_back(keyb);							// Include the point in the list of open points

			}
			/* And if exist but is not open, open the point again */
			else if (Surface._pointStatComplete[keyb] == true)
			{
				Surface.openPoint(keyb);								// Open the point
				//Surface._LDallPoints[keyb].push_back(LDTot_b);		// Add the value to the end of the values it has (if any)
				//Surface._pointCompletion[keyb] = totLength_f;			// Sets the completion status of the point
				openPoints_b.push_back(keyb);							// Include the point in the list of open points
			}

#pragma endregion
		}

#pragma endregion

		/* Once we reach the end of the growing state we can start closing points*/
#pragma region Finish points

		if (growingState == false)
		{
#pragma region Close the Initial Point
			/* If the initial point is not closed close it as soon as we leave the growing state */
			if (Surface._pointCompletion[zeroPoint] == false)
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
			else
			{
#pragma region Forward Closing

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

					bufferCounter += 1;
					totLength += timeStep;
				}

				/* Close the point */


#pragma endregion
			}
		}


#pragma endregion




#if 0

		/* Update the forward and backard buffer variables that save the total value
	 of the LD and will be updated as the trjs advances */
#pragma region Update Buffer

	 /* If we reach tau from one side we have reached it from the other and
	 LDTot does not have to grow more but it has to Update its value */
		if (growingState)
		{
			//TODO: change this to be not dependent of the size.
			for (size_t i = 0; i < 2; i++)
			{
				bufferCounter += 1;											// Increase the counter to remember the position you have to remove from each buffer

				LDTot_f[i] += LDList.back()[i];								// Add the last value to the forward Dyn		
				LDTot_f[i] -= LDList[bufferCounter][i];						// and remove the first value added to the buffer

				LDTot_b[i] += LDList[0][i];									// Add the first value to the backward Dyn					
				LDTot_b[i] -= LDList[LDList.size - (1 + bufferCounter)][i];	// and remove the last value added the buffer			
			}

		}
		/* If not the LDTot is in growing state and needs to be Increased */
		else
		{
			for (size_t i = 0; i < 2; i++)
			{
				LDTot_f[i] += LDList.back()[i];				// Add the back value to the forward Dyn		
				LDTot_f[i] += LDList[0][i];					// and the first value
				totLength_f += 2 * timeStep;				// Increase the length by the two values added						

				LDTot_b[i] += LDList.back()[i];				// Add the back value to the backward Dyn			
				LDTot_b[i] += LDList[0][i];					// and the first value					
				totLength_b += 2 * timeStep;				// Increase the length by the two values added						

			}

		}

#pragma endregion

		/* Update points while in growing state */
#pragma region Update Points
		if (growingState)
		{
			/* Forward Updating */
			for (size_t k = 0; k < openPoints_f.size(); k++)			//Actualize every opened point
			{
				key = openPoints_f[k];								//Get the Key of an open point
				Surface._pointCompletion[key] += 2 * timeStep;		//Increase the length of the point

				/* Add the new LD value */
				for (size_t i = 0; i < 2; i++)
				{
					Surface._LDallPoints[key].back()[i] += LDList.back()[i];	// Add the last one 
					Surface._LDallPoints[key].back()[i] += LDList[0][i];		// and the first one 
				}

			}

			/* Backward Updating */
			for (size_t k = 0; k < openPoints_b.size(); k++)			//Actualize every opened point
			{
				key = openPoints_b[k];								//Get the Key of an open point
				Surface._pointCompletion[key] += 2 * timeStep;		//Increase the length of the point

				/* Add the new LD value */
				for (size_t i = 0; i < 2; i++)
				{
					Surface._LDallPoints[key].back()[i] += LDList.back()[i];	// Add the last one 
					Surface._LDallPoints[key].back()[i] += LDList[0][i];		// and the first one 
				}

			}

		}
		else
		{

		}

#pragma endregion


#pragma region Close Points

		/* If the two oldest open point is finished then close it */
		if (Surface._pointCompletion[openPoints_f[0]] >= (2 * tau))
		{
			Surface.closePoint(openPoints_f[0]);
			openPoints_f.erase(openPoints_f.begin());
		}

#pragma endregion


#endif // 0




	}


	//TODO:Solve the problem of saving the value of each point in both directions

	return 0;
}