/* This class creates the dynamics object that allows us to 
 calculate the dynamic process and its properties during
 the trajectories */

#pragma once
#include <iostream>					// cout	
#include <sstream>					// stringstream
#include <iomanip>					// setprecision
#include <fstream>					// ofstream
#include <vector>					// vector
#include "utilities.h"				// absol
#include "Oscillator.h"				// Oscillator

class Dynamics
{
	//double _time;							//Total Time of the Dynamic
	utilities ut;							//Utilities class for different functions

	//Langevin Variables
	double A, B;							//Temporary variables used for Langevin calculation
	double sigma;							//

	double _previousValue;					// Tool variable to keep track of one of the DoF to check if it crosses zero

public:
	int _numberStep;						//Total number of steps
	double _timeStep;						//Time step of the Dynamic
	double _halfTimeStep;					//Half timestep for the leapfrog algorithm
	
	//Langevin Variables
	std::vector<double> _gasDeviation;
	double _kineticBath;

	Dynamics();
	
	void setTimeStep(double);				// Sets the time step of the Dynamic 
	void setTime(double t);					// Set the total time of the Dynamic 
	void setLangevin(double mu, double beta, double gamma); // Sets the Langeving parameters used for the stochatic Dynamics
	
	bool doesItCross(double value, double crossing);

#if KEY_DETERM

	//Deterministic Dynamic Step
	void DynamicStep(Oscillator&);

#else
	//Langevin Dynamic Production Time Step
	void DynamicStep(Oscillator & osc, /*double B, double A,*/ double rtherm);

#endif // DETERM

};

Dynamics::Dynamics() { _previousValue = 0; }

/* Functions to initialize the Dynamics */
#pragma region InitializeDynamics

/// <summary>
/// Sets the time step of the Dynamic 
/// </summary>
void Dynamics::setTimeStep(double timeStep)
{
	_timeStep = timeStep;
	_halfTimeStep = timeStep / 2.;
}

/// <summary>
/// Set the total time of the Dynamic 
/// </summary>
void Dynamics::setTime(double t)
{
	_numberStep = ut.absol(int(t / _timeStep));
}


/// <summary>
/// Sets the Langeving parameters used for the stochatic Dynamics
/// </summary>
void Dynamics::setLangevin(double mu /*bath mass*/, double beta /*Effective Temperature*/, double gamma /*Dissipation Factor*/)
{
	sigma = sqrt(2 * gamma*fabs(_timeStep) / beta);
	A = gamma * ut.absol(_timeStep) / (2 * mu);
	B = 1 / (1 + A);
	A = (1 - A) / (1 + A);
	_kineticBath = 1 / (2 * beta);
}

#pragma endregion


/* Perfom a step in the Dynamic with two possible versions
 Deterministic or Stochastich */
#pragma region DynamicStep

#if KEY_DETERM

 /// <summary>
 /// Performs a Deterministic Dynamic Step using the 
 /// Leapfrog algorithm
 /// </summary>
void Dynamics::DynamicStep(Oscillator & osc)
{
	for (int unsigned i = 0; i < osc._size; i++)
	{
		osc._position[i] += (osc._velocity[i] += osc._acceleration[i] * _halfTimeStep) * _timeStep;
	}
	osc.calcAcceleration();
	for (int unsigned i = 0; i < osc._size; i++)
	{
		osc._velocity[i] += osc._acceleration[i] * _halfTimeStep;
	}
}

#else

 /// <summary>
 /// Performs a single molecular dynamics time step via the algorithm of 
 ///	Grønbech-Jensen & Farago
 /// </summary>
void Dynamics::DynamicStep(Oscillator & osc, double rtherm)
{

	//initial acceleration
	std::vector<double> initialAccel = osc._acceleration;

	//Set the new position
	for (int unsigned i = 0; i < osc._size - 1; i++)
	{
		osc._position[i] += _timeStep * (osc._velocity[i] + osc._acceleration[i] * _halfTimeStep);
	}

	//Bath has to be always the last defined degree of freedom
	osc._acceleration.back() += B * _timeStep*(osc._velocity.back() + osc._acceleration.back() * _halfTimeStep + rtherm / (2 * osc._redMass.back()));
	osc.calcAcceleration();

	//Set the new Velocity
	for (int unsigned i = 0; i < osc._size - 1; i++)
	{
		osc._velocity[i] += _halfTimeStep * (initialAccel[i] + osc._acceleration[i]);
	}
	osc._velocity.back() = A * osc._velocity.back() + _halfTimeStep * (initialAccel.back() + osc._acceleration.back()) + B * rtherm / osc._redMass.back();
}


#endif

#pragma endregion

/* Check if the followed property has crossed its zero point */
bool Dynamics::doesItCross(double value, double crossing)
{
	double value0 = _previousValue - crossing;
	double value1 = value - crossing;
	
	if (value0*value1 < 0)				// If the value changes sign we crossed the crossing point
	{
		_previousValue = value;
		return true;
	}
	else								// If not we didnt or the previous value is zero so we are at the first calling of the fucntion
	{
		_previousValue = value;
		return false;
	}
}
