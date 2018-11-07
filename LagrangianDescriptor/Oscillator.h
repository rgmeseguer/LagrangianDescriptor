/* This class creates the oscillator that will be used in the Dynamics
code to run a dynamic. It requires an external function that defines 
the way to calculate the potential and gradient of the Oscillator */

#pragma once
#include <vector>					// vector
#include <iostream>					// std::cout, std::fixed
#include <math.h>					// sqrt fabs 	
#include <stdlib.h>					// abs srand rand
#include <time.h>					// time	

#include "utilities.h"				// RandomFloat

class Oscillator
{
	//utilities functions
	utilities ut; 

	//Necessary functions of the gradient and the potential function
	typedef double(*PotFunct)(std::vector<double>, std::vector<double>);
	PotFunct potentialFunction;
	typedef std::vector<double>(*GradFunct)(std::vector<double>, std::vector<double>);
	GradFunct gradientFunction;

public:
	std::vector<double> _coeff;											// coefficients of the system

	Oscillator(std::vector<double> coeff, std::vector<double> mass, \
		       double(PotFunc)(std::vector<double>, std::vector<double>), \
			   std::vector<double>(GradFunc)(std::vector<double>, std::vector<double>));
	std::vector<double> _redMass;										// Reduced Mass of the 2 oscillators
	size_t _size;														// Size of the system

	//Storage Variables
	std::vector<double> _position;										// Position
	std::vector<double> _velocity;										// velocity
	std::vector<double> _kinEnergy;										// Kinetic Energy
	std::vector<double> _gradient;										// Gradient
	std::vector<double> _acceleration;									// Acceleration
	double _potEnergy;													// Potential Energy
		
	/* Sets the initial conditions of the oscillator */
	void setInitP(std::vector<double>, bool);
	void setInitV(std::vector<double>, bool);

	/* Calculate the properties of the Oscillator */
	void calcAcceleration();											// Recalculate the Acceleration
	std::vector<double> calcMomenta();									// Calculate the momenta
	double calcTotalEnergy();											// Total energy
	void calcGradient();												// Recalculate the gradient

	/* Changing initial conditions to keep Energy constant */
	
	//Sets the velocity of the two motions to keep the energy constant I=[0:3999] is the index for the velocity of the bath
	void setInitialVel_NVE(double Energy, int I, int scale);

	//Finds the point for the bath knowing the position of the remaining DF and the KE being 0
	void findBathPos_NVE(double Energy);
	
	//Sets a random velocity for the two motions keeping the energy constant
	void randomInitialVel_NVE(double Energy, int Syst, int Bath);

};

Oscillator::Oscillator\
(std::vector<double> coeff, std::vector<double> mass,\
 double(PotFunc)(std::vector<double>, std::vector<double>),\
 std::vector<double>(GradFunc)(std::vector<double>, std::vector<double>))
{
	_coeff = coeff;
	_redMass = mass;
	_size = _redMass.size();
	std::vector<double> zerotmp(_size, 0.);
	_position = zerotmp;
	_velocity = zerotmp;
	_acceleration = zerotmp;
	_kinEnergy = zerotmp;
	_gradient = zerotmp;
	_potEnergy = 0.;
	
	potentialFunction = PotFunc;
	gradientFunction  = GradFunc;

}

/* Sets the initial conditions of the oscillator */
#pragma region SetInitCond

///<sumary>
///Set the position of the system and recalculate the energy if true
///</sumary>
void Oscillator::setInitP(std::vector<double> coord, bool Energy)
{
	_position = coord;
	_gradient = gradientFunction(_coeff, _position);//the gradient
	for (auto i = 0; i < _size; i++)//the acceleration
	{
		_acceleration[i] = -_gradient[i] / _redMass[i];
	}

	if (Energy)
	{
		_potEnergy = potentialFunction(_coeff, _position);//set the potential energy
	}

}

///<sumary>
/// Set the velocities of the system and recalculate the energy if true
///</sumary>
void Oscillator::setInitV(std::vector<double> veloc, bool Energy)
{
	_velocity = veloc;
	if (Energy)
	{
		//also set the kinetic energy
		for (auto i = 0; i < _size; i++)
		{
			_kinEnergy[i] = _redMass[i] / 2 * (_velocity[i] * _velocity[i]);
		}
	}
}

#pragma endregion

/* Calculate the properties of the Oscillator */
#pragma region CalcProp

///<sumary>
/// Recalculate the Acceleration
///</sumary>
void Oscillator::calcAcceleration()
{
	_gradient = gradientFunction(_coeff, _position);		// The gradient
	for (size_t i = 0; i < _size; i++)						// The acceleration
	{
		_acceleration[i] = -_gradient[i] / _redMass[i];
	}
}

///<sumary>
/// calculate and returns the momenta
///</sumary>
std::vector<double> Oscillator::calcMomenta()
{
	std::vector<double> mom(_size, 0.);
	for (auto i = 0; i < _size; i++)
	{
		mom[i] = _velocity[i] * _redMass[i];
	}
	return mom;
}

///<sumary>
/// calculate and returns the total Energy
///</sumary>
double Oscillator::calcTotalEnergy()
{
	double totalKinEnergy = 0.;
	_potEnergy = potentialFunction(_coeff, _position);					// Calculates the potential energy

	for (auto i = 0; i < _size; i++)
	{
		_kinEnergy[i] = _redMass[i] / 2 * (_velocity[i] * _velocity[i]); // and the kinetic energy
		totalKinEnergy += _kinEnergy[i];
	}

	return _potEnergy + totalKinEnergy;
}

///<sumary>
/// recalculate the gradient
///</sumary>
void Oscillator::calcGradient()
{
	// Recalculate the gradient
	_gradient = gradientFunction(_coeff, _position);
}

#pragma endregion

/* Changing initial conditions to keep Energy constant */
#pragma region KeepEnergy

///<sumary>
///Sets the velocity of the two motions to keep the energy constant I=[0:3999] is the index for the velocity of the bath
///</sumary>
void Oscillator::setInitialVel_NVE(double Energy, int I, int scale)
{
	int precission = ut.powerd(10, scale);

	
	//TODO: Make this not size dependent
#pragma region Select the direction of the system and the bath
	int sysdir, bathdir;
	if (I > 4 * precission)
	{
		std::cout << "Error:Index too high, Max number of calculations " << 4 * precission << std::endl;
		std::terminate();
	}
	else if (I > 3 * precission)
	{
		I = I - 3 * precission;
		sysdir = -1;
		bathdir = -1;
	}
	else if (I > 2 * precission)
	{
		I = I - 2 * precission;
		sysdir = 1;
		bathdir = -1;

	}
	else if (I > 1 * precission)
	{
		I = I - 1 * precission;
		sysdir = -1;
		bathdir = 1;
	}
	else
	{
		sysdir = 1;
		bathdir = 1;
	}

#pragma endregion

	double kineticTot = Energy - _potEnergy;				//Calculate the total kinetic energy
	double dek = kineticTot / precission;					//Divide the energy depending on the precission
	double kinSys = dek * I;								//Set the KE of the system depending on the index
	double kinBath = kineticTot - kinSys;					//Set the KE of the bath

	/*Set the velocities of the system */
	_velocity[0] = sqrt(2 * kinSys / _redMass[0]) * sysdir;
	_velocity[1] = sqrt(2 * kinBath / _redMass[1])* bathdir;

	/* Also set the kinetic energy */
	for (auto i = 0; i < _size; i++)
	{
		_kinEnergy[i] = _redMass[i] / 2 * (_velocity[i] * _velocity[i]);
	}

}

///<sumary>
///Finds the point for the bath knowing the position of the remaining DF and the KE being 0
///</sumary>
void Oscillator::findBathPos_NVE(double Energy)
{
	//We need to searh for the root, the coordinate in the bath that keeps
	//the energy constant with the given position of the system
	double postition0, position1;		//Inital and final Positions
	double stepSize;					//Step size
	double E0, E1;//Energies


#pragma region Initial Conditions
	//Set the initial Velocities to 0
	for (auto vel : _velocity) { vel = 0.; }
	setInitV(_velocity, true);
	//Set the initial 3 points to search for the root
	postition0 = 3.0;//Beggining
	position1 = 1.5;//End
	stepSize = (position1 - postition0) / 100.;//Divide the range

	//Energy on the Beggining
	std::vector<double> init = _position;
	init.back() = postition0;
	setInitP(init, true);
	E0 = calcTotalEnergy() - Energy; //Energy differece

	position1 = postition0;//Keep track of the Beggining Point
	E1 = E0;//Keep track of the Beggining EnergyDiff
#pragma endregion

	//Get close to the root
	while (fabs(E0) > 1.e-7)
	{
		//Do another step
		position1 += stepSize;
		init.back() = position1;
		setInitP(init, true);
		E1 = calcTotalEnergy() - Energy;
		//Did you miss the root?
		if (/*you cross it*/ ((E1*E0 < 0.)
			|| (fabs(E1) > fabs(E0)))
			&& /*and you are close enough to the root*/ (E0 < 10))
		{
			//You missed the root get back and reduce the step
			position1 -= stepSize; stepSize *= 0.5;
			init.back() = position1;
			setInitP(init, true);
			E1 = calcTotalEnergy() - Energy;
		}
		else //Reduce the step so it cannot end in a infinite loop
		{
			stepSize *= 0.999;
		}
		//If no Root do not calculate
		if (fabs(stepSize) < 1.e-10) { std::cout << "something went wrong " << position1 << " " << fabs(E0) << std::endl; break; }
		//Actualize the back position
		E0 = E1;
		postition0 = position1;
	}
	init.back() = postition0;
	setInitP(init, true);//Set the final position of the bath


}

///<sumary>
///Sets the velocity of the two motions to keep the energy constant I=[0:3999] is the index for the velocity of the bath
///</sumary>
void Oscillator::randomInitialVel_NVE(double Energy, int Syst, int Bath)
{
	srand(time(NULL));
	int one[2] = { -1,1 };
	double kineticTot = Energy - _potEnergy;//Calculate the total kinetic energy
	
											//randomly divide it between the system and the bath
	double kinBath = ut.RandomFloat(0., kineticTot);
	double kinSys = kineticTot - kinBath;

	int ranIndex = rand() % 2;
	_velocity[Syst] = sqrt(2 * kinSys / _redMass[Syst]) * one[ranIndex];

	ranIndex = rand() % 2;
	_velocity[Bath] = sqrt(2 * kinBath / _redMass[Bath])* one[ranIndex];
	setInitV(_velocity, true);

	std::cout << Energy << " " << calcTotalEnergy() << " " << _kinEnergy[0] << " " << _kinEnergy[1] << std::endl;
}

#pragma endregion
