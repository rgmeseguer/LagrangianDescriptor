#pragma once
#include <vector>					//vector	
#include <math.h>					//sqrt

#include "utilities.h"				//powerd, absol

class LagDesc
{
	utilities ut;
public:
	LagDesc();

/* Lag Descriptor functions */
#if KEY_M0
	double MArc(std::vector<double> p);						
	double M5(std::vector<double> v, std::vector<double> a);
#endif

#if KEY_MF
	std::vector<double> MAction(std::vector<double>, std::vector<double>);
#endif 

};

LagDesc::LagDesc(){}

#if KEY_MF

///<surface>
/// Lagrangian Descriptor that calculates the action (pdq)
///</surface>
std::vector<double> LagDesc::MAction(std::vector<double> momenta, std::vector<double> rdot)
{
	std::vector<double> M(rdot.size());			// Variable to store the value of the LD

	for (size_t i = 0; i < rdot.size(); i++)
	{
		M[i] = sqrt(momenta[i] * rdot[i]);		// Calculate the action
	}
	return M;
}

#endif

#if KEY_M0
double LagDesc::MArc(std::vector<double> velocities) 
{
	double sqV = 0.;
	for (size_t i = 0; i < velocities.size(); i++)
	{
		sqV += (velocities[i] * velocities[i]);
	}
	return sqrt(sqV);
}

double LagDesc::M5(std::vector<double> velocity, std::vector<double> acceleration)
{
	/* Temporary Variables for storing calculations */
	double sqV = 0.; 
	double sqA = 0.;
	double sqB = 0.;

	/**/
	for (size_t i = 0; i < velocity.size(); i++)
	{
		sqV += (velocity[i] * velocity[i]);
		sqA += (acceleration[i] * acceleration[i]);
		sqB += (acceleration[i] * velocity[i]);
	}

	/**/
	return 1. / (1. + sqrt(ut.absol((sqV*sqA - sqB * sqB) / ut.powerd(sqV, 3))));
}

#endif
