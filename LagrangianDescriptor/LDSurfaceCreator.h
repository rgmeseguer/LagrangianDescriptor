/* This Class is designed to store the methods and 
the data used to create the points that will
define out LD surface. */

#pragma once
#include <fstream>				// std::ofstream, print
#include <sstream>				// stringstream
#include <vector>				// std::vector
#include <string>				// string
#include <map>					// map
#include <iostream>				// std::cout, std::fixed
#include <iomanip>				 // std::setprecision

class LDSurfaceCreator
{
	std::vector<std::vector<std::string>> _location;		// List of locations created in the map
	
	bool selectiveSaving;									// Decide if we want to save every point or those that cross zero in the followed DoF
	std::vector<std::vector<std::string>> _savingLoc;		// List of saving locations


public:
	LDSurfaceCreator(bool);
	
	/* Variables for storing values */
	std::map<std::vector<std::string>, std::vector<std::vector<double>>> _LDallPoints;	// Vector of LD values of a point 
	std::map<std::vector<std::string>, std::vector<double>> _LDavPoints;				// Average value of the LD of a point
	std::map<std::vector<std::string>, bool> _pointStatComplete;						// Status of the point
	std::map<std::vector<std::string>, double> _pointCompletion;						// Completion of the point
	std::map<std::vector<std::string>, int> _pointOpenedTimes;							// Times a point has been opened

	/* Variables for storing the Surface */
	void savingPointAdd(std::vector<std::string>);										// Adds one point to the list of saving points

	/* Functions to create the points */
	void addPoint(std::vector<std::string>);														// Create a new point using its "location" as a key
	void openPoint(std::vector<std::string>);														// Opens an existing point
	void closePoint(std::vector<std::string>);														// Closes an opened point
	void SavePointAver(std::ofstream&, std::vector<std::vector<std::string>> savingSelection);		// Calculates the average of the points and saves them
	bool doesKeyExist(std::vector<std::string>);													// Check the existence of a key

	std::vector<std::string> keyWrite(std::vector<double>, std::vector<double>);		// Write the key in the correct format



};



LDSurfaceCreator::LDSurfaceCreator(bool selective) 
{ 
	_location = {}; 
	selectiveSaving = selective;
}

///<surface>
/// Adds one point to the list of saving points.
///</surface>
void LDSurfaceCreator::savingPointAdd(std::vector<std::string> point)
{
	_savingLoc.push_back(point);
}

///<surface>
/// Opens an existing point.
///</surface>
void LDSurfaceCreator::openPoint(std::vector<std::string>loc)
{
	_pointStatComplete[loc] = false;		//Set the status of the point as not completed
	_pointOpenedTimes[loc] += 1;			//Increase the count of the times this point has been opened
} 

///<surface>
/// Closes an existing point.
///</surface>
void LDSurfaceCreator::closePoint(std::vector<std::string>loc) 
{
	_pointStatComplete[loc] = true;			//Set the status of the point as completed
}					

///<surface>
/// Create a new point using its "location" as a key
///</surface>
void LDSurfaceCreator::addPoint(std::vector<std::string>location)
{
	_LDallPoints[location] = {};			//Initialize the point
	_pointCompletion[location] = 0;			//Length of the point begins at 0
	_pointStatComplete[location] = true;	//The point begins completed/closed
	_pointOpenedTimes[location] = 0;		//Number of times the point is opened
	_location.push_back(location);			//Add its key to the list of keys
}

///<surface>
/// Create a new point using its "location" as a key
///</surface>
bool LDSurfaceCreator::doesKeyExist(std::vector<std::string>key)
{
	for (size_t i = 0; i < _location.size(); i++)
	{
		if (key == _location[i]) { return true; }
	}
	return false;
}

///<surface>
///Saves the average of the point in a file
///</surface>
void LDSurfaceCreator::SavePointAver(std::ofstream &sfile, std::vector<std::vector<std::string>> savingSelection)
{
	/* If selectiveSaving is true only save those points that are marked to be saved*/
	if (selectiveSaving)
	{
		for (size_t i = 0; i < savingSelection.size(); i++)				//Run over all the keys created
		{
			std::vector<std::string> key = savingSelection[i];			//Get the key					

#pragma region Average Calculation

		/* Calculate the Average of the Point */
			_LDavPoints[key] = { 0.,0. };

			/* Add all the LD values of the Point */
			for (size_t j = 0; j < _LDallPoints[key].size(); j++)
			{
				_LDavPoints[key][0] += _LDallPoints[key][j][0];
				_LDavPoints[key][1] += _LDallPoints[key][j][1];
			}

			/* Divide it by the total number of values */
			_LDavPoints[key][0] /= _pointOpenedTimes[key];
			_LDavPoints[key][1] /= _pointOpenedTimes[key];

#pragma endregion
			/* Saves the point in the file */
			for (size_t j = 0; j < key.size(); j++)
			{
				sfile << key[j] << ' ';
			}
			sfile << _LDavPoints[key][0] << ' ';
			sfile << _LDavPoints[key][1] << ' ';
			sfile << _LDavPoints[key][0] + _LDavPoints[key][1] << ' ';
			sfile << std::endl;
		}

	}
	else
	{
		for (size_t i = 0; i < _location.size(); i++)				//Run over all the keys created
		{
			std::vector<std::string> key = _location[i];			//Get the key					

#pragma region Average Calculation

		/* Calculate the Average of the Point */
			_LDavPoints[key] = { 0.,0. };

			/* Add all the LD values of the Point */
			for (size_t j = 0; j < _LDallPoints[key].size(); j++)
			{
				_LDavPoints[key][0] += _LDallPoints[key][j][0];
				_LDavPoints[key][1] += _LDallPoints[key][j][1];
			}

			/* Divide it by the total number of values */
			_LDavPoints[key][0] /= _pointOpenedTimes[key];
			_LDavPoints[key][1] /= _pointOpenedTimes[key];

#pragma endregion
			/* Saves the point in the file */
			for (size_t j = 0; j < key.size(); j++)
			{
				sfile << key[j] << ' ';
			}
			sfile << _LDavPoints[key][0] << ' ';
			sfile << _LDavPoints[key][1] << ' ';
			sfile << _LDavPoints[key][0] + _LDavPoints[key][1] << ' ';
			sfile << std::endl;
		}

	}
}

///<surface>
///Writes the key in the correct format
///</surface>
std::vector<std::string> LDSurfaceCreator::keyWrite(std::vector<double> position, std::vector<double> momenta)
{
	int size = position.size();
	std::stringstream stream;													//Temporary string
	std::vector<std::string> stringCoord(size, " "), stringMoment(size, " ");	// Strings needed for the key of the point
	std::vector<std::string> key(size * 2, " ");								// The value to return

	for (size_t i = 0; i < size; i++)
	{

		/*Print the new key in correct string format*/
		stream << std::fixed << std::setprecision(3) << position[i];	// Position
		stringCoord[i] = stream.str(); stream.str(std::string());

		stream << std::fixed << std::setprecision(5) << momenta[i];		// Momenta
		stringMoment[i] = stream.str(); stream.str(std::string());

		key[i] = stringCoord[i];
		key[i + size] = stringMoment[i];									//The location of the new point
	}

	return key;
}
