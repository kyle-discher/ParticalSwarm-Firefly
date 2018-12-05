/** \file firefly.h
\brief Definition for the Firefly class

Firefly Algorithm behaves like fireflies do while they mate. 
This class runs the firefly algorithm
*/
#pragma once
#include "populationV2.h"

/** \class Firefly
	This class handles the execution of the firefly algorithm
*/
class Firefly {
private:
	int ns; ///< Solution Size
	int dim; ///< Number of Dimensions
	int fitnessFunc; ///< Fitness Function to be operated on
	bound bounds; ///< Function Bounds
	double beta0;
	double betaMin;
	double gamma;
	double alpha;
	int iterations; ///< Number of iterations to run the algorithm
public:
	/*! Run an instance of firefly and return the best fitness
	\return double
	*/
	double runFirefly();
	/*! Initialize Class in code
	\param fitnessFunc int
	\param ns int
	\param dim int
	\param upper double
	\param lower double
	\param iterations int
	\param alpha double
	\param beta0 double
	\param betaMin double
	\param gamma double
	*/
	Firefly(int fitnessFunc, int ns, int dim, double upper, double lower, int iterations, double alpha, double beta0, double betaMin, double gamma);
	/*! Initialize Class using file param
	\param param string
	*/
	Firefly(string params);
	/*! Default Destructor
	*/
	~Firefly() {};
};