/** \file population.h
	\brief Definition of the Population Struct

	Defines the functions and variables for the population class. 
	Also defines the bounds struct.
*/

#pragma once
#include <string>
#include <random>
#include <time.h>
#include <iostream>
#include "benchmarksV2.h"

using namespace std;


/** \struct bound
	Defines bounds for use in fitness function calculations.
*/
struct bound {
	double upper;
	double lower;
};

/** \struct Population
	\Brief Defines a population populationSize by dim.

	Creates a population to be used in Differential Evolution and Genetic Algorithm.
	The population has populationSize members each with dim values in them.
*/
struct Population {
	int populationSize; ///< Population Size
	int dim; ///< Dimensions
	int fitnessFunc; ///< Which fitness function will be used.
	double * fitness; ///< Fitnesses of the members in the population
	double** population = nullptr; ///< Population
	/*! Create a population with bounds min and max
		\param min double
		\param max double
		\return void
	*/
	void generatePopulation(double min, double max);
	int getBestSolutionIndex(); ///< Returns index of the best particle
	void emptyPopulation(); ///< Create a new blank population
	void updatePopulationFitness(); ///< Update the population fitness
	void printPopulation(); ///< Output the population to the console.
	void deletePopulation(); ///< Delete the population
	void sortByFitnessDecending(); ///< Sort Highest Fitness First
	/*! The constructor for Population. Sets the population size, dim and fitness function.
		\param populationSize int
		\param dim int
		\param fitnessFunc int
	*/
	void copyPopulation(Population* pop);
	Population(int populationSize, int dim, int fitnessFunc); ///< Create the struct with population size, dim and f
};