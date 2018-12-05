/** \file particleSwarm.h
	\brief Defines ParticleSwarm and Thread Pointers
*/

#pragma once
#include "populationV2.h"
#include <fstream>
#include <thread>
#include <condition_variable>

/** \struct ThreadPointers
Defines data sent to threads in ParticleSwarm::runParticleSwarm()
*/
struct ThreadPointers {
	Population * p; ///< Pointer to the current population
	Population * pBest; ///< Pointer to personalBest
	double** v; ///< Particle Velocity
	double* gBest; ///< Global Best
	double* gBestFitness; ///< Global best Fitness
	int iterations; ///< Number of Iterations
	int row; ///< Row thread will be working on
	mt19937_64 * rng; ///< Same number generator for rand functions
};

/** \class ParticleSwarm
	The particle swarm class handles the data and threads for the Particle Swarm Algorithm
*/
class ParticleSwarm {
private:
	int ns; ///< Population Size
	int dim; ///< Dimensions
	int fitnessFunc; ///< Fitness Function
	bound bounds; ///< Fitnes Function Bounds
	int iterations; ///< Number of Iterations
	double c1, c2; ///< Personal Best, Global Best Modifiers
	mutex lock; ///< Thread Lock
	condition_variable cv; ///< Locks all threads to stay on same iteration and update Global best
	int activeThreads; ///< Counts threads still running
	/*! Update a single row of particles in the swarm.
	\param tp ThreadPointers
	\param this ParticleSwarm
	*/
	void particleRow(ThreadPointers tp);
	/*! Create Velocity vector
	\param v double**&
	*/
	void allocateVelocity(double**& v);
public:
	/*! Run an instance of ParticleSwarm returns best fitness
	\return double
	*/
	double runParticleSwarm();

	/*! Initialize Class in code
	\param fitnessFunc int
	\param ns int
	\param dim int
	\param upper double
	\param lower double
	\param iterations int
	\param c1 double
	\param c2 double
	*/
	ParticleSwarm::ParticleSwarm(int fitnessFunc, int ns, int dim, double upper, double lower, int iterations, double c1, double c2);
	/*! Initialize Class using file param
	\param param string
	*/
	ParticleSwarm(string param);
	/*! Default Destructor
	*/
	~ParticleSwarm();
};