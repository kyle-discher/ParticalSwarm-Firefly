/** \file firefly.cpp
\brief Defines functions for the firefly class
*/
#include "firefly.h"
#include <fstream>

double Firefly::runFirefly() {
	mt19937_64 rng(time(NULL));
	uniform_real_distribution<double> dist(0, 1);
	Population pop(this->ns, this->dim, this->fitnessFunc);
	pop.generatePopulation(this->bounds.lower, this->bounds.upper);
	Population tempPop(this->ns, this->dim, this->fitnessFunc);
	for (int t = 0; t < this->iterations; t++) {
		tempPop.copyPopulation(&pop);
		for (int i = 0; i < this->ns; i++) {
			for (int j = i + 1; j < this->ns; j++) {
				double r = 0.0;
				for (int l = 0; l < this->dim; l++) {
					r += pow(tempPop.population[i][l] - tempPop.population[j][l], 2) ;
				}
				r = sqrt(r);
				if (tempPop.fitness[i] > tempPop.fitness[j]) {
					double beta = this->betaMin  + this->beta0 * exp((-this->gamma) * pow(r, 2));
					for (int z = 0; z < dim; z++) {
						double rand = dist(rng);
						double temp = this->alpha * (rand - .5);
						pop.population[i][z] = tempPop.population[i][z] + beta * (tempPop.population[j][z] - tempPop.population[i][z]) + temp;
					}
				}
			}
		}
		pop.updatePopulationFitness();
		pop.sortByFitnessDecending();
	}
	return pop.fitness[pop.getBestSolutionIndex()];
}

Firefly::Firefly(int fitnessFunc, int ns, int dim, double upper, double lower, int iterations, double alpha, double beta0, double betaMin, double gamma) {
	this->fitnessFunc = fitnessFunc;
	this->ns = ns;
	this->dim = dim;
	this->bounds.upper = upper;
	this->bounds.lower = lower;
	this->iterations = iterations;
	this->alpha = alpha;
	this->beta0 = beta0;
	this->betaMin = betaMin;
	this->gamma = gamma;
}

Firefly::Firefly(string params) {
	fstream fin;
	fin.open(params, ios::in);
	if (fin.is_open()) {
		fin >> this->fitnessFunc;
		fin >> this->ns;
		fin >> this->dim;
		fin >> this->bounds.upper;
		fin >> this->bounds.lower;
		fin >> this->iterations;
		fin >> this->alpha;
		fin >> this->beta0;
		fin >> this->betaMin;
		fin >> this->gamma;
	}
	else {
		cout << "Cannot Open Param" << endl;
	}
}