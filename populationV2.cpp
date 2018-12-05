/** \file population.cpp
\brief Defines functions for the Population struct
*/
#include "populationV2.h"
#include <algorithm>

void Population::generatePopulation(double min, double max) {
	mt19937_64 rng(time(NULL));
	uniform_real_distribution<double> dist(min, max);
	this->population = new double*[this->populationSize];
	for (int i = 0; i < this->populationSize; i++) {
		this->population[i] = new double[this->dim];
	}
	for (int i = 0; i < this->populationSize; i++) {
		for (int j = 0; j < this->dim; j++) {
			this->population[i][j] = dist(rng);
		}
	}

	this->fitness = new double[this->populationSize];
	this->updatePopulationFitness();
}

void Population::emptyPopulation() {
	this->population = new double*[this->populationSize];
	for (int i = 0; i < this->populationSize; i++) {
		this->population[i] = new double[this->dim];
	}
	this->fitness = new double[this->populationSize];
}
void Population::copyPopulation(Population* pop) {
	if (this->dim == pop->dim && this->populationSize == pop->populationSize && this->fitnessFunc == pop->fitnessFunc) {
		if (this->population != nullptr) {
			this->deletePopulation();
		}
		this->emptyPopulation();
		for (int i = 0; i < this->populationSize; i++) {
			for (int j = 0; j < this->dim; j++) {
				this->population[i][j] = pop->population[i][j];
			}
			this->fitness[i] = pop->fitness[i];
		}
	}
	else {
		cout << "Cannot Copy Population" << endl;
	}
}
int Population::getBestSolutionIndex() {
	int best = 0;
	for (int i = 1; i < this->populationSize; i++) {
		if (this->fitness[i] < this->fitness[best]) {
			best = i;
		}
	}
	return best;
}
void Population::updatePopulationFitness() {
	for (int i = 0; i < this->populationSize; i++) {
		this->fitness[i] = (benchmarkFunctions[this->fitnessFunc])(this->population[i], this->dim);
	}
}

void Population::sortByFitnessDecending() {
	for (int i = 0; i < this->populationSize - 1; i++) {
		int index = i;
		for (int j = i + 1; j < this->populationSize; j++) {
			if (this->fitness[j] > this->fitness[index]) {
				index = j;
			}
		}
		if (index != i) {
			swap(this->population[i], this->population[index]);
			swap(this->fitness[i], this->fitness[index]);
		}
	}
}

void Population::printPopulation() {
	cout << endl;
	for (int i = 0; i < this->populationSize; i++) {
		for (int j = 0; j < this->dim; j++) {
			cout << this->population[i][j] << " ";
		}
		cout << "Fitness: " << this->fitness[i] << endl;
	}
}
Population::Population(int populationSize, int dim, int fitnessFunc) {
	this->populationSize = populationSize;
	this->dim = dim;
	this->fitnessFunc = fitnessFunc;
}
void Population::deletePopulation() {
	for (int i = 0; i < this->populationSize; i++) {
		delete this->population[i];
	}
	delete this->population;
	delete this->fitness;
}