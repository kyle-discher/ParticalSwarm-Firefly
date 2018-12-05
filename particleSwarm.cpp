#include "particleSwarm.h"
#include <windows.h>
/** \file particleSwarm.cpp
\brief Defines functions for the ParticleSwarm class
*/

void ParticleSwarm::allocateVelocity(double**& v) {
	mt19937_64 rng(time(NULL));
	uniform_real_distribution<double> dist(0, (this->bounds.upper - this->bounds.lower / 2));
	v = new double*[this->ns];
	for (int i = 0; i < this->ns; i++) {
		v[i] = new double[this->dim];
	}
	for (int i = 0; i < this->ns; i++) {
		for (int j = 0; j < this->dim; j++) {
			v[i][j] = dist(rng);
		}
	}
}

void ParticleSwarm::particleRow(ThreadPointers tp) {
	unique_lock<mutex> iterationLock(this->lock);
	uniform_real_distribution<double> dist(0, 1);
	for (int t = 0; t < this->iterations; t++) {
		for (int j = 0; j < this->dim; j++) {
			//Particle Velocity
			tp.v[tp.row][j] = tp.v[tp.row][j] + this->c1 * dist(*tp.rng) * (tp.pBest->population[tp.row][j] - tp.p->population[tp.row][j]) +
				this->c2 * dist(*tp.rng) * (tp.gBest[j] - tp.p->population[tp.row][j]);
			//Update Particles in row
			tp.p->population[tp.row][j] = tp.p->population[tp.row][j] + tp.v[tp.row][j];
			if (tp.p->population[tp.row][j] > this->bounds.upper) {
				tp.p->population[tp.row][j] = this->bounds.upper;
			}
			else if(tp.p->population[tp.row][j] < this->bounds.lower){
				tp.p->population[tp.row][j] = this->bounds.lower;
			}
		}
		//Check fitness for row
		tp.p->fitness[tp.row] = (benchmarkFunctions[tp.p->fitnessFunc])(tp.p->population[tp.row], tp.p->dim);
		if (tp.p->fitness[tp.row] < tp.pBest->fitness[tp.row]) {
			for (int i = 0; i < tp.p->dim; i++) {
				tp.pBest->population[tp.row][i] = tp.p->population[tp.row][i];
			}
			tp.pBest->fitness[tp.row] = tp.p->fitness[tp.row];
		}

		this->activeThreads--;
		if (activeThreads > 0) {
			this->cv.wait(iterationLock);
		}
		else {
			int index = tp.pBest->getBestSolutionIndex();
			if (tp.pBest->fitness[index] < *tp.gBestFitness) {
				for (int i = 0; i < tp.p->dim; i++) {
					tp.gBest[dim] = tp.pBest->population[index][i];
				}
				*tp.gBestFitness = tp.pBest->fitness[index];
			}
			this->activeThreads = this->ns;
			this->cv.notify_all();
		}
	}
}

double ParticleSwarm::runParticleSwarm() {
	mt19937_64 rng(time(NULL));
	Population p(this->ns, this->dim, this->fitnessFunc);
	p.generatePopulation(this->bounds.lower, this->bounds.upper);
	Population pBest(this->ns, this->dim, this->fitnessFunc);
	pBest.copyPopulation(&p);
	//Random Velocity half the distance between bounds
	double ** v = nullptr;
	this->allocateVelocity(v);
	double * gBest = new double[this->dim + 1];
	double gBestFitness;

	int bestIndex = p.getBestSolutionIndex();
	for (int i = 0; i < dim; i++) {
		gBest[i] = p.population[bestIndex][i];
	}
	gBestFitness = p.fitness[bestIndex];
	
	
	vector<thread> threads;

	for (int i = 0; i < this->ns; i++) {
		//Prep Pointers for Threads
		ThreadPointers tp;
		tp.p = &p;
		tp.row = i;
		tp.pBest = &pBest;
		tp.gBest = gBest;
		tp.gBestFitness = &gBestFitness;
		tp.v = v;
		tp.iterations = this->iterations;
		tp.rng = &rng;
		threads.push_back(thread(&ParticleSwarm::particleRow, this, tp));
	}
	for (int i = 0; i < this->ns; i++) {
		threads.at(i).join();
	}
	/*p.printPopulation();
	cout << endl;
	for (int i = 0; i < this->dim; i++) {
		cout << gBest[i] << " ";
	}
	cout << "Fitness: " << gBestFitness << endl;*/
	//Deallocate
	for (int i = 0; i < this->ns; i++) {
		delete v[i];
	}
	delete v;
	delete gBest;
	p.deletePopulation();
	pBest.deletePopulation();
	return gBestFitness;
}

ParticleSwarm::ParticleSwarm(int fitnessFunc, int ns, int dim, double upper, double lower, int iterations, double c1, double c2) {
	this->fitnessFunc = fitnessFunc;
	this->ns = ns;
	this->dim = dim;
	this->bounds.upper = upper;
	this->bounds.lower = lower;
	this->iterations = iterations;
	this->c1 = c1;
	this->c2 = c2;
	this->activeThreads = this->ns;
}

ParticleSwarm::ParticleSwarm(string param) {
	fstream fin;
	fin.open(param, ios::in);
	if (fin.is_open()) {
		fin >> this->fitnessFunc;
		fin >> this->ns;
		fin >> this->dim;
		fin >> this->bounds.upper;
		fin >> this->bounds.lower;
		fin >> this->iterations;
		fin >> this->c1;
		fin >> this->c2;
		this->activeThreads = this->ns;
	}
	else {
		cout << "Cannot Open Param" << endl;
	}
}

ParticleSwarm::~ParticleSwarm() {

}