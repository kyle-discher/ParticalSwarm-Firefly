/** \file benchmarksV2.cpp
\brief Function definitions for the function headers in benchmarksV2.h
*/

#include "benchmarksV2.h"
#include <iostream>
#include <fstream>

double schwefels(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim; i++) {
		retVal += sin(sqrt(abs(input[i]))) * -input[i];
	}
	return retVal;
}
double deJong1st(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim; i++) {
		retVal += input[i] * input[i];
	}
	return retVal;
}
double rosenbrock(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		retVal += 100 * pow((pow(input[i], 2)) - input[i+1], 2) + pow(1 - input[i], 2);
	}
	return retVal;
}
double rastrigin(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim; i++) {
		retVal += pow(input[i], 2) - 10 * cos(2 * M_PI * input[i]);
	}
	return 2 * dim * retVal;
}
double griewangk(double input[], int dim) {
	double sum = 0.0;
	double capitalPI = cos(input[0]);
	for (int i = 0; i < dim; i++) {
		sum += pow(input[i], 2) / 4000;
	}
	for (int j = 1; j < dim; j++) {
		capitalPI *= cos(input[j] / sqrt(j + 1));
	}
	return 1 + sum - capitalPI;
}
double sineEnvelopSine(double input[], int dim) {
	double retVal = 0.0;
	double redundant = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		redundant = pow(input[i], 2) + pow(input[i+1], 2);
		retVal += .5 + (pow(sin(redundant - .5), 2) / pow(1 + .001*redundant, 2));
	}
	return -1 * retVal;
}
double stretchVSine(double input[], int dim) {
	double retVal = 0.0;
	double rootVal = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		rootVal = pow(input[i], 2) + pow(input[i+1], 2);
		retVal += pow(rootVal, .25) * sin(50 * pow(rootVal, .1) + 1);
	}
	return retVal;
}
double ackleys1(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		retVal += (1 / pow(M_E, .2)) * sqrt(pow(input[i], 2) + pow(input[i + 1], 2)) + 3 * (cos(2 * input[i]) + sin(2 * input[i + 1]));
	}
	return retVal;
}
double ackleys2(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		retVal += 20 + M_E - 20 / (pow(M_E, .2 * sqrt((pow(input[i], 2) + pow(input[i + 1], 2)) / 2))) - pow(M_E, .5 * (cos(2 * M_PI * input[i]) + cos(2 * M_PI * input[i + 1])));
	}
	return retVal;
}
double eggHolder(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		retVal -= input[i] * sin(sqrt(abs(input[i] - input[i+1] - 47))) + (input[i+1] + 47) * sin(sqrt(abs(input[i+1] + 47 + input[i] / 2)));
	}
	return retVal;
}
double rana(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		retVal += input[i] * sin(sqrt(abs(input[i+1] - input[i] + 1))) * cos(sqrt(abs(input[i+1] + input[i] + 1)));
		retVal += (input[i+1] + 1) * cos(sqrt(abs(input[i+1] - input[i] + 1))) * sin(sqrt(abs(input[i+1] + input[i] + 1)));
	}
	return retVal;
}
double pathological(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		retVal += .5 + (pow(sin(sqrt(100 * pow(input[i], 2) + pow(input[i+1], 2))), 2) - .5) / (1 + .001 * pow(pow(input[i], 2) - 2 * input[i] * input[i+1] + pow(input[i+1], 2), 2));
	}
	return retVal;
}
double michalewicz(double input[], int dim) {
	double retVal = 0.0;
	for (int i = 0; i < dim; i++) {
		retVal += sin(input[i]) * sin(((i + 1) * pow(input[i], 2)) / M_PI);
	}
	return retVal * -1;
}
double mastersCosine(double input[], int dim) {
	double retVal = 0.0;
	double rootValue = 0.0;
	for (int i = 0; i < dim - 1; i++) {
		rootValue = pow(input[i], 2) + pow(input[i+1], 2) + .5 * input[i+1] * input[i];
		retVal += pow(M_E, -(1 / 8) * rootValue) * cos(4 * sqrt(rootValue));
	}
	return retVal;
}
double shekelsFoxholes(double input[], int dim) {
	double retVal = 0.0;
	double nSum = 0.0;
	for (int i = 0; i < 30; i++) {
		for (int j = 0; j < dim; j++) {
			nSum += pow(input[j] - FOXHOLEA[i][j], 2);
		}
		retVal += 1 / (FOXHOLEC[i] + nSum);
		nSum = 0.0;
	}
	return -1 * retVal;
}

double (*benchmarkFunctions[15])(double input[], int dim) = { schwefels, deJong1st, rosenbrock, rastrigin, griewangk, sineEnvelopSine,
stretchVSine, ackleys1, ackleys2, eggHolder, rana, pathological,
michalewicz, mastersCosine, shekelsFoxholes };