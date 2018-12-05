#include "particleSwarm.h"
#include "firefly.h"

int main() {
	int choice = 0;
	Firefly firefly("fireflyParams.txt");
	ParticleSwarm particleSwarm("swarmParams.txt");
	cout << "Which would you like to run?" << endl;
	cout << "1. Particle Swarm" << endl;
	cout << "2. Firefly" << endl;
	cin >> choice;
	switch (choice) {
	case 1:
		cout << "Running Particle Swarm" << endl;
		cout << "Particle Swarm returned: " << particleSwarm.runParticleSwarm();
		break;
	case 2:
		cout << "Running Firefly" << endl;
		cout << "Firefly Returned: " << firefly.runFirefly();
		break;
	default:
		cout << "Invalid entry" << endl;
	}
	return 0;
}