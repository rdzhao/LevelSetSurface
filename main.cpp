#include "Density.h"

 
int main(int argc, char* argv[])
{
	Density den;
	int mode;

	mode = stoi(argv[1]);

	if (mode == 0) { // for cryoem
		//den.readXYZR(argv[1], stod(argv[2]), stod(argv[3]), stod(argv[4]));

		//den.readDataset(argv[1], stod(argv[2]), stod(argv[3]), stod(argv[4]));
		//den.createGrid();
		//den.fillDensity();
		//den.writeDensity();
		//den.createLevelSetSurface();
		//den.writeMesh();

		den.readVox(argv[2], stod(argv[3]));
		den.createLevelSetSurface();
		den.createSamples();
		den.writeMesh();
	}
	
	if (mode == 1) { // for xyzr
		den.readXYZR(argv[2], stod(argv[3]), stod(argv[4]), stod(argv[5]));
		den.createGrid();
		den.fillDensity();
		//den.writeDensity();
		den.createLevelSetSurface();
		den.writeMesh();
	}

	return 1;
}