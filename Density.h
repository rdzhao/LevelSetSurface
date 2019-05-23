#pragma once

#ifndef _CAVITY_DETECTION_H_ 
#define _CAVITY_DETECTION_H_
 
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <ctime>

#include <openvdb/openvdb.h>
#include <openvdb/Grid.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/VolumeToMesh.h>

typedef openvdb::Vec3d vec3;

using namespace std;

class Density
{
public:
	int readXYZR(string fn, double gs, double ex, double iso);
	int readDataset(string fn, double gs, double ex, double iso);
	int readVox(string fn, double iso);

	int createGrid();
	int fillDensity();

	int writeDensity();
	int writeMesh();

	int createLevelSetSurface();
	void createSamples();
	double getDensity(vec3 pos, vec3 anchor, double miu);

	vector<vec3> atoms;
	vector<double> radius;
	double gridsize;
	double extension;
	double isovalue;


	double xmin, ymin, zmin;
	double xmax, ymax, zmax;
	int xdim, ydim, zdim;

	vector<double> voxs;

	openvdb::DoubleGrid::Ptr grid;

	vector<openvdb::Vec3s> vertices;
	vector<openvdb::Vec4I> quads;
};

#endif // !_CAVITY_DETECTION_H_ 
