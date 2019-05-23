#include "Density.h"

int Density::readXYZR(string fn, double gs, double ex, double iso)
{
	ifstream in(fn.c_str());

	string str;
	while (getline(in, str)) {
		stringstream ss(str);

		double x, y, z, r;
		ss >> x >> y >> z >> r;
		//cout << "Atoms: " << x << " " << y << " " << z << " " << r << endl;
		vec3 v(x, y, z);
		
		atoms.push_back(v);
		radius.push_back(r);
	}

	gridsize = gs;
	extension = ex;
	isovalue = iso;

	return 1;
}

int Density::readDataset(string fn, double gs, double ex, double iso)
{
	ifstream in(fn.c_str());

	string str;
	while (getline(in, str)) {
		stringstream ss(str);

		double x, y, z;
		x = stod(str.substr(31, 8));
		y = stod(str.substr(39, 8));
		z = stod(str.substr(47, 8));

		vec3 v(x, y, z);

		atoms.push_back(v);
		radius.push_back(1);
	}

	gridsize = gs;
	extension = ex;
	isovalue = iso;

	return 1;
}

int Density::readVox(string fn, double iso)
{
	// read vox data
	ifstream in(fn.c_str());
	
	string str;

	getline(in, str);
	stringstream ss0(str);
	ss0 >> xmin >> ymin >> zmin;

	getline(in, str);
	stringstream ss1(str);
	ss1 >> xmax >> ymax >> zmax;

	getline(in, str);
	stringstream ss2(str);
	ss2 >> xdim >> ydim >> zdim;

	gridsize = (xmax - xmin) / xdim;

	while (getline(in, str)) {
		stringstream ss(str);
		double val;
		ss >> val;

		voxs.push_back(val);
	}

	isovalue = iso;

	return 1;
}

int Density::createGrid()
{
	// find bounding box;
	xmin = atoms[0].x();
	ymin = atoms[0].y();
	zmin = atoms[0].z();
	xmax = atoms[0].x();
	ymax = atoms[0].y();
	zmax = atoms[0].z();

	for (int i = 1; i < atoms.size(); ++i) {
		if (xmin > atoms[i].x())
			xmin = atoms[i].x();
		if (ymin > atoms[i].y())
			ymin = atoms[i].y();
		if (zmin > atoms[i].z())
			zmin = atoms[i].z();

		if (xmax < atoms[i].x())
			xmax = atoms[i].x();
		if (ymax < atoms[i].y())
			ymax = atoms[i].y();
		if (zmax < atoms[i].z())
			zmax = atoms[i].z();
	}

	xmin -= extension;
	ymin -= extension;
	zmin -= extension;
	xmax += extension;
	ymax += extension;
	zmax += extension;
	cout << xmin << " " << ymin << " " << zmin << endl;
	cout << xmax << " " << ymax << " " << zmax << endl;
	cout << gridsize << endl;

	cout << xmax - xmin << endl;
	xdim = ceil((xmax - xmin) / gridsize);
	ydim = ceil((ymax - ymin) / gridsize);
	zdim = ceil((zmax - zmin) / gridsize);
	cout << xdim << " " << ydim << " " << zdim << endl;

	xmax = xmin + xdim * gridsize;
	ymax = ymin + ydim * gridsize;
	zmax = zmin + zdim * gridsize;

	voxs.resize(xdim*ydim*zdim, 0);

	return 1;
}

int Density::fillDensity()
{

	for (int l = 0; l < atoms.size(); ++l) {
		for (int k = 0; k < zdim; ++k) {
			for (int j = 0; j < ydim; ++j) {
				for (int i = 0; i < xdim; ++i) {
					vec3 pos(xmin + (i)*gridsize, ymin + (j)*gridsize, zmin + (k)*gridsize);
					voxs[i + j * xdim + k * xdim*ydim] += getDensity(pos, atoms[l], 2*radius[l]);
					//cout << voxs[i + j * xdim + k * xdim*ydim] << endl;
				}
			}
		}
	}

	return 1;
}

int Density::writeDensity()
{
	ofstream out("density.vox");

	out << xmin << " " << ymin << " " << zmin << endl;
	out << xmax << " " << ymax << " " << zmax << endl;
	out << xdim << " " << ydim << " " << zdim << endl;

	for (int i = 0; i < voxs.size(); ++i) {
		if (voxs[i] > 1)
			voxs[i] = 1;
		if (voxs[i] < 0)
			voxs[i] = 0;

		out << voxs[i] << endl;
	}

	return 1;
}

int Density::createLevelSetSurface()
{
	openvdb::initialize();
	
	grid = openvdb::DoubleGrid::create();

	//openvdb::Mat4d transformation;
	//transformation.setZero();
	//transformation(0, 0) = gridsize;
	//transformation(1, 1) = gridsize;
	//transformation(2, 2) = gridsize;
	//transformation(3, 3) = 1;
	//transformation(0, 3) = xmin;
	//transformation(1, 3) = ymin;
	//transformation(2, 3) = zmin;

	//grid->setTransform(openvdb::math::Transform::createLinearTransform(transformation));

	openvdb::DoubleGrid::Accessor accessor = grid->getAccessor();
	
	for (int k = 0; k < zdim; ++k) {
		for (int j = 0; j < ydim; ++j) {
			for (int i = 0; i < xdim; ++i) {
				openvdb::Coord xyz(i, j, k);
				//cout << i<<" " << j << " " << k<< endl;
				//cout << xdim<<" " << ydim << " " << zdim << endl;
				//cout << voxs[i + j*xdim + k*xdim*ydim] << endl;
				accessor.setValue(xyz, voxs[i + j*xdim + k*xdim*ydim]);
			}
		}
	}

	openvdb::tools::volumeToMesh<openvdb::DoubleGrid>(*grid, vertices, quads, isovalue);
	//cout << "123123"<< endl;
	for (int i = 0; i < vertices.size(); ++i) {
		vertices[i] *= gridsize;
		vertices[i] += openvdb::Vec3s(xmin, ymin, zmin);
	}
}

void Density::createSamples()
{
	ofstream out("sample.xyzr");
	ofstream obj("sample.obj");

	int spacing = 4;
	for (int k = 0; k < zdim; ++k) {
		for (int j = 0; j < ydim; ++j) {
			for (int i = 0; i < xdim; ++i) {
				//openvdb::Coord xyz(i, j, k);
				//cout << i<<" " << j << " " << k<< endl;
				//cout << xdim<<" " << ydim << " " << zdim << endl;
				//cout << voxs[i + j*xdim + k*xdim*ydim] << endl;
				if (i%spacing != 0)
					continue;
				if (j%spacing != 0)
					continue;
				if (k%spacing != 0)
					continue;
				
				if (voxs[i + j*xdim + k*xdim*ydim] > isovalue) {
					double x = xmin + gridsize*i;
					double y = ymin + gridsize*j;
					double z = zmin + gridsize*k;
					out << x << " " << y << " " << z << " 1" << endl;
					obj << "v " << x << " " << y << " "<<z << endl;
				}
			}
		}
	}

	out.close();
}

int Density::writeMesh()
{
	//ofstream out("mesh.obj");

	//for (int i = 0; i < vertices.size(); ++i) {
	//	out << "v " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << endl;
	//}

	//for (int i = 0; i < quads.size(); ++i) {
	//	out << "f " << quads[i].x() + 1 << " " << quads[i].y() + 1 << " " << quads[i].z() + 1 << " " << quads[i].w() + 1 << endl;
	//}

	ofstream cloudout("cloud.txt");
	for (int i = 0; i < vertices.size(); ++i) {
		cloudout << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << endl;
	}

	ofstream out("mesh.off");
	out << "COFF" << endl;
	out << vertices.size() << " " << quads.size() * 2 << " 0" << endl;
	
	for (int i = 0; i < vertices.size(); ++i) {
		out << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << " "
			<< "255 255 255 255" << endl;
	}

	for (int i = 0; i < quads.size(); ++i) {
		out << "3 " << quads[i].x() << " " << quads[i].y() << " " << quads[i].z() << endl;
		out << "3 " << quads[i].x() << " " << quads[i].z() << " " << quads[i].w() << endl;
	}

	return 1;
}

double Density::getDensity(vec3 pos, vec3 anchor, double radius)
{
	vec3 d = pos - anchor;

	double val = exp(-pow((sqrt(d.dot(d))) / radius, 2));

	return val;
}