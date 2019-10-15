#ifndef _VECTOR1_H_
#define _VECTOR1_H_

class vector{
public:
	// Constructors
	vector(): x(0.0), y(0.0), z(0.0){}
	vector(double x1, double y1, double z1): x(x1), y(y1), z(z1){}

	// Fields
	double x;
	double y;
	double z;
};

#endif