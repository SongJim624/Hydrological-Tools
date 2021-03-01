#include <mkl.h>
#include <string>
#include <map>

class States
{
public:
	static long length;

	float* potential;
	float* theta;
	float* temperature;
	float* flux, * flux_middle;
	std::map<std::string, float*> concentrations;

public:
	States();
	States(const States& states);
	~States();
};