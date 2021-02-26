#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOLUTE_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOLUTE_
class Solute
{
private:
	const long length;

public:

	float* concentration;

	Solute(const long& length);
	virtual ~Solute() {}

	virtual void difussion(float* vocelity, float * theta, float* disperation) = 0;
};
#endif // !_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOLUTE_