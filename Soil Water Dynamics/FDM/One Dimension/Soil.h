#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOIL_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOIL_
class Soil
{
public:
    virtual ~Soil() {};

    virtual void theta(float* h, float* res) = 0;
    virtual void capacity(float* h, float* theta, float* res) = 0;
    virtual void conductivity(float* theta, float* res) = 0;
    virtual void potential(float* h, float* res) = 0;

    virtual long size() = 0;
};
#endif //!_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOIL_