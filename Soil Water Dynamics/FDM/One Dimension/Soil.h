#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOIL_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOIL_

#include <list>
#include <vector>
#include <array>
#include <mkl.h>

class Soil
{
public:
    virtual void theta(float* h, float* res) = 0;
    virtual void capacity(float* h, float* theta, float* res) = 0;
    virtual void conductivity(float* theta, float* res) = 0;
    virtual void potential(float* h, float* res) = 0;

    virtual long size() = 0;
};

class VGM : public Soil
{
private:
    long length;
    float* thetas, * thetar, * alpha, * n, * ksat, * hs;
    float* m, * thetam, * inter, * F1;
    bool* label;

private:
    void Convert(const std::list<float>& compart, const std::list<float>& value, float*& res);

public:

    VGM(const std::array<std::list<float>, 7>& soil);
    ~VGM();

    virtual void theta(float* h, float* res);
    virtual void capacity(float* h, float* theta, float* res);
    virtual void conductivity(float* theta, float* res);
    virtual void potential(float* theta, float* res);

    virtual long size();
};
#endif //!_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOIL_