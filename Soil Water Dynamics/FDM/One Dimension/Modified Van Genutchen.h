#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_VGM_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_VGM_

#include "Soil.h"
#include <mkl.h>
#include <list>
#include <array>
#include <vector>

class VGM : public Soil
{
private:
    long length;
    float* thetas, * thetar, * alpha, * n, * ksat, * hs;
    float* m, * thetam, * inter, * F1, * temp, * positive_ones, * negative_ones;
    bool* label;

private:
    void Convert(const std::list<float>& compart, const std::list<float>& value, float*& res);

public:

    VGM(const std::array<std::list<float>, 7>& soil);
    virtual ~VGM();

    virtual void theta(float* h, float* res);
    virtual void capacity(float* h, float* theta, float* res);
    virtual void conductivity(float* theta, float* res);
    virtual void potential(float* theta, float* res);

    virtual long size();
};
#endif//! _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_VGM_

/*
class VGM_Interpolation : public Soil
{
private:
    const long length, N;

    long* indices;

    std::map<long, std::array<float*, 4>> dictionaries;
    std::map<long, std::array<float, 6>> soils;

private:
    //inter varaible when the potential exceeds the table;

private:
    void Convert(const std::list<float>& compart, const std::list<float>& value, float*& res);

public:
    VGM_Interpolation();
    ~VGM_Interpolation();

    virtual void theta(float* h, float* res);
    virtual void capacity(float* h, float* theta, float* res);
    virtual void conductivity(float* theta, float* res);
    virtual void potential(float* theta, float* res);

    virtual long size();
};
*/