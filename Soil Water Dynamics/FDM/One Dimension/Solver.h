#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOLVER_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOLVER_

//Using the MKL to speed up the computation
#include <mkl.h>

//for solving the equations
#include "Tridiagonal.h"

//for the soil wate retention
#include "Soil.h"

//for the process of boudary condition
#include "Boundary.h"

//for the process of water sink 
//especially the plant uptake
#include "Sink.h"

#include "States.h"

#include <string>
#include <map>
#include <array>

class Water
{
public:
    virtual bool Solve(States& Old, States& New, const float& dt,
        Sink* sink, Boundary* upper, Boundary* lower) = 0;
    virtual void Flux(States& Old, States& New, const float& dt) = 0;
    virtual ~Water() {};
};

class Ion
{
private:
    float lambda, Do, a, b;
public:
    Ion(const std::array<float, 4>& parameters);
    void diffusion(float* velocity, float* theta, float* dsh, const long& length);
};

class Solute
{
protected:
    std::map<std::string, Ion*> ions;

public:

    virtual ~Solute() {};
    
    virtual bool Solve(States& Old, States& New, const float& dt, 
        Sink* sink, Boundary* upper, Boundary* lower) = 0;
};

class Heat
{

};
#endif // !_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOLVER_