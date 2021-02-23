#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOLVER_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOLVER_

//Using the MKL to speed up the computation
//#include <mkl.h>

//for solving the equations
#include "Tridiagonal.h"

//for the soil wate retention
#include "Soil.h"

//for the process of boudary condition
#include "Boundary.h"

//for the process of water sink 
//especially the plant uptake
#include "Sink.h"

class Solver
{
private:
    Soil* soil;
    
    long length;

    const float maxIter;
    const float tol;

    Tridiagonal A;
    float* f, * theta, * capacity, * conductivity, * s;

public:
    Solver(Soil* soil, const long& maxIter = 20, const float& tol = 0.01);
    ~Solver();

    bool Solve(float* h, float* theta0, const float& dt,
        Sink* sink, Boundary* upper, Boundary* lower);
    void Flux(const float& theta0, const float& dt, float* flux);
};

#endif // !_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SOLVER_