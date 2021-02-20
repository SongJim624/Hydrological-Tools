#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_BOUNDARY_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_BOUNDARY_

#include "Tridiagonal.h"
#include <map>


enum class Type
{
    first = 1,
    second = 2,
    third = 3
};

class Boundary
{
protected:
    float* h, * theta, * theta0, * capacity, * conductivity, * sink;

public:
    Boundary();
    ~Boundary();

    virtual void Modify(Tridiagonal& A, float* f) = 0;

    virtual void Attach(float* h,  float* capacity, float* conductivity, 
        float* theta, float* theta0, float* sink);
};

/*
* 
* 
class Rain : public Boundary
{

};

class Irrigation : public Boundary
{

};


class Evaporation : public Boundary
{


};
*/


class FixUpper : public Boundary
{
private:
    const float head;
public:
    FixUpper(const float& h);
    virtual void Modify(Tridiagonal& A, float *f);
};

class FixLower : public Boundary
{
private:
    const float head;
    const long length;

public:
    FixLower(const float& h, const long& length);
    virtual void Modify(Tridiagonal& A, float* f);
};
#endif // !_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_BOUNDARY_
