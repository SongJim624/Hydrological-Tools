#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_BOUNDARY_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_BOUNDARY_
/*
* The basic type to prcess the boundary condition is the class Boundary
* In this file, two basic types was provided
* Moreover, the length in the constructor of Boudnary is the indicator of the boudanry location
* length = 0 means the upper boundary, length != 0 means the lower boudnary
* The third type boudnary is not provided because the realization of it involves many custom functions
* The more complex boundary should be realized by the composition of the three types
*/

#include "Tridiagonal.h"
#include <map>
//
class Boundary
{
protected:
    float* h, * theta, * theta0, * capacity, * conductivity, * sink;
    const long length;
public:
    Boundary(const long& length = 0);
   virtual  ~Boundary();

    virtual void Modify(Tridiagonal& A, float* f, const float& dt) = 0;

    virtual void Attach(float* h,  float* capacity, float* conductivity, 
        float* theta, float* theta0, float* sink);
};

class First : public Boundary
{
protected:
    const float head;

protected:
    void Upper(Tridiagonal& A, float* f, const float& dt);
    void Lower(Tridiagonal& A, float* f, const float& dt);

public:
    First(const long& length, const float& head);
    virtual void Modify(Tridiagonal& A, float* f, const float& dt);
};

class Second : public Boundary
{
private:
    const float flux;

private:
    void Upper(Tridiagonal& A, float* f, const float& dt);
    void Lower(Tridiagonal& A, float* f, const float& dt);

public:
    Second(const long& length, const float& flux);
    virtual void Modify(Tridiagonal& A, float* f, const float& dt);
};
#endif // !_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_BOUNDARY_
