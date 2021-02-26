#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION__TRIDIAGONAL_MATRIX_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION__TRIDIAGONAL_MATRIX_
/*
The index of row and column begins with 0
The structure of the matrix is :
    b, c, 0, 0
    a, b, c, 0
    0, a, b, c
    0, 0, a, b
*/
#include <math.h>
#include <mkl.h>
#include <vector>

class Vector
{
private:
    const long length;
    float* data;

public:
    Vector(const std::vector<float>& data);
    Vector(const Vector& vec);
    ~Vector();

    float& operator [](const long& index);
};

class Tridiagonal
{
private:
    long length;
    float* B, * A, * C;

public:
    Tridiagonal(const long& length);
    ~Tridiagonal();

    float& operator ()(const long& r, const long& c);
    float* Diagonal(const int& offset = 0);
    long size();

    void multiple(float * x, float * res);
};



float Norm(float* A, float* B, const long& length);
void Chasing(Tridiagonal& A, float* f, const long& length);
#endif // !_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION__TRIDIAGONAL_MATRIX_