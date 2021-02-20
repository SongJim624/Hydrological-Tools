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
};

void Chasing(Tridiagonal& A, float* f, const long& length);
#endif // !_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION__TRIDIAGONAL_MATRIX_