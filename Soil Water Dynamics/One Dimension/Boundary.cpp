#include "Boundary.h"

Boundary::Boundary()
    : h(nullptr), capacity(nullptr), conductivity(nullptr), 
    theta(nullptr), theta0(nullptr), sink(nullptr){}

Boundary::~Boundary()
{
    h = capacity = conductivity = theta = theta0 = sink = nullptr;
}

void Boundary::Attach(float* h, float* capacity, float* conductivity,
    float* theta, float* theta0, float* sink)
{
    this->h = h == this->h ? this->h : h;
    this->capacity = capacity == this->capacity ? this->capacity : capacity;
    this->conductivity = conductivity == this->conductivity ? this->conductivity : conductivity;

    this->theta = theta == this->theta ? this->theta : theta;
    this->theta0 = theta0 == this->theta0 ? this->theta0 : theta0;
    this->sink = sink == this->sink ? this->sink : sink;
}
/*
float Evporation(const float& Kex, const float& ETo, const float& CC,
    const float &h, const float& hc, const float& hcc)
{
//The unit of the ETo is mm, and theresult should be changed to cm in this process
    if (h > hc)
    {
        return  ETo * Kex * (1 - CC) / 240;
    }
    else
    {
        return ETo * max(0, log(hcc / h) / log(hcc / hc)) / 240;
    }
}

//first
inline void Upper(const float& h, Tridiagonal& A, float* f)
{
    A(0, 0) = 1;
    f[0] = h;

    f[1] -= f[0] * matrix(1, 0);
    A(1, 0) = 0;
}

//second
inline void Upper(const float& flux, float* ks, const float& cap, const float& theta,
    const float& theta0, const float& sink, const float &dt, Triangle& A, float* f)
{
    A(0, 0) = 0.5 * (cap + dt * (ks[0] + ks[1]));
    A(0, 1) = -0.5 * (ks[0] + ks[1]) * dt;
    f[0] = 0.5 * (h[0] * cap - (theta - theta0) - dt * (ks[1] + ks[0]))
        + flux * dt - 0.5 * sink[0] * dt;
}

void Upper_Third()
{

}

//first
inline void Lower(const float& h, const long& length, Triangle& A, float* f)
{
    A(length - 1, length - 1) = 1;
    f[length - 1] = h;

    f[length - 2] -= f[length - 1] * A(length - 2, length - 1);
    A(length - 2, length - 1) = 0;
}

inline void Lower_Second()
{

}

inline void Lower_Third()
{

}

inline void Free_Drainage(const long& length, Triangle& A, float* f)
{
    A(length - 1, length - 1) = -A(length - 2, length - 1);
    A(length - 1, length - 2) = A(length - 2, length - 1);
    f[length - 1] = 0;
}
*/