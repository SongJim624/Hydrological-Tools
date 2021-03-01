#include "Boundary.h"

Boundary::Boundary(const long& length)
    : h(nullptr), capacity(nullptr), conductivity(nullptr), 
    theta(nullptr), theta0(nullptr), sink(nullptr), length(length){}

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

First::First(const long& length, const float& head)
    : Boundary(length), head(head){}

void First::Upper(Tridiagonal& A, float* f, const float& dt)
{
    A(0, 0) = 1;
    f[0] = head;

//    f[1] -= f[0] * A(1, 0);
//    A(1, 0) = 0;
}

void First::Lower(Tridiagonal& A, float* f, const float& dt)
{
    A(length - 1, length - 1) = 1;
    f[length - 1] = head;

//    f[length - 2] -= f[length - 1] * A(length - 2, length - 1);
//    A(length - 2, length - 1) = 0;
}

void First::Modify(Tridiagonal& A, float* f, const float& dt)
{
    length == 0 ? Upper(A, f, dt) : Lower(A, f, dt);
}

Second::Second(const long& length, const float& flux)
    : Boundary(length), flux(flux) {}

void Second::Upper(Tridiagonal& A, float* f, const float& dt)
{
    A(0, 0) = 0.5 * (capacity[0] + dt * (conductivity[0] + conductivity[1]));
    A(0, 1) = -0.5 * (conductivity[0] + conductivity[1]) * dt;

    f[0] = 0.5 * (h[0] * capacity[0] - (theta[0] - theta0[0]) - dt * (conductivity[1] + conductivity[0]))
        + flux * dt - 0.5 * sink[0] * dt;
}
//Here need the further modify
void Second::Lower(Tridiagonal& A, float* f, const float& dt)
{
    A(length - 1, length - 1) = 0.5 * dt * (conductivity[length - 1] + conductivity[length - 2]);
    A(length - 1, length - 2) = -0.5 * (conductivity[length - 1] + conductivity[length - 2]) * dt;

    f[length - 1] = 0.5 * dt * (conductivity[length - 1] + conductivity[length - 2]) + flux * dt;
}

void Second::Modify(Tridiagonal& A, float* f, const float& dt)
{
    length == 0 ? Upper(A, f, dt) : Lower(A, f, dt);
}