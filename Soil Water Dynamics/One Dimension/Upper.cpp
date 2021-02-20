#include "Boundary.h"
/*
void Upper::Modify(Tridiagonal& A, float* f)
{
    switch (type)
    {
    case Type::first:
        First(h, A, f);
        break;
    case Type::second:
        Second(flux, Tridiagonal & A, float* f)
        break;
    case Type::third:
        break;
    }
}

void Upper::First(const float& h, Tridiagonal&A, float* f)
{
    A(0, 0) = 1;
    f[0] = h;

    f[1] -= f[0] * A(1, 0);
    A(1, 0) = 0;
}

inline void Upper::Second(const float& flux, const float& dt, Tridiagonal& A, float* f)
{
    A(0, 0) = 0.5 * (capacity[0] + dt * (conductivity[0] + conductivity[1]));
    A(0, 1) = -0.5 * (conductivity[0] + conductivity[1]) * dt;
    f[0] = 0.5 * (h[0] * capacity[0] - (theta[0] - theta0[0]) - dt * (conductivity[1] + conductivity[0]))
        + flux * dt - 0.5 * sink[0] * dt;
}
*/