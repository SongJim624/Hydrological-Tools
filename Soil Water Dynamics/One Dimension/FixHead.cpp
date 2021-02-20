#include "Boundary.h"

FixUpper::FixUpper(const float& h) : head(h)
{}

void FixUpper::Modify(Tridiagonal& A, float* f)
{
    A(0, 0) = 1;
    f[0] = head;

    f[1] -= f[0] * A(1, 0);
    A(1, 0) = 0;
}