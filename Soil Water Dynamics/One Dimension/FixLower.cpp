#include "Boundary.h"

FixLower::FixLower(const float& h, const long& length)
	: head(h), length(length){}

void FixLower::Modify(Tridiagonal& A, float* f)
{
    A(length - 1, length - 1) = 1;
    f[length - 1] = head;

    f[length - 2] -= f[length - 1] * A(length - 2, length - 1);
    A(length - 2, length - 1) = 0;
}