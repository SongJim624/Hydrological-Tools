#include "Tridiagonal.h"

Tridiagonal::Tridiagonal(const long& length)
    : length(length), B(new float[length]),
      A(new float[length - 1]), C(new float[length - 1])
{}

Tridiagonal::~Tridiagonal()
{
    delete[] A;
    delete[] B;
    delete[] C;
}

float& Tridiagonal:: operator ()(const long& r, const long& c)
{
    return r == c ? B[r] :
        r > c ? A[c] : C[r];
}

float* Tridiagonal::Diagonal(const int& offset)
{
    switch(offset)
    {
    case 0:
    {
        return B;
    }
    case 1:
    {
        return C;
    }
    case -1:
    {
        return A;
    }
    }
}

long Tridiagonal::size()
{
    return length;
}

void Tridiagonal::multiple (float* x, float * res)
{
/*
    res[0] = B[0] * x[0] + C[0] * x[1];

    for(long i = 1; i < length - 1; ++i)
    {
        res[i] = A[i - 1] * x[i - 1] + B[i] * x[i] + C[i] * x[i + 1];
    }

    res[length - 1] = A[length - 2] * x[length - 2] + B[length - 1] * x[length - 1];
*/
    float* temp = (float*)MKL_calloc(length, sizeof(float), 64);

    vsMul(length, B, x, res);

    vsMul(length - 1, C, x + 1, temp);
    vsAdd(length - 1, temp, res, res);

    vsMul(length - 1, A, x, temp);
    vsAdd(length - 1, res + 1, temp, res + 1);

    MKL_free(temp);
}

float Norm(float* A, float* B, const long& length)
{
    float norm = -INFINITY;

    for (long i = 0; i < length; ++i)
    {
        norm = fmaxf(norm, fabsf((A[i] - B[i]) / B[i]));
    }

    return norm;
}

void Chasing(Tridiagonal& A, float* b, const long& length)
{
    float* l = (float *)MKL_calloc(length - 1, sizeof(float), 64);
    cblas_scopy(length - 1, A.Diagonal(-1), 1, l, 1);

    float* m = (float *)MKL_calloc(length, sizeof(float), 64);
    cblas_scopy(length, A.Diagonal(0), 1, m, 1);
    
    float* u = nullptr;
    u = A.Diagonal(1);

    for (long i = 1; i < length; ++i)
    {
        l[i - 1] /= m[i - 1];
        m[i] -= u[i - 1] * l[i - 1];
    }

    for (long i = 1; i < length; ++i)
    {
        b[i] -= l[i - 1] * b[i - 1];
    }

    b[length - 1] /= m[length - 1];

    for (long i = length - 2; i != 0; --i)
    {
        b[i] = (b[i] - u[i]/*A(i, i + 1)*/ * b[i + 1]) / m[i];
    }

    b[0] = (b[0] - u[0] * b[1]) / m[0];

    MKL_free(m);
    MKL_free(l);
}