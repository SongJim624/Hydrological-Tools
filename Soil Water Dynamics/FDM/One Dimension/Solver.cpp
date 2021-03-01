#include "Solver.h"

Ion::Ion(const std::array<float, 4>& parameters)
    : lambda(parameters[0]), Do(parameters[1]), a(parameters[2]), b(parameters[3])
{
}

void Ion::diffusion(float* velocity, float* theta, float* dsh, const long& length)
{
    float* temp = (float*)MKL_calloc(length, sizeof(float), 64);
    for (long i = 0; i < length; ++i)
    {
        temp[i] = b;
    }

    vsAbs(length, velocity, dsh);
    
    vsMul(length, theta, temp, temp);
    vsExp(length, temp, temp);

    cblas_saxpby(length, a * Do, temp, 1, lambda, dsh, 1);
    
    MKL_free(temp);
}