#include "Solver.h" 

inline float Norm(float* A, float* B, const long& length)
{
    float norm = -INFINITY;

    for (long i = 0; i < length; ++i)
    {
        norm = fmaxf(norm, abs((A[i] - B[i]) / B[i]));
    }

    return norm;
}

Solver::Solver(Soil* soil, const long& maxIter, const float& tol)
    : soil(soil), length(soil->size()), A(Tridiagonal(length)),
    f(new float[length]), theta(new float[length]), capacity(new float[length]), 
    conductivity(new float[length]), s(new float[length]), maxIter(maxIter), tol(tol){}

Solver::~Solver()
{
    delete[] f;
    delete[] theta;
    delete[] capacity;
    delete[] conductivity;
    delete[] s;
}


bool Solver::Solve(float* h, float* theta0, const float &dt, 
    Sink* sink, Boundary* upper, Boundary* lower)
{
    size_t count = 0;

    upper->Attach(h, capacity, conductivity, theta, theta0, s);
    lower->Attach(h, capacity, conductivity, theta, theta0, s);


    while (++count < maxIter)
    {
        soil->theta(h, theta);
        soil->capacity(h, theta, capacity);
        soil->conductivity(theta, conductivity);

        sink->Modify(h, s, length);

        A(0, 1) = 0.0;
        A(length - 1, length - 2) = 0.0;

        for (size_t i = 1; i < length - 1; ++i)
        {
            A(i, i - 1) = -0.5 * dt * (conductivity[i - 1] + conductivity[i]);
            A(i, i + 1) = -0.5 * dt * (conductivity[i] + conductivity[i + 1]);
            A(i, i) = capacity[i] - A(i, i - 1) - A(i, i + 1);
        }

        for (size_t i = 1; i < length - 1; ++i)
        {
            f[i] = h[i] * capacity[i] - (theta[i] - theta0[i]) - 
                0.5 * dt * (conductivity[i + 1] - conductivity[i - 1]) - s[i] * dt;
        }

        upper->Modify(A, f);
        lower->Modify(A, f);

        Chasing(A, f, length);

        if (Norm(h, f, length) < tol)
        {
            return true;
        }

        for (long i = 0; i < length; ++i)
        {
            h[i] = f[i];
        }
    }

    return false;
}