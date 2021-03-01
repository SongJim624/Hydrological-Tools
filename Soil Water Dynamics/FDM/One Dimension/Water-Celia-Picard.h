#include "Solver.h"

class Water_Celia_Picard : public Water
{
private:
    Soil* soil;

    long length;

    const float maxIter;
    const float tol;

    Tridiagonal A;
    float* f, * theta, * capacity, * conductivity, * s;
    float* temp;
//        *
//        positive_ones, * negative_ones, * half, * dt;

public:
    Water_Celia_Picard(Soil* soil, const long& maxIter = 20, const float& tol = 0.01);
    ~Water_Celia_Picard();

    virtual bool Solve(States& Old, States& New, const float& dt,
        Sink* sink, Boundary* upper, Boundary* lower);
    virtual void Flux(States& Old, States& New, const float& dt);
};