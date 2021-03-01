#include "Water-Celia-Picard.h" 

Water_Celia_Picard::Water_Celia_Picard(Soil* soil, const long& maxIter, const float& tol)
    : soil(soil), length(soil->size()), A(Tridiagonal(length)), maxIter(maxIter), tol(tol),
    f((float*)MKL_calloc(length, sizeof(float), 64)), 
    theta((float*)MKL_calloc(length, sizeof(float), 64)),
    capacity((float*)MKL_calloc(length, sizeof(float), 64)),
    conductivity((float*)MKL_calloc(length, sizeof(float), 64)),
    s((float*)MKL_calloc(length, sizeof(float), 64)), 
    temp((float*)MKL_calloc(length, sizeof(float), 64)){}

Water_Celia_Picard::~Water_Celia_Picard()
{
    MKL_free(f);
    MKL_free(theta);
    MKL_free(capacity);
    MKL_free(conductivity);
    MKL_free(s);

    MKL_free(temp);
}

bool Water_Celia_Picard::Solve(States& Old, States& New, const float &dt,
    Sink* sink, Boundary* upper, Boundary* lower)
{
    size_t count = 0;

    upper->Attach(New.potential, capacity, conductivity, New.theta, Old.theta, s);
    lower->Attach(New.potential, capacity, conductivity, New.theta, Old.theta, s);

    while (++count < maxIter)
    {
        soil->theta(New.potential, New.theta);
        soil->capacity(New.potential, New.theta, capacity);
        soil->conductivity(New.theta, conductivity);

        sink->Modify(New.potential, s, length);

// A
        cblas_scopy(length - 1, conductivity, 1, A.Diagonal(-1), 1);
        cblas_saxpby(length - 1, -0.5f * dt, conductivity + 1, 1, -0.5f * dt, A.Diagonal(-1), 1);

// C
        cblas_scopy(length - 1, A.Diagonal(-1), 1, A.Diagonal(1), 1);

// B
        vsSub(length - 2, capacity + 1, A.Diagonal(-1), A.Diagonal(0) + 1);
        vsSub(length - 2, A.Diagonal(0) + 1, A.Diagonal(1) + 1, A.Diagonal(0) + 1);

// F
        vsMul(length - 2, New.potential + 1, capacity + 1, f + 1);

        vsSub(length - 2, New.theta + 1, Old.theta + 1, temp + 1);
        vsSub(length - 2, f + 1, temp + 1, f + 1);

        cblas_scopy(length - 2, conductivity, 1, temp + 1, 1);
        cblas_saxpby(length - 2, 0.5f * dt, conductivity + 2, 1, -0.5f * dt, temp + 1, 1);     
        vsSub(length - 2, f + 1, temp + 1, f + 1);

        cblas_saxpy(length - 2, -dt, s + 1, 1, f + 1, 1);

        A(0, 1) = 0.0;
        A(length - 1, length - 2) = 0.0;

        upper->Modify(A, f, dt);
        lower->Modify(A, f, dt);

        Chasing(A, f, length);

        if (Norm(New.potential, f, length) < tol)
        {
            return true;
        }

        cblas_scopy(length, f, 1, New.potential, 1);
    }

    return false;
}

void Water_Celia_Picard::Flux(States& Old, States& New, const float& dt)
{
//After solving one time successfully, the f euqals to the water potential h
//to reduce the input variable, the f is used instead of h.
    if (&Old == &New)
    {
        soil->theta(New.potential, theta);
        soil->conductivity(theta, conductivity);
    }

//flux of the middle node
    cblas_scopy(length, conductivity, 1, temp, 1);
    vsSub(length - 1, New.potential, New.potential + 1, New.flux_middle);    
    cblas_saxpby(length - 1, 0.5f, conductivity, 1, 0.5f, temp + 1, 1);
    vsMul(length - 1, New.flux_middle, temp + 1, New.flux_middle);
    cblas_saxpy(length - 1, 1.0f, temp + 1, 1, New.flux_middle, 1);

//flux of the nodes
    if (dt == 0)
    {
        New.flux[0] = -0.5f * (conductivity[0] + conductivity[1]) * (New.potential[1] - New.potential[0] - 1.0f) + s[0];
    }
    else
    {
        New.flux[0] = -0.5f * (conductivity[0] + conductivity[1]) * (New.potential[1] - New.potential[0] - 1.0f) + (New.theta[0] - Old.theta[0]) / dt + s[0];
    }
    cblas_scopy(length - 2, New.flux_middle + 1, 1, New.flux + 1, 1);
    cblas_saxpby(length - 2, 0.5f, New.flux_middle, 1, 0.5f, New.flux + 1, 1);

	New.flux[length - 1] = -0.5f * (conductivity[length - 1] + conductivity[length - 2]) * (New.potential[length - 1] - New.potential[length - 2] - 1.0f);
}