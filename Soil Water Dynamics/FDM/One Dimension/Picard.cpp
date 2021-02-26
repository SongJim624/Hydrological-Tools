#include "Solver.h" 

Picard::Picard(Soil* soil, const long& maxIter, const float& tol)
    : soil(soil), length(soil->size()), A(Tridiagonal(length)),
    f((float*)MKL_calloc(length, sizeof(float), 64)), 
    theta((float*)MKL_calloc(length, sizeof(float), 64)),
    capacity((float*)MKL_calloc(length, sizeof(float), 64)),
    conductivity((float*)MKL_calloc(length, sizeof(float), 64)),
    s((float*)MKL_calloc(length, sizeof(float), 64)), 
    temp((float*)MKL_calloc(length, sizeof(float), 64)), 
    positive_ones((float*)MKL_calloc(length, sizeof(float), 64)), 
    negative_ones((float*)MKL_calloc(length, sizeof(float), 64)), 
    half((float*)MKL_calloc(length, sizeof(float), 64)),
    dt((float*)MKL_calloc(length, sizeof(float), 64)), maxIter(maxIter), tol(tol)
{
    for (long i = 0; i < length; ++i)
    {
        positive_ones[i] = 1;
        negative_ones[i] = -1;
        half[i] = 0.5;
    }
}

Picard::~Picard()
{
    MKL_free(f);
    MKL_free(theta);
    MKL_free(capacity);
    MKL_free(conductivity);
    MKL_free(s);

    MKL_free(temp);
    MKL_free(positive_ones);
    MKL_free(negative_ones);
    MKL_free(half);
    MKL_free(dt);
}

bool Picard::Solve(float* h, float* theta0, const float &dt, 
    Sink* sink, Boundary* upper, Boundary* lower)
{
    size_t count = 0;

    for (long i = 0; i < length; ++i)
    {
        this->dt[i] = dt;
    }

    upper->Attach(h, capacity, conductivity, theta, theta0, s);
    lower->Attach(h, capacity, conductivity, theta, theta0, s);

    while (++count < maxIter)
    {
        soil->theta(h, theta);
        soil->capacity(h, theta, capacity);
        soil->conductivity(theta, conductivity);

        sink->Modify(h, s, length);

// A
        vsAdd(length - 1, conductivity, conductivity + 1, A.Diagonal(-1));
        vsMul(length - 1, half, A.Diagonal(-1), A.Diagonal(-1));
        vsMul(length - 1, this->dt, A.Diagonal(-1), A.Diagonal(-1));
        vsMul(length - 1, negative_ones, A.Diagonal(-1), A.Diagonal(-1));

// C
        cblas_scopy(length - 1, A.Diagonal(-1), 1, A.Diagonal(1), 1);

// B
        vsSub(length - 2, capacity + 1, A.Diagonal(-1), A.Diagonal(0) + 1);
        vsSub(length - 2, A.Diagonal(0) + 1, A.Diagonal(1) + 1, A.Diagonal(0) + 1);

// F
        vsMul(length - 2, h + 1, capacity + 1, f + 1);

        vsSub(length - 2, theta + 1, theta0 + 1, temp + 1);
        vsSub(length - 2, f + 1, temp + 1, f + 1);

        vsSub(length - 2, conductivity + 2, conductivity, temp + 1);
        vsMul(length - 2, temp + 1, this->dt, temp + 1);
        vsMul(length - 2, temp + 1, half + 1, temp + 1);
        vsSub(length - 2, f + 1, temp + 1, f + 1);

        vsMul(length - 2, s + 1, this->dt + 1, temp + 1);
        vsSub(length - 2, f + 1, temp + 1, f + 1);
        
        A(0, 1) = 0.0;
        A(length - 1, length - 2) = 0.0;

        upper->Modify(A, f, dt);
        lower->Modify(A, f, dt);

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

void Picard::Flux(const float& theta0, const float& dt, float* flux)
{
    //After solving one time successfully, the f euqals to the water potential h
    //to reduce the input variable, the f is used instead of h.

	flux[0] = -0.5 * (conductivity[0] + conductivity[1]) * (f[1] - f[0] - 1) - 0.5 * (theta[0] - theta0)/ dt - 0.5 * s[0];

	for (long i = 1; i < length - 1; ++i)
	{
		flux[i] = 0.5 * (
			-0.5 * (conductivity[i] + conductivity[i + 1]) * (f[i + 1] - f[i] - 1)
			- 0.5 * (conductivity[i] + conductivity[i - 1]) * (f[i] - f[i - 1] - 1));
	}

	flux[length - 1] = -0.5 * (conductivity[length - 1] + conductivity[length - 2]) * (f[length - 1] - f[length - 2] - 1);
}