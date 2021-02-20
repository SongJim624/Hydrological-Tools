#define _Solver_Test_

#ifdef _Tridiagonal_Test_

#include "Tridiagonal.h"

void main()
{
    Tridiagonal A(4);

    A(0, 0) =  2; A(0, 1) = -1;
    A(1, 0) = -1; A(1, 1) =  3; A(1, 2) = -2;
                  A(2, 1) = -1; A(2, 2) =  2; A(2, 3) = -1;
                                A(3, 2) = -3; A(3, 3) =  5;

    float * f = new float[4]{6, 1, 0, 1};

    Chasing(A, f, 4);
//the correct answer is : [5, 4, 3, 2]
    delete[] f;
}
#endif

#ifdef _VGM_Test_

#include "Soil.h"

void main()
{
    std::list<float> compart(6, 20.0f);
	std::list<float> thetas{ 0.360f, 0.370f,  0.376f, 0.375f, 0.371f, 0.370f };
	std::list<float> thetar{ 0.044f, 0.044f, 0.050f, 0.050f, 0.039f, 0.029f };
	std::list<float> alpha{ 0.007f, 0.007f, 0.006f, 0.007f, 0.009f, 0.043f };
	std::list<float> n{ 1.59f, 1.40f, 1.48f, 1.45f, 1.32f, 1.30f };
	std::list<float> ksat{ 20.9f, 24.6f, 25.9f, 17.0f, 25.3f, 32.2f };
	std::list<float> hs(6, -2.0f);

    Soil* soil = new VGM(std::array<std::list<float>, 7>{compart, thetas, thetar, alpha, n, ksat, hs});

    float* theta = new float[120];
    float* capacity = new float[120];
    float* conductivity = new float[120];
    float* potential = new float[120];
    float* h = new float[120];

    for (long i = 0; i < 120; ++i)
    {
        h[i] = -i * 10.0 - 2.0;
    }

    soil->theta(h, theta);
    soil->capacity(h, theta, capacity);
    soil->conductivity(theta, conductivity);
    soil->potential(theta, potential);

    delete soil;
    delete[] theta;
    delete[] conductivity;
    delete[] capacity;
    delete[] potential;
    delete[] h;
}
#endif

#ifdef _Solver_Test_

#include "Solver.h"
#include <iostream>
#include <fstream>

void main()
{
    std::list<float> compart(6, 200.0f);
    std::list<float> thetas{ 0.360f, 0.370f,  0.376f, 0.375f, 0.371f, 0.370f };
    std::list<float> thetar{ 0.044f, 0.044f, 0.050f, 0.050f, 0.039f, 0.029f };
    std::list<float> alpha{ 0.007f, 0.007f, 0.006f, 0.007f, 0.009f, 0.043f };
    std::list<float> n{ 1.59f, 1.40f, 1.48f, 1.45f, 1.32f, 1.30f };
    std::list<float> ksat{ 20.9f, 24.6f, 25.9f, 17.0f, 25.3f, 32.2f };
    std::list<float> hs(6, -2.0f);

    Soil* soil = new VGM(std::array<std::list<float>, 7>{compart, thetas, thetar, alpha, n, ksat, hs});

    Boundary* upper = new FixUpper(-2);
    Boundary* lower = new FixLower(-3500, soil->size());

    Sink* sink = new NoSink();

    float* h = new float[soil->size()];
    float* theta = new float[soil->size()];

    for (long i = 0; i < soil->size(); ++i)
    {
        h[i] = -3500.0f;
    }

    h[0] = -2.0f;
    soil->theta(h, theta);

    Solver solver(soil);

    std::ofstream THETA("theta.txt");
    std::ofstream H("potential.txt");

    //Adaptive time control
    float  t = 0;
    while (t < 6.0)
    {
        float dt = fminf(1.0f, 6.0f - t);
        
        while (!solver.Solve(h, theta, dt, sink, upper, lower))
        {
            dt *= 0.5;
        }

        soil->theta(h, theta);

        for (long i = 0; i < soil->size(); ++i)
        {
            THETA << theta[i] << "\t";
            H << h[i] << "\t";
        }

        THETA << std::endl;
        H << std::endl;

        t += dt;
    }

    H.close();
    THETA.close();

    delete soil;
    delete upper;
    delete lower;
    delete sink;

    delete[] theta;
    delete[] h;
}
#endif