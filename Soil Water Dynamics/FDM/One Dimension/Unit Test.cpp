//#define _Tridiagonal_Test_
#define _Solver_Test_

#ifdef _Tridiagonal_Test_

#include <vector>
#include "Tridiagonal.h"

void main()
{
    Tridiagonal A(4);

    A(0, 0) =  2; A(0, 1) = -1;
    A(1, 0) = -1; A(1, 1) =  3; A(1, 2) = -2;
                  A(2, 1) = -1; A(2, 2) =  2; A(2, 3) = -1;
                                A(3, 2) = -3; A(3, 3) =  5;

    float * f = new float[4]{6, 1, 0, 1};
    float* ans = new float[4]{ 5, 4, 3, 2 };
    float* res = new float[4];

    Chasing(A, f, 4);
    A.multiple(ans, res);

//the correct answer is : [5, 4, 3, 2]
    delete[] f;
    delete[] ans;
    delete[] res;
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
	std::list<float> ksat{ 20.9f, 24.6f, 25.9f, 17.0f, 25.3f, 32.2f }; //ks might be a little big
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

#include "Modified Van Genutchen.h"
#include "Solver.h"
#include "Water-Celia-Picard.h"
#include "Solute-Bresler-Picard.h"
#include <iostream>
#include <fstream>

void main()
{
/*
    std::list<float> compart(6, 20.0f);
    std::list<float> thetas{ 0.360f, 0.370f,  0.376f, 0.375f, 0.371f, 0.370f };
    std::list<float> thetar{ 0.044f, 0.044f, 0.050f, 0.050f, 0.039f, 0.029f };
    std::list<float> alpha{ 0.007f, 0.007f, 0.006f, 0.007f, 0.009f, 0.043f };
    std::list<float> n{ 1.59f, 1.40f, 1.48f, 1.45f, 1.32f, 1.30f };
    std::list<float> ksat{ 20.9f, 24.6f, 25.9f, 17.0f, 25.3f, 32.2f };
    std::list<float> hs(6, -2.0f);
*/
    std::list<float> compart{ 100.0 };
    std::list<float> thetas{ 0.41 };
    std::list<float> thetar{ 0.065 };
    std::list<float> alpha{ 0.075 };
    std::list<float> n{ 1.89 };
    std::list<float> ksat{ 4.34 };
    std::list<float> hs{ -2 };


    Soil* soil = new VGM(std::array<std::list<float>, 7>{compart, thetas, thetar, alpha, n, ksat, hs});

    Boundary* upper = new First(0, -2.0f);
    Boundary* lower = new First(soil->size(), -350.0f);

    Boundary* upper_salt = new First(0, 100.0f);
    Boundary* lower_salt = new First(soil->size(), 0.0f);

    Sink* sink = new Sink();

    States::length = soil->size();

    States Old = States();

    for (long i = 0; i < soil->size(); ++i)
    {
        Old.potential[i] = -350.0f;
    }

    soil->theta(Old.potential, Old.theta);
    Old.concentrations["Na"] = (float*)MKL_calloc(soil->size(), sizeof(float), 64);

    std::map<std::string, std::array<float, 4>> ions;
    ions["Na"] = std::array<float, 4>{0.55, 2.4, 0.28, 0.39};


    Water* solver = new Water_Celia_Picard(soil);
    Solute* solver_solute = new Solute_Bresler_Picard(soil, ions);


    solver->Flux(Old, Old, 0.0f);
    States New = States(Old);

    std::ofstream THETA("theta.txt");
    std::ofstream H("potential.txt");
    std::ofstream FLUX("flux.txt");
    std::ofstream SALT("salt.txt");

    //Adaptive time control

    float  t = 0;
    while (t < 6.0)
    {
        float dt = fminf(1.0f, 6.0f - t);
        
        while (!solver->Solve(Old, New, dt, sink, upper, lower))
        {
            dt *= 0.5;
        }

        soil->theta(New.potential, New.theta);

        solver->Flux(Old, New, dt);

        solver_solute->Solve(Old, New, dt, sink, upper_salt, lower_salt);

        for (long i = 0; i < soil->size(); ++i)
        {
            THETA << New.theta[i] << "\t";
            FLUX << New.flux[i] << "\t";
            SALT << New.concentrations["Na"][i] << "\t";
        }

        THETA << std::endl;
        FLUX << std::endl;
        SALT << std::endl;

        cblas_scopy(soil->size(), New.potential, 1, Old.potential, 1);
        cblas_scopy(soil->size(), New.theta, 1, Old.theta, 1);
        cblas_scopy(soil->size(), New.concentrations["Na"], 1, Old.concentrations["Na"], 1);

        t += dt;
    }

    H.close();
    THETA.close();
    FLUX.close();
    SALT.close();
    delete soil;
    delete upper;
    delete lower;
    delete sink;
    delete solver;
    delete solver_solute;
}
#endif