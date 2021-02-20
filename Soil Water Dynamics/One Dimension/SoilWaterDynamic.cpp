// SoilWaterDynamic.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
/*
#include <iostream>
#include "Soil.h"
#include "Soil.cpp"
int main()
{    
    std::list<std::list<float>> soil(0);
    
    std::list<float> compart(6, 20.0f);
    std::list<float> thetas{ 0.360f, 0.370f,  0.376f, 0.375f, 0.371f, 0.370f };
    std::list<float> thetar{ 0.044f, 0.044f, 0.050f, 0.050f, 0.039f, 0.029f };
    std::list<float> alpha{ 0.007f, 0.007f, 0.006f, 0.007f, 0.009f, 0.043f };
    std::list<float> n{ 1.59f, 1.40f, 1.48f, 1.45f, 1.32f, 1.30f};
    std::list<float> ksat{ 20.9f, 24.6f, 25.9f, 17.0f, 25.3f, 32.2f };
    std::list<float> hs(6, -2.0f);

    soil.push_back(compart);
    soil.push_back(thetas);
    soil.push_back(thetar);
    soil.push_back(alpha);
    soil.push_back(n);
    soil.push_back(ksat);
    soil.push_back(hs);

    Soil s(soil);


    std::list<float> h(120, -100);
    h[0] = -2;

    float dt = 0.01;

    Soil<float> soil(compart, thetas, thetar, alpha, n, ksat, hs);

    Body<float>* body = new Body<float>();

    body->Update();

    delete body;

    double T = 0;
    while (T < 24)
    {
        std::list<float> theta = Solver(h, theta0, dt, soil, depend);
        theta0 = theta;
        T += dt;
    }

    std::cout << "Hello World!\n";
}
*/