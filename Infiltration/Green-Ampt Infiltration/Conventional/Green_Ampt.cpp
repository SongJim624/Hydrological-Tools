#include <math>

/*
* Parameters:
* layer     : the water layer on the surface of the soil
* thetas    : the saturation of the soil
* theta0    : the initial water content of the soil
* Ksat      : the hydraulic conductivity of the soil
* sf        : shui xi li
* z         : the depth of the wet front 
*
* Equation:
* t = (thetas - theta0) / K * 
* (z - (sf + H） * ln((z + sf + H) / (sf + H ))
*
*/

inline double F(const double& z, const double& A, const double& B)
{
    return z - A * log(z / A + 1) - B;
}

inline double dF(const double&z, const double& A)
{
    return A / (z + A);
}

//Newton method is used to solve the Green-Ampt equation
double Green_Ampt(const double& Ks, const double& sf,
const double& H, const double& thetas, const double& theta0,
const double& t)
{
    double z = 0.0;

    double A = (sf + H);
    double B = (thetas - theta0) / Ks;

    while(true)
    {
        double norm = F(z, A, B) / dF(z, A);

        if(abs(norm) < 1e-6)
        {
            return z;
        }

        z -= norm;
    }
}