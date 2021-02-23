#include "Solute.h"
#include "Tridiagonal.h"
#include "Soil.h"

class Solution
{
private:
	const long length;
	Soil* soil;
	Tridiagonal A, B; 

	float * theta, * flux, * velocity, * Dsh;
	float * thetan, * fluxn, * velocityn, * Dshn;

	float * conductivity;

	void Flux(float * h, const float& dt, float * conductivity);

public:
	Solution(const long& legnth);
	void Solve(Solute& solution, float* flux, float* theta, float* theta0, const float& dt);
};