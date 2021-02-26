#include "Solute.h"
#include "Tridiagonal.h"
#include "Soil.h"

struct Water
{
	float* theta;
	float* flux;
};

class Solution
{
private:
	const long length;
	Soil* soil;
	Tridiagonal A, B; 

	float * theta, * flux, * velocity, * Dsh;
	float * thetan, * fluxn, * velocityn, * Dshn;

	float * conductivity;

private:
	void Velocity(float* flux, float * theta, float* velocity);
	void Flux(float * h, const float& dt, float * conductivity);

public:
	Solution(const long& legnth);

	void Solve(Solute& solution, Water& water_o, Water& water_n, const float& dt);
};