#include "Solution.h"

Solution::Solution(const long& length) 
	: length(length), A(Tridiagonal(length)), B(Tridiagonal(length))
{}

void Solution::Flux(float * h, float * conductivity, float * flux)
{
	flux[0] = -0.5 * (conductivity[0] + conductivity[1]) * (h[1] - h[0] - 1) - 0.5 * (theta[0] - theta0)/ dt - 0.5 * s[0];

	for (long i = 1; i < length - 1; ++i)
	{
		flux[i] = 0.5 * (
			-0.5 * (conductivity[i] + conductivity[i + 1]) * (h[i + 1] - h[i] - 1)
			- 0.5 * (conductivity[i] + conductivity[i - 1]) * (h[i] - h[i - 1] - 1));
	}

	flux[length - 1] = -0.5 * (conductivity[length - 1] + conductivity[length - 2]) * (h[length - 1] - h[length - 2] - 1);
}

void Solution::Solve(Solute& solution, float* h, float* hn, const float& dt)
{
	soil->theta(h, theta);
	soil->conductivity(theta, conductivity);
	solution.difussion(velocity, theta, Dsh);
	
	soil->theta(hn, thetan);
	soil->conductivity(thetan, conductivityn);
	solution.difussion(velocityn, thetan, Dshn);

	for (long i = 1; i < length - 1; ++i)
	{
		float qa = 0.5 * (
			0.5 * (conductivity[i - 1] + conductivity[i]) * (h[i] - h[i - 1] + 1) + 
			0.5 * (conductivityn[i - 1] + conductivityn[i]) * (hn[i] - hn[i - 1] + 1));

		float qb = 0.5 * (
			0.5 * (conductivity[i] + conductivity[i + 1]) * (h[i + 1] - h[i] + 1) + 
			0.5 * (conductivityn[i] + conductivityn[i + 1]) * (hn[i + 1] - hn[i] + 1));

		float q = 0.5 * (
			0.5 * conductivity[i] * (h[i + 1] - 2 * h[i] + h[i - 1] + 2) +
			0.5 * conductivityn[i] * (hn[i + 1] - 2 * hn[i] + hn[i - 1] + 2)			
		);

		float va = qa / 0.25 * (theta[i - 1] + theta[i] + thetan[i - 1] + thetan[i]);
		float vb = qb / 0.25 * (theta[i] + theta[i + 1] + thetan[i] + thetan[i + 1]);

		float Na = 0.5 * qa - 0.125 * 0.5 * dt * q[i] / (theta[i] + thetan[i]) * va * (thetan[i] - theta[i]);
		float Nb = 0.5 * qb - 0.125 * 0.5 * dt * q[i] / (theta[i] + thetan[i]) * vb * (thetan[i] - theta[i]);

		A(i, i - 1) = Na - 0.5 * (Dsh[i - 1] + Dsh[i]) - qa;
		A(i, i) = theta[i] / r + Dsh[i] - Na - Nb;
		A(i, i + 1) = Nb - 0.5 * (Dsh[i] + Dsh[i + 1]);

		B(i, i - 1) = -A(i, i - 1);
		B(i, i) = c[i + 1] - (flux[] + Dsh[i] - N - N - theta0[i] / r)
		B(i, i + 1) = -A(i, i + 1);
	}

	B.multiple(solution.concentration, f);
	Chasing(A, f, length);

	if(Norm(f, solution.concentration, length) < 0.01)
	{
		return;
	}

	for(long i = 0; i < length; ++i)
	{
		solution.concentration[i] = f[i];
	}
}