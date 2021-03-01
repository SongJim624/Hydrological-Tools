#include "Solver.h"

class Solute_Bresler_Picard : public Solute
{
private:
	Soil* soil;
	const long length;

//M is the basic Matrix stores the varaibles that not changes with the ions
//A is the basic Matrix stores the varaibles that not changes when solving the concentration
	Tridiagonal M, A, B;

//length N - 1
	float* N, * flux_middle, * theta_middle, * velocity_middle, * diffusion;
//length N
	float* temp, * f, * velocity;

public:
	Solute_Bresler_Picard(Soil* soil, const std::map<std::string, std::array<float, 4>>& ions);
	~Solute_Bresler_Picard();

	virtual bool Solve(States& Old, States& New, const float& dt,
		Sink* sink, Boundary* upper, Boundary* lower);
};