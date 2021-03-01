#include "Solute-Bresler-Picard.h"

Solute_Bresler_Picard::Solute_Bresler_Picard(Soil* soil, const std::map<std::string, std::array<float, 4>>& ions)
	: soil(soil), length(soil->size()),
	M(Tridiagonal(length)), A(Tridiagonal(length)), B(Tridiagonal(length)),
	N((float*)MKL_calloc(length - 1, sizeof(float), 64)),
	flux_middle((float*)MKL_calloc(length - 1, sizeof(float), 64)),
	theta_middle((float*)MKL_calloc(length - 1, sizeof(float), 64)),
	velocity_middle((float*)MKL_calloc(length - 1, sizeof(float), 64)),
	diffusion((float*)MKL_calloc(length - 1, sizeof(float), 64)),
	temp((float*)MKL_calloc(length, sizeof(float), 64)),
	f((float*)MKL_calloc(length, sizeof(float), 64)),
	velocity((float*)MKL_calloc(length, sizeof(float), 64))
{
	for (auto iter : ions)
	{
		Solute::ions[iter.first] = new Ion(iter.second);
	}
}

bool Solute_Bresler_Picard::Solve(States& Old, States& New, const float& dt,
	Sink* sink, Boundary* upper, Boundary* lower)
{
/*
//flux of the middle node at t + 1/2 
	cblas_scopy(length - 1, New.flux_middle, 1, flux_middle, 1);
	cblas_saxpby(length - 1, 0.5f, Old.flux_middle, 1, 0.5f, flux_middle, 1);

//theta of the middle node at t + 1/2 
	cblas_scopy(length - 1, New.theta, 1, theta_middle, 1);
	cblas_saxpby(length - 1, 0.25f, Old.theta, 1, 0.25f, theta_middle, 1);
	cblas_saxpy(length - 1, 0.25f, Old.theta + 1, 1, theta_middle, 1);
	cblas_saxpy(length - 1, 0.25f, New.theta + 1, 1, theta_middle, 1);

//velocity of node i + 1/2 at t + 1/2 
	vsDiv(length - 1, flux_middle, theta_middle, velocity_middle);

//velocity of the node at i at t + 1/2
	vsDiv(length, Old.flux, Old.theta, velocity);
	vsDiv(length, New.flux, New.theta, temp);
	cblas_saxpby(length, 0.5f, temp, 1, 0.5, velocity, 1);

//N
	vsSub(length - 1, New.theta, Old.theta, N);
	vsMul(length - 1, N, velocity_middle, N);
	vsMul(length - 1, N, velocity, N);
	cblas_saxpby(length - 1, 0.5f, flux_middle, 1, -0.125f * dt, N, 1);
	
//A
	vsSub(length - 2, N, flux_middle, M.Diagonal(-1));
//B	
	vsAdd(length - 2, N, N + 1, M.Diagonal(0) + 1);
	vsSub(length - 2, flux_middle + 1, M.Diagonal(0) + 1, M.Diagonal(0) + 1);
//	cblas_saxpy(length - 2, 0.5 * dt, New.theta + 1, 1, M.Diagonal(0) + 1, 1);
//C
	cblas_scopy(length - 2, N + 1, 1, M.Diagonal(1) + 1, 1);

	for (auto iter : ions)
	{
		float* concentration = nullptr;
		concentration = New.concentrations[iter.first];

		iter.second->diffusion(velocity_middle, theta_middle, diffusion, length - 1);
//Coefficent Matrix of c(t+1)
		vsSub(length - 2, M.Diagonal(-1), diffusion, A.Diagonal(-1));
		
		vsAdd(length - 2, M.Diagonal(0) + 1, diffusion, A.Diagonal(0) + 1);
		vsAdd(length - 2, A.Diagonal(0) + 1, diffusion + 1, A.Diagonal(0) + 1);
		cblas_saxpy(length - 2, 0.5 * dt, New.theta + 1, 1, A.Diagonal(0) + 1, 1);
		
		vsSub(length - 2, M.Diagonal(-1), diffusion + 1, A.Diagonal(-1));
//Coefficent Matrix of c(t)
		cblas_saxpby(length - 2, -1, A.Diagonal(-1), 1, 0, B.Diagonal(-1), 1);

		vsAdd(length - 2, M.Diagonal(0) + 1, diffusion, B.Diagonal(0) + 1);
		vsAdd(length - 2, A.Diagonal(0) + 1, diffusion + 1, B.Diagonal(0) + 1);
		cblas_saxpy(length - 2, 0.5 * dt, Old.theta + 1, -1, B.Diagonal(0) + 1, 1);

		cblas_saxpby(length - 2, -1, A.Diagonal(1), 1, 0, B.Diagonal(1), 1);
*/
	for (long i = 0; i < length - 1; ++i)
	{
		theta_middle[i] = 0.25f * (Old.theta[i] + Old.theta[i + 1] + New.theta[i] + New.theta[i + 1]);
		flux_middle[i] = 0.5f * (Old.flux_middle[i] + New.flux_middle[i]);
		velocity_middle[i] = flux_middle[i] / theta_middle[i];
	}

	for (long i = 0; i < length; ++i)
	{
		velocity[i] = 0.5f * (Old.flux[i] / Old.theta[i] + New.flux[i] / New.theta[i]);
	}

	for (long i = 0; i < length - 1; ++i)
	{
		N[i] = 0.5f * flux_middle[i] - 0.125f * velocity[i] * velocity_middle[i] * (New.theta[i] - Old.theta[i]);
	}

	for (auto iter : ions)
	{
		float* concentration = nullptr;
		concentration = New.concentrations[iter.first];

		iter.second->diffusion(velocity_middle, theta_middle, diffusion, length - 1);

		for (long i = 1; i < length - 1; ++i)
		{
			B(i, i - 1) = N[i] - diffusion[i] - flux_middle[i];
			B(i, i) = 2 * New.theta[i] / dt + diffusion[i] + diffusion[i + 1] - N[i] - N[i + 1] + flux_middle[i + 1];
			B(i, i + 1) = N[i + 1] - diffusion[i + 1];
		}

		while (true)
		{
			for (long i = 1; i < length - 1; ++i)
			{
				f[i] = -B(i, i - 1) * concentration[i - 1] - B(i, i + 1) * concentration[i + 1] - concentration[i] *
					(flux_middle[i] + diffusion[i] + diffusion[i + 1] - N[i] - N[i + 1] - 2.0f * Old.theta[i] / dt);
			}

			B(0, 1) = 0.0;
			B(length - 1, length - 2) = 0.0;

			upper->Modify(B, f, dt);
			lower->Modify(B, f, dt);

			Chasing(B, f, length);

			if (Norm(f, concentration, length) < 0.01)
			{
				break;
			}

			cblas_scopy(length, f, 1, concentration, 1);
		}
	}

	return true;
}

Solute_Bresler_Picard::~Solute_Bresler_Picard()
{
	MKL_free(N);
	MKL_free(temp);
	MKL_free(f);

	MKL_free(flux_middle);
	MKL_free(theta_middle);
	MKL_free(velocity_middle);
	MKL_free(velocity);
	MKL_free(diffusion);

	for (auto iter : Solute::ions)
	{
		delete iter.second;
	}
}
