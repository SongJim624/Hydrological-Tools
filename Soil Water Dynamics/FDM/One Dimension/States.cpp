#include "States.h"

long States::length = 0;

States::States()
	: potential((float*)MKL_calloc(length, sizeof(float), 64)),
	theta((float*)MKL_calloc(length, sizeof(float), 64)),
	temperature((float*)MKL_calloc(length, sizeof(float), 64)),
	flux((float*)MKL_calloc(length, sizeof(float), 64)),
	flux_middle((float*)MKL_calloc(length - 1, sizeof(float), 64))
{
}

States::States(const States& states)
	:potential((float *)MKL_calloc(length, sizeof(float), 64)),
	theta((float *)MKL_calloc(length, sizeof(float), 64)),
	temperature((float *)MKL_calloc(length, sizeof(float), 64)),
	flux((float *)MKL_calloc(length, sizeof(float), 64)),
	flux_middle((float *)MKL_calloc(length - 1, sizeof(float), 64))
{
	cblas_scopy(length, states.potential, 1, potential, 1);
	cblas_scopy(length, states.theta, 1, theta, 1);
	cblas_scopy(length, states.temperature, 1, temperature, 1);
	cblas_scopy(length, states.flux, 1, flux, 1);
	cblas_scopy(length - 1, states.flux_middle, 1, flux_middle, 1);

	for (auto iter : states.concentrations)
	{
		concentrations[iter.first] = (float*)MKL_calloc(length, sizeof(float), 64);
		cblas_scopy(length, iter.second, 1, concentrations[iter.first], 1);
	}
}

States::~States()
{
	MKL_free(potential);
	MKL_free(theta);
	MKL_free(temperature);
	MKL_free(flux);
	MKL_free(flux_middle);

	for (auto iter : concentrations)
	{
		MKL_free(iter.second);
	}
}