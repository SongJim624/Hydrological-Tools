#include "Body.h"
/*
template<typename T>
Body<T>::Body()
	: depend(Dependency<T>()),
	h(std::vector<T>(0)), theta(std::vector<T>(0)), flux(std::vector<T>(0))
{
	std::vector<T> compart(6, 20.0f);
	std::vector<T> thetas{ 0.360f, 0.370f,  0.376f, 0.375f, 0.371f, 0.370f };
	std::vector<T> thetar{ 0.044f, 0.044f, 0.050f, 0.050f, 0.039f, 0.029f };
	std::vector<T> alpha{ 0.007f, 0.007f, 0.006f, 0.007f, 0.009f, 0.043f };
	std::vector<T> n{ 1.59f, 1.40f, 1.48f, 1.45f, 1.32f, 1.30f };
	std::vector<T> ksat{ 20.9f, 24.6f, 25.9f, 17.0f, 25.3f, 32.2f };
	std::vector<T> hs(6, -2.0f);

	h = std::vector<T>(120, -100);
	h[0] = -2;
	flux = std::vector<T>(h.size());
	depend.ETo = 5;
	depend.Pre = 50;
	depend.Irr = 0;
	depend.CN = 71;
	Kex = 1.10;
	hc = -22.5;
	hcc = -1200;

	depend.sink = std::vector<T>(h.size(), 0.00025);
	sink = depend.sink;
	sink[0] = 0;

	soil = new Soil<T>(compart, thetas, thetar, alpha, n, ksat, hs);
	theta = soil->theta(h);
	state = new State<T>(theta);
}

template<typename T>
Body<T>::~Body()
{
	if (soil) { delete soil; soil = nullptr; }
	if (state) { delete state; state = nullptr; }
}

template<typename T>
T Transpiration(const std::vector<T>& sink, const T& dt)
{
	T trp = 0;

	for (size_t i = 0; i < sink.size() - 1; ++i)
	{
		trp += 0.5 * (sink[i] + sink[i + 1]);
	}

	return trp * dt;
}

template<typename T>
void State<T>::Item(const std::vector<T>& sink, const std::vector<T>& h, const std::vector<T>& theta0,
	const BoundaryCondition<T>& upper, const BoundaryCondition<T>& lower, Soil<T>& soil, const T& dt)
{
	std::vector<T> theta = soil.theta(h);
	std::vector<T> ks = soil.conductivity(theta);

	switch (upper.type)
	{
	case Type::first:
	{
		T flux = -0.5 * (ks[0] + ks[1]) * (h[1] - h[0] - 1) - 0.5 * (theta[0] - theta0[0]) / dt - 0.5 * sink[0];
		infil += flux * dt;
		break;
	}
	case Type::second:
	{
		upper.value > 0 ? infil += upper.value * dt : evp -= upper.value * dt;
		break;
	}
	}

	switch (lower.type)
	{
	case Type::first:
	{
		size_t num = ks.size() - 1;
		T flux = -0.5 * (ks[num] + ks[num - 1]) * (h[num] - h[num - 1] - 1);
		flux > 0 ? dr += flux * dt : cr -= flux * dt;
		break;
	}
	case Type::second:
	{
		lower.value > 0 ? dr += lower.value * dt : cr += lower.value * dt;
		break;
	}
	case Type::third:
	{//Only the condition of free drainage is considered
		size_t num = ks.size() - 1;
		dr += ks[num] * dt;
		break;
	}
	}

	trp += Transpiration(sink, dt);
}

template<typename T>
void State<T>::WaterBalance(const std::vector<T>& theta)
{
	T temp = 0;
	for (size_t i = 0; i < theta.size() - 1; ++i)
	{
		temp += 0.5 * (theta[i] + theta[i + 1]);
	}

	err = 10 * (temp - (wc - evp - trp - dr + infil + cr));
	wc = temp;
}

template<typename T>
void State<T>::Zero() { evp = trp = infil = dr = cr = 0; }

template<typename T>
std::vector<T> Flux(const std::vector<T>& h, const std::vector<T>& theta0, const std::vector<T>& sink, Soil<T>& soil, const T& dt)
{
	std::vector<T> flux(h.size());

	std::vector<T> theta = soil.theta(h);
	std::vector<T> ks = soil.conductivity(theta);

	flux[0] = -0.5 * (ks[0] + ks[1]) * (h[1] - h[0] - 1) - 0.5 * (theta[0] - theta0[0])/ dt - 0.5 * sink[0];

	for (size_t i = 1; i < h.size() - 1; ++i)
	{
		flux[i] = 0.5 * (
			-0.5 * (ks[i] + ks[i + 1]) * (h[i + 1] - h[i] - 1)
			- 0.5 * (ks[i] + ks[i - 1]) * (h[i] - h[i - 1] - 1));
	}

	size_t num = flux.size() - 1;
	flux[num] = -0.5 * (ks[num] + ks[num - 1]) * (h[num] - h[num - 1] - 1);

	return flux;
}

template<typename T>
T SCS(const T& Pre, const unsigned long& CN)
{
	T S = 254.0 * (100.0 / CN - 1);
	T Ia = 0.05 * S;
	//The initial abstraction is fixed at 0.05S
	//There may be some better choice in the future

//Return the run off
	return Pre > Ia ? pow((Pre - Ia), 2) / (Pre + S - Ia) : 0.0f;
}

template<typename T>
void Body<T>::Update()
{
	T time = 0;
	T dt = 0.01;

	state->Zero();
	state->Run() = SCS(depend.Pre, depend.CN);
	infil = 0.1 * (depend.Pre + depend.Irr - state->Run());
	trp = depend.trp;

	while (time < 24)
	{
		Solver(h, theta, dt, *soil, depend);
	// flux is needed to decide wh
		flux = Flux(h, theta, depend.sink, *soil, dt);

		state->Item(sink, h, theta, upper, lower, *soil, dt);
		theta = soil->theta(h);

		//state->WaterBalance(theta);
		time += dt;
	}

	state->Run() += infil - state->Infil();
	state->WaterBalance(theta);
}
*/