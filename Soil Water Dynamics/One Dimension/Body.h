#pragma once

#include "Soil.h"
/*
enum class Type
{
	first = 1, second = 2, third = 3
};

template<typename T>
struct BoundaryCondition
{
	Type type;
	T value;
};

template<typename T>
struct Dependency
{
	T ETo, Pre;
	unsigned long CN;
	T Irr, reduce_percnet, irr_percent;
	T CC, trp;
	std::vector<T> sink;
};

template<typename T>
class State
{
private:
	T evp, trp, pre, irr, infil, run, dr, cr, wc, err;
//	const std::vector<T>& theta;

public:
	State(const std::vector<T>& theta) : evp(0), trp(0), pre(0), irr(0), infil(0), run(0), dr(0), cr(0), wc(0), err(0)
	{
		for (size_t i = 0; i < theta.size() - 1; ++i)
		{
			wc += 0.5 * (theta[i] + theta[i + 1]);
		}
	}

	void Zero();
	void Item(const std::vector<T>& sink, const std::vector<T>& h, const std::vector<T>& theta0, 
		const BoundaryCondition<T>& upper, const BoundaryCondition<T>& lower, Soil<T>& soil, const T& dt);
	void WaterBalance(const std::vector<T>& theta);
	T& Infil() { return infil; }
	T& Trp() { return trp; }
	T& Run() { return run; }
};

template<typename T>
class Body
{
	Soil<T> *soil;
	State<T> *state;
	Dependency<T> depend;
	std::vector<T> h, theta, flux;

	BoundaryCondition<T> upper, lower;


These three variables are not the same the variables in the State Class
These are to record the process of infiltration and tanspiration
and would change during the update process
the sink in this class presents the sink speed
the sink in the dependency class represents the sink ability, the max sink speed.

	T infil, trp;
	std::vector<T> sink;

	T Kex, hc, hcc;

private:
//	void Lower();
//	void Upper(const dependency<T> depend);

	std::vector<T> Solver(std::vector<T>& h, const std::vector<T>& theta0, const T& dt, Soil<T>& soil,
		const Dependency<T>& depend);

public:
	Body();
	~Body();

public:
	void Update();

	void Query();
	void Description();
	void SetUp();
	void Convert();
	void Report();
	void Clock();
};
*/