#include "Soil.h"

VGM::VGM(const std::array<std::list<float>, 7>& soil)
{
	length = 0;
	
	for (auto it : *soil.begin())
	{
		length += it;
	}

	m = new float[length];
	thetam = new float[length];
	inter = new float[length];
	F1 = new float[length];

	label = new bool[length];

	Convert(soil[0], soil[1], thetas);
	Convert(soil[0], soil[2], thetar);
	Convert(soil[0], soil[3], alpha);
	Convert(soil[0], soil[4], n);
	Convert(soil[0], soil[5], ksat);
	Convert(soil[0], soil[6], hs);

	for (size_t i = 0; i < length; ++i)
	{
		inter[i] = pow(-alpha[i] * hs[i], n[i]);
		m[i] = 1 - 1 / n[i];
		thetam[i] = (thetas[i] - thetar[i]) * pow(1 + inter[i], m[i]) + thetar[i];
		F1[i] = pow(inter[i] / (1 + inter[i]), m[i]);
	}
}

VGM::~VGM()
{
	delete[] thetas;
	delete[] thetar;
	delete[] alpha;
	delete[] n;
	delete[] ksat;
	delete[] hs;

	delete[] m;
	delete[] thetam;
	delete[] inter;
	delete[] F1;

	delete[] label;
}

void VGM::Convert(const std::list<float>& compart, const std::list<float>& value, float*& res)
{
	std::vector<float> r(0);

	for (size_t i = 0; i < compart.size(); ++i)
	{
		std::vector<float> temp(*std::next(compart.begin(), i), *std::next(value.begin(), i));
		r.insert(r.end(), temp.begin(), temp.end());
	}

	res = new float[r.size()];
	
	for (long i = 0; i < r.size(); ++i)
	{
		res[i] = r[i];
	}
}

void VGM::theta(float * h, float * res)
{
	for (long i = 0; i < length; ++i)
	{
		label[i] = h[i] < hs[i] ? true : false;
		inter[i] = label[i] ? pow(-alpha[i] * h[i], n[i]) : 0;
		res[i] = label[i] ? (thetam[i] - thetar[i]) * pow(1 + inter[i], -m[i]) + thetar[i] : thetas[i];
	}
}

void VGM::capacity(float* h, float* theta, float* res)
{
	for (long i = 0; i < length; ++i)
	{
		res[i] = label[i] ? (1 - n[i]) * (theta[i] - thetar[i]) * inter[i] / (inter[i] + 1) / h[i] : 0;
	}
}

void VGM::conductivity(float* theta, float* res)
{
	for (long i = 0; i < length; ++i)
	{
		res[i] = label[i] ? ksat[i] * sqrt((theta[i] - thetar[i]) / (thetas[i] - thetar[i])) *
			pow((1 - pow(inter[i] / (1 + inter[i]), m[i])) / (1 - F1[i]), 2) : ksat[i];
	}
}

void VGM::potential(float* theta, float * h)
{
	float* th = new float[length];
	float* cap = new float[length];

	for (long i = 0; i < length; ++i)
	{
		h[i] = -100;
	}

	while (true)
	{
		this->theta(h, th);
		this->capacity(h, th, cap);

		float norm = -INFINITY;

		for (long i = 0; i < length; ++i)
		{
			th[i] = (th[i] - theta[i]) / cap[i];
			norm = fmax(norm, abs(th[i]));

			h[i] -= th[i];
		}

		if (norm < 1e-3)
		{
			break;
		}
	}

	delete[] th;
	delete[] cap;
}

long VGM::size()
{
	return length;
}