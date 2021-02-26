#include "Modified Van Genutchen.h"

VGM::VGM(const std::array<std::list<float>, 7>& soil)
{
	length = 0;

	for (auto it : *soil.begin())
	{
		length += it;
	}

	m = (float*)MKL_calloc(length, sizeof(float), 64);

	inter = (float*)MKL_calloc(length, sizeof(float), 64);
	F1 = (float*)MKL_calloc(length, sizeof(float), 64);
	thetam = (float*)MKL_calloc(length, sizeof(float), 64);

	positive_ones = (float*)MKL_calloc(length, sizeof(float), 64);
	negative_ones = (float*)MKL_calloc(length, sizeof(float), 64);
	temp = (float*)MKL_calloc(length, sizeof(float), 64);
	
	for (long i = 0; i < length; ++i)
	{
		positive_ones[i] = 1;
		negative_ones[i] = -1;
	}

	label = new bool[length];

	Convert(soil[0], soil[1], thetas);
	Convert(soil[0], soil[2], thetar);
	Convert(soil[0], soil[3], alpha);
	Convert(soil[0], soil[4], n);
	Convert(soil[0], soil[5], ksat);
	Convert(soil[0], soil[6], hs);

// inter
	vsMul(length, alpha, hs, inter);
	vsMul(length, inter, negative_ones, inter);
	vsPow(length, inter, n, inter);
	
// m
	vsDiv(length, negative_ones, n, m);
	vsAdd(length, positive_ones, m, m);

// thetam
	vsAdd(length, positive_ones, inter, thetam);
	vsPow(length, thetam, m, thetam);
	vsSub(length, thetas, thetar, temp);
	vsMul(length, thetam, temp, thetam);
	vsAdd(length, thetam, thetar, thetam);

//F1
	vsAdd(length, positive_ones, inter, F1);
	vsDiv(length, inter, F1, F1);
	vsPow(length, F1, m, F1);
	vsSub(length, positive_ones, F1, F1);
}

VGM::~VGM()
{
	MKL_free(thetas);
	MKL_free(thetar);
	MKL_free(alpha);
	MKL_free(n);
	MKL_free(ksat);
	MKL_free(hs);

	MKL_free(m);
	MKL_free(thetam);
	MKL_free(inter);
	MKL_free(F1);

	MKL_free(positive_ones);
	MKL_free(negative_ones);
	MKL_free(temp);

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

	res = (float*) MKL_calloc(r.size(), sizeof(float), 64);
	
	for (long i = 0; i < r.size(); ++i)
	{
		res[i] = r[i];
	}
}

void VGM::theta(float * h, float * res)
{
	vsMul(length, alpha, h, inter);
	vsMul(length, negative_ones, inter, inter);
	vsPow(length, inter, n, inter);

	vsMul(length, negative_ones, m, temp);
	vsAdd(length, positive_ones, inter, res);
	vsPow(length, res, temp, res);

	vsSub(length, thetam, thetar, temp);
	vsMul(length, temp, res, res);
	vsAdd(length, thetar, res, res);

	for (long i = 0; i < length; ++i)
	{
		label[i] = h[i] < hs[i] ? true : false;

		if (!label[i])
		{
			res[i] = thetas[i];
		}
	}
}

void VGM::capacity(float* h, float* theta, float* res)
{
	vsSub(length, positive_ones, n, temp);
	vsSub(length, theta, thetar, res);
	vsMul(length, res, temp, res);
	vsMul(length, res, inter, res);
	vsAdd(length, inter, positive_ones, temp);
	vsDiv(length, res, temp, res);
	vsDiv(length, res, h, res);

	for (long i = 0; i < length; ++i)
	{
		if (!label[i])
		{
			res[i] = 0;
		}
	}
}

void VGM::conductivity(float* theta, float* res)
{
	vsSub(length, theta, thetar, res);
	vsSub(length, thetas, thetar, temp);
	vsDiv(length, res, temp, res);
	vsSqr(length, res, res);
	vsMul(length, res, ksat, res);

	vsAdd(length, positive_ones, inter, temp);
	vsDiv(length, inter, temp, temp);
	vsPow(length, temp, m, temp);
	vsSub(length, positive_ones, temp, temp);
	vsDiv(length, temp, F1, temp);
	vsMul(length, temp, temp, temp);
	vsMul(length, res, temp, res);

	for (long i = 0; i < length; ++i)
	{
		if (!label[i])
		{
			res[i] = ksat[i];
		}
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