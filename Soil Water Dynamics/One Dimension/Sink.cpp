#include "Sink.h"

void Sink::Modify(float* h, float* sink, const long& length)
{
	for (long i = 0; i < length; ++i)
	{
		sink[i] = 0.0;
	}
}