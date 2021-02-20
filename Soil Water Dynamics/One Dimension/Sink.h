#ifndef _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SINK_
#define _HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SINK_
class Sink
{
public:
	virtual void Modify(float* h, float* sink, const long& length) = 0;
};

class NoSink : public Sink
{
public:
	virtual void Modify(float* h, float* sink, const long& length);
};

#endif // !_HYDROLOGICAL_TOOLS_SOIL_WATER_DYNAMICS_ONE_DIMENSION_SINK_