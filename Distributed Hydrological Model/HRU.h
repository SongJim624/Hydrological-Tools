#include <list>

struct Coordinate
{
    float x, y, z;
};

class Soil
{

};

class Climate
{
private:
    std::list<float> Tmaxs, Tmins, ETos, Rads;
};

class Solute
{
private:
    const long length;

public:    
    float * concentrations;
};

class HRU
{
//Components    
private:
    Climate climate;
    Soil soil;

//States
private:
    float* theta;
    std::list<Solute> solutes;
    float snow;

public:

};