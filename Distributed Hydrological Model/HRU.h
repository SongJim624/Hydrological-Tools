#include <list>

struct Coordinate
{
    float x, y, z;
};

struct Climate
{
    float Tmax, Tmin, Pre, U2, Rad;
};

class Soil
{

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
private:
    float Area;
    std::list<Coordinate> boundary;
//Components    
private:
    Climate climate;
    Soil soil;

    
    
    
//States
private:

public:

};
