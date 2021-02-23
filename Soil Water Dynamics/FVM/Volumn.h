#include <list>
#include <array>

struct Node
{
	float x, y, z;
};

class Volumn
{
private:
	
	float theta, temperature;	
	std::list<float> solutes;

// tetrahedrold Volumn
	std::array<Node, 4> nodes;

	std::list<Volumn*> conjunctions;
public:

};