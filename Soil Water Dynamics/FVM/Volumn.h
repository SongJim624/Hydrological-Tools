#include <list>
#include <array>

class Matrix
{

};

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

	Node center;
	std::list<Volumn*> contiguities;
public:
	Volumn(const Node&);
};

//Mesh is a dense matrix
class Mesh
{
private:
	Node* nodes;
	const long length;

public:
	Mesh();
	~Mesh();

	Node& operator()(const long& x, const long& n y, const long& z);
	Volumn& operator()(const long& x, const long& n y, const long& z);
};