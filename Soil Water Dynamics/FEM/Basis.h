class Basis
{

};

class Linear : public Basis
{
public:
	virtual float operator()(const float& x, const float& y, const float& z);
};

class Hermite : public Basis
{

};