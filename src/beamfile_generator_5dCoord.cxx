#include<beamfile_generator_5dCoord.h>

#include<cmath>


int antok::beamfileGenerator::fiveDimCoord::_orderDim = 0;

antok::beamfileGenerator::fiveDimCoord::fiveDimCoord()
{
	_coords.resize(5, 0.);
	_eventNumber = -1;
}

antok::beamfileGenerator::fiveDimCoord::fiveDimCoord(double x0, double x1, double x2, double x3, double x4, long eventNumber)
{
	_coords.resize(5, 0.);
	_coords[0] = x0;
	_coords[1] = x1;
	_coords[2] = x2;
	_coords[3] = x3;
	_coords[4] = x4;
	_eventNumber = eventNumber;
}

bool antok::beamfileGenerator::fiveDimCoord::operator<(const antok::beamfileGenerator::fiveDimCoord& rhs) const
{
	return _coords[_orderDim] < rhs._coords[_orderDim];
}

antok::beamfileGenerator::fiveDimCoord& antok::beamfileGenerator::fiveDimCoord::operator+=(const antok::beamfileGenerator::fiveDimCoord& rhs)
{
	for(unsigned int i = 0; i < 5; ++i) {
		_coords[i] += rhs._coords[i];
	}
	return *this;
}

antok::beamfileGenerator::fiveDimCoord& antok::beamfileGenerator::fiveDimCoord::operator-=(const antok::beamfileGenerator::fiveDimCoord& rhs)
{
	for(unsigned int i = 0; i < 5; ++i) {
		_coords[i] -= rhs._coords[i];
	}
	return *this;
}

antok::beamfileGenerator::fiveDimCoord& antok::beamfileGenerator::fiveDimCoord::operator*=(const double& factor)
{
	for(unsigned int i = 0; i < 5; ++i) {
		_coords[i] *= factor;
	}
	return *this;
}

antok::beamfileGenerator::fiveDimCoord& antok::beamfileGenerator::fiveDimCoord::operator/=(const double& factor)
{
	for(unsigned int i = 0; i < 5; ++i) {
		_coords[i] /= factor;
	}
	return *this;
}

const antok::beamfileGenerator::fiveDimCoord antok::beamfileGenerator::fiveDimCoord::operator+(const antok::beamfileGenerator::fiveDimCoord& rhs) {

	return antok::beamfileGenerator::fiveDimCoord(*this) += rhs;

}

const antok::beamfileGenerator::fiveDimCoord antok::beamfileGenerator::fiveDimCoord::operator-(const antok::beamfileGenerator::fiveDimCoord& rhs) {

	return antok::beamfileGenerator::fiveDimCoord(*this) -= rhs;

}

double antok::beamfileGenerator::fiveDimCoord::distance(const fiveDimCoord& point) const
{
	double dist = 0.;
	for(unsigned int i = 0; i < 5; ++i) {
		dist += (point._coords[i] - _coords[i])*(point._coords[i] - _coords[i]);
	}
	return std::sqrt(dist);
}

double antok::beamfileGenerator::fiveDimCoord::distance2(const fiveDimCoord& point) const
{
	double dist = 0.;
	for(unsigned int i = 0; i < 5; ++i) {
		dist += (point._coords[i] - _coords[i])*(point._coords[i] - _coords[i]);
	}
	return dist;
}

std::ostream& antok::beamfileGenerator::fiveDimCoord::print(std::ostream& out) const
{
	out << "five Dim Coord: [" << _coords[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _coords[i];
	}
	out << "] (event number " << _eventNumber << ")" << std::endl;
	return out;
}
