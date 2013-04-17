#include<beamfile_generator_5dBin.h>

#include<algorithm>
#include<assert.h>
#include<cmath>
#include<limits>


const double antok::beamfileGenerator::fiveDimBin::EPSILON = 5. * std::numeric_limits<double>::epsilon();
long antok::beamfileGenerator::fiveDimBin::nExistingBins = 0;

antok::beamfileGenerator::fiveDimBin::fiveDimBin(double a0,
                                                 double a1,
                                                 double a2,
                                                 double a3,
                                                 double a4,
                                                 double b0,
                                                 double b1,
                                                 double b2,
                                                 double b3,
                                                 double b4)
	: _a(5, 0.),
	  _b(5, 0.),
	  _neighbors(),
	  _onLowerEdge(5, false),
	  _onUpperEdge(5, false),
	  _sigmaCache(0)
{
	_a[0] = a0 < b0 ? a0 : b0;
	_a[1] = a1 < b1 ? a1 : b1;
	_a[2] = a2 < b2 ? a2 : b2;
	_a[3] = a3 < b3 ? a3 : b3;
	_a[4] = a4 < b4 ? a4 : b4;
	_b[0] = a0 < b0 ? b0 : a0;
	_b[1] = a1 < b1 ? b1 : a1;
	_b[2] = a2 < b2 ? b2 : a2;
	_b[3] = a3 < b3 ? b3 : a3;
	_b[4] = a4 < b4 ? b4 : a4;
	nExistingBins += 1;
}

antok::beamfileGenerator::fiveDimBin::fiveDimBin(const std::vector<double>& a, const std::vector<double>& b)
	: _a(5, 0.),
	  _b(5, 0.),
	  _neighbors(),
	  _onLowerEdge(5, false),
	  _onUpperEdge(5, false),
	  _sigmaCache(0)
{
	assert(a.size() == 5);
	assert(b.size() == 5);
	for(unsigned int i = 0; i < 5; ++i) {
		_a[i] = a[i] < b[i] ? a[i] : b[i];
		_b[i] = a[i] < b[i] ? b[i] : a[i];
	}
	nExistingBins += 1;
}

antok::beamfileGenerator::fiveDimBin::~fiveDimBin()
{
	_neighbors.clear();
	nExistingBins -= 1;
	delete _sigmaCache;
}

void antok::beamfileGenerator::fiveDimBin::set(const std::vector<double>& a, const std::vector<double>& b)
{
	assert(a.size() == 5);
	assert(b.size() == 5);
	for(unsigned int i = 0; i < 5; ++i) {
		_a[i] = a[i] < b[i] ? a[i] : b[i];
		_b[i] = a[i] < b[i] ? b[i] : a[i];
	}
}

bool antok::beamfileGenerator::fiveDimBin::inBin(const std::vector<double>& x) const
{
	assert(x.size() == 5);
	for(unsigned int i = 0; i < 5; ++i) {
		if(x[i] < _a[i] or x[i] >= _b[i]) {
			return false;
		}
	}
	return true;
}

bool antok::beamfileGenerator::fiveDimBin::areWeNeighbors(boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin) const
{
	for(unsigned int i = 0; i < 5; ++i) {
		bool neighbors = true;
		if(doubleEqual(bin->getLowerCorner()[i], _b[i]) or doubleEqual(bin->getUpperCorner()[i], _a[i])) {
			for(unsigned int j = 0; j < 5; ++j) {
				if(i == j) {
					continue;
				}
				if(bin->getLowerCorner()[j] > _b[j] or bin->getUpperCorner()[j] < _a[j]) {
					neighbors = false;
					break;
				}
			}
			if(neighbors) {
				return true;
			}
		}
	}
	return false;

}

void antok::beamfileGenerator::fiveDimBin::removeNeighbor(boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin) {
	_neighbors.erase(boost::const_pointer_cast<antok::beamfileGenerator::fiveDimBin>(bin));
}

std::pair<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin>,
          boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> > antok::beamfileGenerator::fiveDimBin::divide(const int& dim,
                                                                                                                 const double& splitPoint) const
{

	std::vector<double> upper1 = _b;
	std::vector<double> lower2 = _a;
	upper1[dim] = splitPoint;
	lower2[dim] = splitPoint;
	boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> newBin1(new antok::beamfileGenerator::fiveDimBin(_a, upper1));
	boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> newBin2(new antok::beamfileGenerator::fiveDimBin(lower2, _b));
	newBin1->addNeighbor(newBin2);
	newBin2->addNeighbor(newBin1);
	for(std::set<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> >::iterator it = _neighbors.begin();
	    it != _neighbors.end();
	    ++it)
	{
		(*it)->removeNeighbor(shared_from_this());
		if(newBin1->areWeNeighbors(*it)) {
			newBin1->addNeighbor(*it);
			(*it)->addNeighbor(newBin1);
		}
		if(newBin2->areWeNeighbors(*it)) {
			newBin2->addNeighbor(*it);
			(*it)->addNeighbor(newBin2);
		}
	}
	newBin1->setOnLowerEdge(_onLowerEdge);
	newBin2->setOnLowerEdge(_onLowerEdge);
	newBin1->setOnUpperEdge(_onUpperEdge);
	newBin2->setOnUpperEdge(_onUpperEdge);
	if(_onLowerEdge[dim]) {
		newBin2->getOnLowerEdge()[dim] = false;
	}
	if(_onUpperEdge[dim]) {
		newBin1->getOnUpperEdge()[dim] = false;
	}

	for(unsigned int i = 0; i < 5; ++i) {

	}
	return std::pair<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin>,
	                 boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> >(newBin1, newBin2);

}

const std::vector<double>& antok::beamfileGenerator::fiveDimBin::getSigmas(unsigned int binContent,
                                                                           bool forceCalculation) const
{

	if(not _sigmas) {
		_sigmas = new std::vector<double>(5, 0.);
		forceCalculation = true;
	}
	if(forceCalculation) {
		for(unsigned int i = 0; i < 5; ++i) {
			*(_sigmas[i]) = (getUpperCorner()[i] - getLowerCorner()[i]) / (double)binContent;
		}
	}
	return *_sigmas;

}

std::ostream& antok::beamfileGenerator::fiveDimBin::print(std::ostream& out) const {
	out << "Five dimensional bin: " << std::endl;
	out << "    lower Corner ....... [" << _a[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _a[i];
	}
	out << "]" << std::endl;
	out << "    upper Corner ....... [" << _b[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _b[i];
	}
	out << "]" << std::endl;
	out << "    #neighbors ......... " << _neighbors.size() << std::endl;
	out << "    on lower edge ...... [" << _onLowerEdge[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _onLowerEdge[i];
	}
	out << "]" << std::endl;
	out << "    on upper edge ...... [" << _onUpperEdge[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _onUpperEdge[i];
	}
	out << "]" << std::endl;
	return out;
}

bool antok::beamfileGenerator::fiveDimBin::doubleEqual(const double& lhs, const double& rhs)
{

	if(std::fabs(lhs - rhs) < EPSILON) {
		return true;
	}
	return false;

}
