#include<beamfileGeneratorHelpers.h>

#include<algorithm>
#include<assert.h>
#include<cmath>

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
	  _b(5, 0.)
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
}

antok::beamfileGenerator::fiveDimBin::fiveDimBin(const std::vector<double>& a, const std::vector<double>& b)
	: _a(5, 0.),
	  _b(5, 0.)
{
	assert(a.size() == 5);
	assert(b.size() == 5);
	for(unsigned int i = 0; i < 5; ++i) {
		_a[i] = a[i] < b[i] ? a[i] : b[i];
		_b[i] = a[i] < b[i] ? b[i] : a[i];
	}
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
	return out;
}

int antok::beamfileGenerator::fiveDimCoord::_orderDim = 0;

antok::beamfileGenerator::fiveDimCoord::fiveDimCoord()
{
	_coords.resize(5, 0.);
}

antok::beamfileGenerator::fiveDimCoord::fiveDimCoord(double x0, double x1, double x2, double x3, double x4)
{
	_coords.resize(5, 0.);
	_coords[0] = x0;
	_coords[1] = x1;
	_coords[2] = x2;
	_coords[3] = x3;
	_coords[4] = x4;
}

bool antok::beamfileGenerator::fiveDimCoord::operator<(const antok::beamfileGenerator::fiveDimCoord& rhs) const
{
	return _coords[_orderDim] < rhs._coords[_orderDim];
}

namespace {

	struct __compareCoords {
		bool operator ()(antok::beamfileGenerator::fiveDimCoord *lhs, antok::beamfileGenerator::fiveDimCoord *rhs)
		{
			return *lhs < *rhs;
		}
	};

}

void antok::beamfileGenerator::getAdaptiveBins(std::list<std::pair<std::vector<antok::beamfileGenerator::fiveDimCoord*>*,
                                                         antok::beamfileGenerator::fiveDimBin> >& bins,
                                               const antok::beamfileGenerator::fiveDimBin& bin,
                                               std::vector<antok::beamfileGenerator::fiveDimCoord*>* inputVector,
                                               int dim,
                                               bool debug)
{
	long entries = inputVector->size();
	if(debug) {
		std::cout<<"------------------------------------------------------"<<std::endl;
		std::cout<<"called with " << entries << " entries." << std::endl;
		std::cout<<"input bin:"<<std::endl;
		bin.print(std::cout);
	}
	if(entries < (2 * antok::beamfileGenerator::MIN_ENTRIES)) {
		bins.push_back(std::pair<std::vector<antok::beamfileGenerator::fiveDimCoord*>*,
				                 antok::beamfileGenerator::fiveDimBin>(inputVector, bin));
	} else {
		antok::beamfileGenerator::fiveDimBin newBin1;
		antok::beamfileGenerator::fiveDimBin newBin2;
		std::vector<antok::beamfileGenerator::fiveDimCoord*>* inputVector1 = 0;
		std::vector<antok::beamfileGenerator::fiveDimCoord*>* inputVector2 = 0;
		{
			antok::beamfileGenerator::fiveDimCoord::_orderDim = dim;
			std::sort(inputVector->begin(), inputVector->end(), __compareCoords());
			unsigned int half = (unsigned int)((entries / 2) + 0.5);
			const double& middle = (*inputVector)[half]->_coords[dim];
			if(debug) {
				std::cout<<"entries: " << entries << ", half: " << half << std::endl;
				std::cout<<"middle coordinates: ["<<(*inputVector)[half]->_coords[0];
				for(unsigned int i = 1; i < 4; ++i) {
					std::cout<<", "<<(*inputVector)[half]->_coords[i];
				}
				std::cout<<"]"<<std::endl;
				std::cout<<"dim: " << dim << ", middle: " << middle << std::endl;
			}
			std::vector<double> upper1 = bin._b;
			std::vector<double> lower2 = bin._a;
			upper1[dim] = middle;
			lower2[dim] = middle;
			newBin1.set(bin._a, upper1);
			newBin2.set(lower2, bin._b);
			if(debug) {
				std::cout<<"bin1: "<<std::endl;
				newBin1.print(std::cout);
				std::cout<<"bin2: "<<std::endl;
				newBin2.print(std::cout);
			}
			inputVector1 = new std::vector<antok::beamfileGenerator::fiveDimCoord*>();
			inputVector2 = new std::vector<antok::beamfileGenerator::fiveDimCoord*>();
			for(long i = 0; i < entries; ++i) {
				const double& x = (*inputVector)[i]->_coords[dim];
				if(x < middle) {
					inputVector1->push_back((*inputVector)[i]);
				} else {
					inputVector2->push_back((*inputVector)[i]);
				}
			}
			if(debug) {
				std::cout<<"inputTree1 has " << inputVector1->size() << " entries." << std::endl;
				std::cout<<"inputTree2 has " << inputVector2->size() << " entries." << std::endl;
				std::cout<<"Difference: " << std::abs((int)inputVector1->size() - (int)inputVector2->size()) << std::endl;
			}
			if(std::abs((int)inputVector1->size() - (int)inputVector2->size()) > 1) {
				std::cout<<"------------------------------------------------------"<<std::endl;
				std::cout<<"called with " << entries << " entries." << std::endl;
				std::cout<<"input bin:"<<std::endl;
				bin.print(std::cout);
				std::cout<<"entries: " << entries << ", half: " << half << std::endl;
				std::cout<<"middle coordinates: ["<<(*inputVector)[half]->_coords[0];
				for(unsigned int i = 1; i < 4; ++i) {
					std::cout<<", "<<(*inputVector)[half]->_coords[i];
				}
				std::cout<<"]"<<std::endl;
				std::cout<<"dim: " << dim << ", middle: " << middle << std::endl;
				std::cout<<"inputTree1 has " << inputVector1->size() << " entries." << std::endl;
				std::cout<<"inputTree2 has " << inputVector2->size() << " entries." << std::endl;
				std::cout<<"Difference: " << std::abs((int)inputVector1->size() - (int)inputVector2->size()) << std::endl;
//				assert(false);
			}
			delete inputVector;
		}
		if(dim == 4) {
			dim = 0;
		} else {
			++dim;
		}
		antok::beamfileGenerator::getAdaptiveBins(bins, newBin1, inputVector1, dim, debug);
		antok::beamfileGenerator::getAdaptiveBins(bins, newBin2, inputVector2, dim, debug);
	}
}
