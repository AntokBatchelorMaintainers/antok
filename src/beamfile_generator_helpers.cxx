#include<beamfile_generator_helpers.h>

#include<algorithm>
#include<assert.h>
#include<cmath>
#include<iostream>

#include<beamfile_generator_5dBin.h>

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
                                                                   boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> > >& bins,
                                               boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin,
                                               std::vector<antok::beamfileGenerator::fiveDimCoord*>* inputVector,
                                               int dim,
                                               bool debug,
                                               unsigned int depth)
{

	long entries = inputVector->size();
	if(debug) {
		std::cout<<"------------------------------------------------------"<<std::endl;
		std::cout<<"called with " << entries << " entries at depth " << depth << "." << std::endl;
		std::cout<<"input bin:"<<std::endl;
		bin->print(std::cout);
	}
	if(dim == 0 and entries < (32 * antok::beamfileGenerator::MIN_ENTRIES)) {
		bins.push_back(std::pair<std::vector<antok::beamfileGenerator::fiveDimCoord*>*,
		                         boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> >(inputVector, bin));

	} else {
		std::pair<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin>,
		          boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> > newBins;
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
				for(unsigned int i = 1; i < 5; ++i) {
					std::cout<<", "<<(*inputVector)[half]->_coords[i];
				}
				std::cout<<"]"<<std::endl;
				std::cout<<"dim: " << dim << ", middle: " << middle << std::endl;
			}
			newBins = bin->divide(dim, middle);
			if(debug) {
				std::cout<<"bin1: "<<std::endl;
				newBins.first->print(std::cout);
				std::cout<<"bin2: "<<std::endl;
				newBins.second->print(std::cout);
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
				bin->print(std::cout);
				std::cout<<"entries: " << entries << ", half: " << half << std::endl;
				std::cout<<"middle coordinates: ["<<(*inputVector)[half]->_coords[0];
				for(unsigned int i = 1; i < 5; ++i) {
					std::cout<<", "<<(*inputVector)[half]->_coords[i];
				}
				std::cout<<"]"<<std::endl;
				std::cout<<"dim: " << dim << ", middle: " << middle << std::endl;
				std::cout<<"inputTree1 has " << inputVector1->size() << " entries." << std::endl;
				std::cout<<"inputTree2 has " << inputVector2->size() << " entries." << std::endl;
				std::cout<<"Difference: " << std::abs((int)inputVector1->size() - (int)inputVector2->size()) << std::endl;
			}
			delete inputVector;
		}
		if(dim == 4) {
			dim = 0;
		} else {
			++dim;
		}
		antok::beamfileGenerator::getAdaptiveBins(bins, newBins.first, inputVector1, dim, debug, depth + 1);
		antok::beamfileGenerator::getAdaptiveBins(bins, newBins.second, inputVector2, dim, debug, depth + 1);
	}
}
