#include<beamfile_generator_helpers.h>

#include<algorithm>
#include<assert.h>
#include<cmath>
#include<iostream>
#include<sstream>

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

void antok::beamfileGenerator::getAdaptiveBins(std::list<boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> >& bins,
                                               boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin,
                                               int dim,
                                               bool debug,
                                               unsigned int depth)
{

	long entries = bin->getEntries();
	const unsigned int indent = depth * antok::beamfileGenerator::INDENT;
	if(debug) {
		std::cout<<std::string(indent, ' ')<<"------------------------------------------------------"<<std::endl;
		std::cout<<std::string(indent, ' ')<<"called with " << entries << " entries at depth " << depth << "." << std::endl;
		std::cout<<std::string(indent, ' ')<<"input bin:"<<std::endl;
		bin->print(std::cout, depth);
	}
	if(dim == 0 and entries < (32 * antok::beamfileGenerator::MIN_ENTRIES)) {
		bins.push_back(bin);

	} else {
		std::pair<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin>,
		          boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> > newBins;
		{
			std::stringstream out;
			newBins = bin->divide(dim, debug ? &out : 0, depth);
			int entryDifference = std::abs((int)newBins.first->getEntries() - newBins.second->getEntries());
			if(debug) {
				std::cout<<out.str();
				std::cout<<std::string(indent, ' ')<<"bin1: "<<std::endl;
				newBins.first->print(std::cout, depth);
				std::cout<<std::string(indent, ' ')<<"bin2: "<<std::endl;
				newBins.second->print(std::cout, depth);
				std::cout<<std::string(indent, ' ')<<"inputTree1 has " << newBins.first->getEntries() << " entries." << std::endl;
				std::cout<<std::string(indent, ' ')<<"inputTree2 has " << newBins.second->getEntries() << " entries." << std::endl;
				std::cout<<std::string(indent, ' ')<<"Difference: " << entryDifference << std::endl;
			}
			if(entryDifference > 1) {
				std::cout<<"------------------------------------------------------"<<std::endl;
				std::cout<<"called with " << entries << " entries." << std::endl;
				std::cout<<"input bin:"<<std::endl;
				bin->print(std::cout);
				std::cout<<out.str();
				std::cout<<"inputTree1 has " << newBins.first->getEntries() << " entries." << std::endl;
				std::cout<<"inputTree2 has " << newBins.second->getEntries() << " entries." << std::endl;
				std::cout<<"Difference: " << entryDifference<< std::endl;
			}
		}
		if(dim == 4) {
			dim = 0;
		} else {
			++dim;
		}
		antok::beamfileGenerator::getAdaptiveBins(bins, newBins.first, dim, debug, depth + 1);
		antok::beamfileGenerator::getAdaptiveBins(bins, newBins.second, dim, debug, depth + 1);
	}
}
