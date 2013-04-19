#include<beamfile_generator_5dBin.h>

#include<algorithm>
#include<assert.h>
#include<cmath>
#include<limits>
#include<sstream>

#include<beamfile_generator_helpers.h>


namespace {

	struct __compareCoords {
		bool operator ()(antok::beamfileGenerator::fiveDimCoord* lhs, antok::beamfileGenerator::fiveDimCoord* rhs)
		{
			return *lhs < *rhs;
		}
	};

	struct __compareCoordsDistance {

		__compareCoordsDistance(const antok::beamfileGenerator::fiveDimCoord& ref) : _ref(ref) { }
		bool operator()(antok::beamfileGenerator::fiveDimCoord* lhs, antok::beamfileGenerator::fiveDimCoord* rhs)
		{
			return lhs->distance(_ref) < rhs->distance(_ref);
		}
		const antok::beamfileGenerator::fiveDimCoord& _ref;
	};

}


const double antok::beamfileGenerator::fiveDimBin::EPSILON = 5. * std::numeric_limits<double>::epsilon();
long antok::beamfileGenerator::fiveDimBin::_nExistingBins = 0;
bool antok::beamfileGenerator::fiveDimBin::_debug = false;
bool antok::beamfileGenerator::fiveDimBin::_printNeighbors = false;
bool antok::beamfileGenerator::fiveDimBin::_differentSigmaCalculationForEdges = true;

antok::beamfileGenerator::fiveDimBin::fiveDimBin(double a0,
                                                 double a1,
                                                 double a2,
                                                 double a3,
                                                 double a4,
                                                 double b0,
                                                 double b1,
                                                 double b2,
                                                 double b3,
                                                 double b4,
                                                 std::vector<antok::beamfileGenerator::fiveDimCoord*>* entries)
	: _a(5, 0.),
	  _b(5, 0.),
	  _entries(entries),
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
	_nExistingBins += 1;
}

antok::beamfileGenerator::fiveDimBin::fiveDimBin(const std::vector<double>& a,
                                                 const std::vector<double>& b,
                                                 std::vector<antok::beamfileGenerator::fiveDimCoord*>* entries)
	: _a(5, 0.),
	  _b(5, 0.),
	  _entries(entries),
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
	_nExistingBins += 1;
}

antok::beamfileGenerator::fiveDimBin::~fiveDimBin()
{
	_neighbors.clear();
	_nExistingBins -= 1;
	delete _sigmaCache;
	delete _entries;
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
		if(bin->getLowerCorner()[i] > _b[i] or bin->getUpperCorner()[i] < _a[i]) {
			return false;
		}
	}
	return true;

}

void antok::beamfileGenerator::fiveDimBin::removeNeighbor(boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin) {
	_neighbors.erase(boost::const_pointer_cast<antok::beamfileGenerator::fiveDimBin>(bin));
}

std::pair<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin>,
          boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> > antok::beamfileGenerator::fiveDimBin::divide(const int& dim,
                                                                                                                 std::stringstream* outPtr,
                                                                                                                 unsigned int depth) const
{
	antok::beamfileGenerator::fiveDimCoord::_orderDim = dim;
	std::sort(_entries->begin(), _entries->end(), __compareCoords());
	const unsigned int half = (unsigned int)((_entries->size() / 2) + 0.5);
	const double splitPoint = (*_entries)[half]->_coords[dim];
	if(outPtr) {
		std::stringstream& out = *outPtr;
		const unsigned int indent = depth * antok::beamfileGenerator::INDENT;
		out << std::string(indent, ' ') << "entries: " << _entries->size() << ", half: " << half << std::endl;
		out << std::string(indent, ' ') << "middle coordinates: [" << (*_entries)[half]->_coords[0];
		for(unsigned int i = 1; i < 5; ++i) {
			out << ", " << (*_entries)[half]->_coords[i];
		}
		out << "]" << std::endl;
		out << std::string(indent, ' ') << "dim: " << dim << ", middle: " << splitPoint << std::endl;
	}

	std::vector<double> upper1 = _b;
	std::vector<double> lower2 = _a;
	upper1[dim] = splitPoint;
	lower2[dim] = splitPoint;
	std::vector<antok::beamfileGenerator::fiveDimCoord*>* entries1 = new std::vector<antok::beamfileGenerator::fiveDimCoord*>();
	std::vector<antok::beamfileGenerator::fiveDimCoord*>* entries2 = new std::vector<antok::beamfileGenerator::fiveDimCoord*>();
	for(unsigned int i = 0; i < _entries->size(); ++i) {
		const double& x = (*_entries)[i]->_coords[dim];
		if(x < splitPoint) {
			entries1->push_back((*_entries)[i]);
		} else {
			entries2->push_back((*_entries)[i]);
		}
	}
	boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> newBin1(new antok::beamfileGenerator::fiveDimBin(_a, upper1, entries1));
	boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> newBin2(new antok::beamfileGenerator::fiveDimBin(lower2, _b, entries2));
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
	return std::pair<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin>,
	                 boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> >(newBin1, newBin2);

}

double antok::beamfileGenerator::fiveDimBin::getVolume() const
{
	double binVolume = 1.;
	for(unsigned int i = 0; i < 5; ++i) {
		binVolume *= getUpperCorner()[i] - getLowerCorner()[i];
	}
	return binVolume;
}

unsigned int antok::beamfileGenerator::fiveDimBin::getEdgeity() const {

	unsigned int edgeity = 0;
	for(unsigned int i = 0; i < 5; ++i) {
		edgeity += _onLowerEdge[i] + _onUpperEdge[i];
	}
	return edgeity;

}

const std::vector<std::vector<double> >& antok::beamfileGenerator::fiveDimBin::getSigmas(int& method, bool forceCalculation) const
{

	if(not _sigmaCache) {
		_sigmaCache = new std::vector<std::vector<double> >(_entries->size(), std::vector<double>(5, 0.));
		forceCalculation = true;
	}
	if(forceCalculation) {

		// get the scaling factors
		std::vector<double> scalingFactors(4, 0.);
		double scalingFactorProduct = 1.;
		double scalingNorm = (getUpperCorner()[0] - getLowerCorner()[0]);
		if(_debug) {
			std::cout<<std::endl<<"calculating sigmas for bin:"<<std::endl;
			print(std::cout);
		}
		for(unsigned int i = 0; i < 5; ++i) {
			if(_debug) {
				std::cout<<"edge "<<i<<" length: "<< getUpperCorner()[i] - getLowerCorner()[i]<<std::endl;
			}
			if(i > 0) {
				scalingFactors[i-1] = (getUpperCorner()[i] - getLowerCorner()[i]) / scalingNorm;
				if(_debug) {
					std::cout<<"scaling factor "<<i<<": "<< scalingFactors[i-1]<<std::endl;
				}
				scalingFactorProduct *= scalingFactors[i-1];
			}
		}

		// we are at the edge, have to be careful
		if(_differentSigmaCalculationForEdges and getEdgeity() > 0) {
			if(_debug) {
				std::cout<<"on an edge bin"<<std::endl;
			}
			std::vector<antok::beamfileGenerator::fiveDimCoord*> allEvents(*_entries);
			std::vector<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> > relevantNeighbors;
			for(std::set<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> >::iterator it = _neighbors.begin();
			    it != _neighbors.end();
			    ++it)
			{
				const std::vector<antok::beamfileGenerator::fiveDimCoord*>* eventsFromBin = (*it)->getEvents();
				for(unsigned int i = 0; i < eventsFromBin->size(); ++i) {
					allEvents.push_back((*eventsFromBin)[i]);
				}
				if((*it)->getEdgeity() == 0) {
					relevantNeighbors.push_back(*it);
				}
			}
			// ok, we got all the info from the neighboring bins
			if(_debug) {
			}

			if(relevantNeighbors.size() > 0) {
			// we have a bin which is not on the edge, take the sigmas from there
				method = 1;
				if(_debug) {
					std::cout<<"found "<<relevantNeighbors.size()<<" neighbors which are"
					         <<" not on the edge (out of "<<_neighbors.size()<<")"<<std::endl;
				}
				std::vector<double> averagedSigmas(5, 0.);
				for(unsigned int i = 0; i < relevantNeighbors.size(); ++i) {
					bool debugCache = false;
					if(_debug) {
						debugCache = true;
						antok::beamfileGenerator::fiveDimBin::setDebug(false);
					}
					int dump;
					const std::vector<double>& sigmas = relevantNeighbors[i]->getSigmas(dump)[0];
					if(debugCache) {
						antok::beamfileGenerator::fiveDimBin::setDebug(true);
					}
					if(_debug) {
						std::cout<<"neighbor "<<i<<" sigmas: ["<<sigmas[0];
					}
					for(unsigned int j = 0; j < 5; ++j) {
						averagedSigmas[j] += sigmas[j];
						if(_debug) {
							std::cout<<", "<<sigmas[j];
						}
					}
					if(_debug) {
						std::cout<<"]"<<std::endl;
					}
				}
				if(_debug) {
					std::cout<<"final sigmas for this bin: ["<<averagedSigmas[0] / relevantNeighbors.size();
				}
				for(unsigned int j = 0; j < 5; ++j) {
					averagedSigmas[j] /= relevantNeighbors.size();
					if(_debug) {
						std::cout<<", "<<averagedSigmas[j];
					}
				}
				if(_debug) {
					std::cout<<"]"<<std::endl;
				}
				for(unsigned int i = 0; i < _entries->size(); ++i) {
								(*_sigmaCache)[i] = averagedSigmas;;
				}
			} else {
			// damn, no luck, we have to calculate the sigmas the hard way
				method = 2;
				static const unsigned int ENTRIES_IN_SPHERE = 15;
				for(unsigned int i = 0; i < _entries->size(); ++i) {
					// sort all events in ascending distance to the current event
					std::sort(allEvents.begin(), allEvents.end(), __compareCoordsDistance(*((*_entries)[i])));

					// get the center of gravity of the closest points
					antok::beamfileGenerator::fiveDimCoord centerOfGravity(*((*_entries)[i]));
					for(unsigned int j = 1; j <= ENTRIES_IN_SPHERE; ++j) {
						centerOfGravity += *(allEvents[j]);
					}
					centerOfGravity /= ENTRIES_IN_SPHERE + 1;

					// get radius
					double radius = 0.;
					for(unsigned int j = 0; j <= ENTRIES_IN_SPHERE; ++j) {
						if(centerOfGravity.distance(*(allEvents[j])) > radius) {
							radius = centerOfGravity.distance(*(allEvents[j]));
						}
					}
					double volume = 5.263789013914324 * pow(radius, 5);

					// now sort everything in ascending distance to the center of gravity
					std::sort(allEvents.begin(), allEvents.end(), __compareCoordsDistance(centerOfGravity));

					// and see how many entries we have in this new sphere
					unsigned int entriesInShiftedSphere = 0;
					while(allEvents[entriesInShiftedSphere]->distance(centerOfGravity) < radius) {
						++entriesInShiftedSphere;
					}

					// and now get the first sigma
					(*_sigmaCache)[i][0] = std::pow(volume / (((double)entriesInShiftedSphere) * scalingFactorProduct), 0.2);

					for(unsigned int j = 0; j < 4; ++j) {
						(*_sigmaCache)[i][j+1] = scalingFactors[j] * (*_sigmaCache)[i][0];
					}
					if(_debug) {
						std::cout<<"sigmas = ["<<(*_sigmaCache)[i][0];
						for(unsigned int j = 1; j < 5; ++j) {
							std::cout<<", "<<(*_sigmaCache)[i][j];
						}
						std::cout<<"]"<<std::endl;
						std::cout<<"calculated with "<<entriesInShiftedSphere<<" events."<<std::endl;
						std::cout<<"ref"<<std::endl;
						(*_entries)[i]->print(std::cout);
						std::cout<<"---------"<<std::endl;
						for(unsigned int j = 0; j < entriesInShiftedSphere; ++j) {
							allEvents[j]->print(std::cout);
							std::cout<<"dist: "<<(*_entries)[i]->distance(*(allEvents[j]))<<std::endl;
						}
						std::cout<<"---------"<<std::endl<<std::endl;
					}
				}
			}

		// we are not at the edge, life is easy
		} else {
			method = 0;
			double binVolume = getVolume();
			if(_debug) {
				std::cout<<"on an inner bin"<<std::endl;
				std::cout<<std::endl;
				this->print(std::cout);
				std::cout<<"binContent: "<<getEntries()<<std::endl;
			}
			std::vector<double> sigmas(5, 0.);
			sigmas[0] = std::pow(binVolume / (((double)getEntries()) * scalingFactorProduct), 0.2);
			if(_debug) {
				std::cout<<"scaling factor product: "<< scalingFactorProduct<<std::endl;
				std::cout<<"binVolume: "<<binVolume<<std::endl;
				std::cout<<"sigma "<<0<<": "<<sigmas[0]<<std::endl;
			}
			for(unsigned int i = 0; i < 4; ++i) {
				sigmas[i+1] = scalingFactors[i] * sigmas[0];
				if(_debug) {
					std::cout<<"sigma "<<i+1<<": "<<sigmas[i+1]<<std::endl;
				}
			}
			if(_debug) {
				std::cout<<std::endl;
			}
			for(unsigned int i = 0; i < _entries->size(); ++i) {
				(*_sigmaCache)[i] = sigmas;
			}
		}
		if(_debug) {
			std::cout<<std::endl;
		}
	}
	return *_sigmaCache;

}

std::ostream& antok::beamfileGenerator::fiveDimBin::print(std::ostream& out, unsigned int indent) const
{
	indent *= antok::beamfileGenerator::INDENT;
	out << std::string(indent, ' ') << "Five dimensional bin: " << std::endl;
	out << std::string(indent, ' ') << "    lower Corner ....... [" << _a[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _a[i];
	}
	out << "]" << std::endl;
	out << std::string(indent, ' ') << "    upper Corner ....... [" << _b[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _b[i];
	}
	out << "]" << std::endl;
	out << std::string(indent, ' ') << "    on lower edge ...... [" << _onLowerEdge[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _onLowerEdge[i];
	}
	out << "]" << std::endl;
	out << std::string(indent, ' ') << "    on upper edge ...... [" << _onUpperEdge[0];
	for(unsigned int i = 1; i < 5; ++i) {
		out << ", " << _onUpperEdge[i];
	}
	out << "]" << std::endl;
	out << std::string(indent, ' ') << "    #neighbors ......... " << _neighbors.size() << std::endl;
	if(_printNeighbors) {
		int neighborIndent = (indent / antok::beamfileGenerator::INDENT) + 1.5;
		_printNeighbors = false;
		for(std::set<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> >::const_iterator it = _neighbors.begin();
		    it != _neighbors.end();
		    ++it)
		{
			(*it)->print(out, neighborIndent);
		}
		_printNeighbors = true;
	}
	return out;
}

bool antok::beamfileGenerator::fiveDimBin::doubleEqual(const double& lhs, const double& rhs)
{

	if(std::fabs(lhs - rhs) < EPSILON) {
		return true;
	}
	return false;

}
