#ifndef ANTOK_BEAMFILE_GENERATOR_HELPERS_H
#define ANTOK_BEAMFILE_GENERATOR_HELPERS_H

#include<list>
#include<vector>

#include<boost/shared_ptr.hpp>

#include<beamfile_generator_5dBin.h>
#include<beamfile_generator_5dCoord.h>

namespace antok {

	namespace beamfileGenerator {

		class fiveDimBin;

		const double MIN_ENTRIES = 10.;
		const unsigned int INDENT = 4;

		void getAdaptiveBins(std::list<boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> >& bins,
		                     boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin,
		                     int dim = 0,
		                     bool debug = false,
		                     unsigned int depth = 0);

		struct eventBookkeeper {

			eventBookkeeper()
			  : binVolume(-1.),
				binContent(-1),
				nNeighbors(-1),
				edgeity(-1),
				sigmaCalculationMethod(-1),
				sigmas(),
				sigmaIndex(0),
				coords(0) { }
	
			double binVolume;
			int binContent;
			int nNeighbors;
			int edgeity;
			int sigmaCalculationMethod;
			boost::shared_ptr<const std::vector<std::vector<double> > > sigmas;
			unsigned int sigmaIndex;
			const antok::beamfileGenerator::fiveDimCoord* coords;
	
		};

		struct compareEventBookkeepers {
			bool operator ()(const eventBookkeeper& lhs, const eventBookkeeper& rhs)
			{
				return lhs.coords->_eventNumber < rhs.coords->_eventNumber;
			}
		};

	}

}

#endif
