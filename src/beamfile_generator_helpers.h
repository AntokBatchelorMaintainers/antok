#ifndef ANTOK_BEAMFILE_GENERATOR_HELPERS_H
#define ANTOK_BEAMFILE_GENERATOR_HELPERS_H

#include<list>
#include<vector>

#include<boost/shared_ptr.hpp>

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

	}

}

#endif
