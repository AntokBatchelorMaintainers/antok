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

		struct fiveDimCoord {

			fiveDimCoord();

			fiveDimCoord(std::vector<double> coords)
				: _coords(coords) { }

			fiveDimCoord(double x0, double x1, double x2, double x3, double x4);

			bool operator<(const fiveDimCoord& rhs) const;

			static int _orderDim;
			std::vector<double> _coords;

		};

		void getAdaptiveBins(std::list<std::pair<std::vector<antok::beamfileGenerator::fiveDimCoord*>*,
		                                         boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> > >& bins,
		                     boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin,
		                     std::vector<antok::beamfileGenerator::fiveDimCoord*>* inputVector,
		                     int dim = 0,
		                     bool debug = false,
		                     unsigned int depth = 0);

	}

}

#endif
