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
			antok::beamfileGenerator::fiveDimCoord& operator+=(const antok::beamfileGenerator::fiveDimCoord& rhs);
			antok::beamfileGenerator::fiveDimCoord& operator-=(const antok::beamfileGenerator::fiveDimCoord& rhs);
			antok::beamfileGenerator::fiveDimCoord& operator*=(const double& factor);
			antok::beamfileGenerator::fiveDimCoord& operator/=(const double& factor);
			const antok::beamfileGenerator::fiveDimCoord operator+(const antok::beamfileGenerator::fiveDimCoord& rhs);
			const antok::beamfileGenerator::fiveDimCoord operator-(const antok::beamfileGenerator::fiveDimCoord& rhs);

			double distance(const fiveDimCoord& point) const;

			std::ostream& print(std::ostream& out) const;

			static int _orderDim;
			std::vector<double> _coords;

		};

		void getAdaptiveBins(std::list<boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> >& bins,
		                     boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin,
		                     int dim = 0,
		                     bool debug = false,
		                     unsigned int depth = 0);

	}

}

#endif
