#ifndef ANTOK_BEAMFILE_GENERATOR_5DCOORD
#define ANTOK_BEAMFILE_GENERATOR_5DCOORD

#include<vector>
#include<iostream>

namespace antok {

	namespace beamfileGenerator {

		class fiveDimCoord {

		public:

			fiveDimCoord();

			fiveDimCoord(std::vector<double> coords, long eventNumber = -1)
				: _coords(coords),
				  _eventNumber(eventNumber) { }

			fiveDimCoord(double x0, double x1, double x2, double x3, double x4, long eventNumber = -1);

			bool operator<(const fiveDimCoord& rhs) const;
			antok::beamfileGenerator::fiveDimCoord& operator+=(const antok::beamfileGenerator::fiveDimCoord& rhs);
			antok::beamfileGenerator::fiveDimCoord& operator-=(const antok::beamfileGenerator::fiveDimCoord& rhs);
			antok::beamfileGenerator::fiveDimCoord& operator*=(const double& factor);
			antok::beamfileGenerator::fiveDimCoord& operator/=(const double& factor);
			const antok::beamfileGenerator::fiveDimCoord operator+(const antok::beamfileGenerator::fiveDimCoord& rhs);
			const antok::beamfileGenerator::fiveDimCoord operator-(const antok::beamfileGenerator::fiveDimCoord& rhs);

			double distance(const fiveDimCoord& point) const;
			double distance2(const fiveDimCoord& point) const;

			std::ostream& print(std::ostream& out) const;

			static int _orderDim;
			std::vector<double> _coords;
			long _eventNumber;

		};

	}

}

#endif
