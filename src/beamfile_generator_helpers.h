#ifndef ANTOK_BEAMFILEGENERATORHELPERS_H
#define ANTOK_BEAMFILEGENERATORHELPERS_H

#include<iostream>
#include<list>
#include<vector>

namespace antok {

	namespace beamfileGenerator {

		const double MIN_ENTRIES = 10.;

		struct fiveDimBin {

			fiveDimBin() : _a(5,0), _b(5,0) { }

			fiveDimBin(double a0, double a1, double a2, double a3, double a4,
					   double b0, double b1, double b2, double b3, double b4);

			fiveDimBin(const std::vector<double>& a, const std::vector<double>& b);

			void set(const std::vector<double>& a, const std::vector<double>& b);

			bool inBin(const std::vector<double>& x) const;

			std::vector<double> _a;
			std::vector<double> _b;

			std::ostream& print(std::ostream& out) const;

		};

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
		                               antok::beamfileGenerator::fiveDimBin> >& bins,
		                     const antok::beamfileGenerator::fiveDimBin& bin,
		                     std::vector<antok::beamfileGenerator::fiveDimCoord*>* inputVector,
		                     int dim = 0,
		                     bool debug = false);

	}

}

#endif
