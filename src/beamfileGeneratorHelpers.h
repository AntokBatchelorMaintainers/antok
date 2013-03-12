#ifndef ANTOK_BEAMFILEGENERATORHELPERS_H
#define ANTOK_BEAMFILEGENERATORHELPERS_H

#include<assert.h>
#include<iostream>
#include<list>
#include<vector>

class TTree;

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

		void getAdaptiveBins(std::list<std::pair<TTree*, antok::beamfileGenerator::fiveDimBin> >& bins,
		                     const antok::beamfileGenerator::fiveDimBin& bin,
		                     TTree* inputTree,
		                     unsigned int dim = 0,
		                     bool debug = false);

	}

}

#endif
