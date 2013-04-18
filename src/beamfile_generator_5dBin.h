#ifndef ANTOK_BEAMFILE_GENERATOR_5DBIN_H
#define ANTOK_BEAMFILE_GENERATOR_5DBIN_H

#include<iostream>
#include<vector>
#include<set>

#include<boost/shared_ptr.hpp>
#include<boost/enable_shared_from_this.hpp>


namespace antok {

	namespace beamfileGenerator {

		class fiveDimBin : public boost::enable_shared_from_this<fiveDimBin> {

		  public:

			fiveDimBin()
				: _a(5,0.),
				  _b(5,0.),
				  _neighbors(),
				  _onLowerEdge(5, false),
				  _onUpperEdge(5, false),
				  _sigmaCache(0) { _nExistingBins += 1; }

			fiveDimBin(double a0, double a1, double a2, double a3, double a4,
					   double b0, double b1, double b2, double b3, double b4);

			fiveDimBin(const std::vector<double>& a, const std::vector<double>& b);

			virtual ~fiveDimBin();

			void set(const std::vector<double>& a, const std::vector<double>& b);

			const std::vector<double>& getLowerCorner() const { return _a; }
			const std::vector<double>& getUpperCorner() const { return _b; }

			bool inBin(const std::vector<double>& x) const;

			std::pair<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin>,
			          boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> > divide(const int& dim, const double& splitPoint) const;

			bool areWeNeighbors(boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin) const;
			void addNeighbor(boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> bin) { _neighbors.insert(bin); }
			void removeNeighbor(boost::shared_ptr<const antok::beamfileGenerator::fiveDimBin> bin);
			const std::set<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> >& getNeighbors() const { return _neighbors; }

			void setOnLowerEdge(const std::vector<bool>& onLowerEdge) { _onLowerEdge = onLowerEdge; }
			void setOnUpperEdge(const std::vector<bool>& onUpperEdge) { _onUpperEdge = onUpperEdge; }
			const std::vector<bool>& getOnLowerEdge() const { return _onLowerEdge; }
			const std::vector<bool>& getOnUpperEdge() const { return _onUpperEdge; }
			std::vector<bool>& getOnLowerEdge() { return _onLowerEdge; }
			std::vector<bool>& getOnUpperEdge() { return _onUpperEdge; }

			double getVolume() const;

			const std::vector<double>& getSigmas(unsigned int binContent, bool forceCalculation = false) const;

			std::ostream& print(std::ostream& out, unsigned int indent = 0) const;

			static const long& getNExistingBins() { return _nExistingBins; }

			static void setDebug(bool debug = true) { _debug = debug; }
			static void setPrintNeighbors(bool printNeighbors = true) { _printNeighbors = printNeighbors; }

		  private:

			std::vector<double> _a; // lower corner
			std::vector<double> _b; // upper corner

			std::set<boost::shared_ptr<antok::beamfileGenerator::fiveDimBin> > _neighbors;
			std::vector<bool> _onLowerEdge;
			std::vector<bool> _onUpperEdge;

			mutable std::vector<double>* _sigmaCache;

			static const double EPSILON;
			static bool doubleEqual(const double& lhs, const double& rhs);

			static long _nExistingBins;

			static bool _debug;
			static bool _printNeighbors;

		};

	}

}

#endif
