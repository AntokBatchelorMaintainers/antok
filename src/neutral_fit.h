#ifndef NEUTRALFIT_H
#define NEUTRALFIT_H

#include <vector>

#include "TMatrixDSym.h"

#include "neutral_problem.h"

class TH1;

namespace antok {

	class NeutralFit {

	public:

		NeutralFit(const TVector3& vertexPosition,
		           const TVector3& cluster1Position,
		           const TVector3& cluster2Position,
		           const TVector3& cluster1PositionError,
		           const TVector3& cluster2PositionError,
		           const double    cluster1Energy,
		           const double    cluster2Energy,
		           const double    cluster1EnergyError,
		           const double    cluster2EnergyError,
		           const double    mass,
		           const double    window = 0.1,
		           const int       whichDeltaE = 0);
		~NeutralFit()
		{
			delete _myFitter;
			delete _myProblem;
		}

		bool isInWindow() const;
		bool doFit();

		double getChi2() const { return _myFitter->getChi2(); }
		double getCL()   const { return _myFitter->getCL();   }

		const TLorentzVector& getLV1()   const { return _lv1;   }
		const TLorentzVector& getLV2()   const { return _lv2;   }
		const TLorentzVector& getLVSum() const { return _lvSum; }

		TH1* gethPull(const size_t i) { return _hPulls[i]; }

		std::vector<double> getPulls() const { return _pulls; }

	private:

		void initPulls();
		void fillPulls(const TVectorD& enhanced);
		TMatrixDSym covMatForCluster(const TVector3& clusterPosition,
		                             const TVector3& clusterPositionError,
		                             const double    clusterEnergy,
		                             const double    clusterEnergyError);

		static bool              _first;
		static std::vector<TH1*> _hPulls;
		std::vector<double>      _pulls;

		const TVector3& _vertexPosition;
		const TVector3& _cluster1Position;
		const TVector3& _cluster2Position;
		const TVector3& _cluster1PositionError;
		const TVector3& _cluster2PositionError;

		const double _cluster1Energy;
		const double _cluster2Energy;
		const double _cluster1EnergyError;
		const double _cluster2EnergyError;

		TLorentzVector _lv1;
		TLorentzVector _lv2;
		TLorentzVector _lvSum;

		const double _mass;
		const double _window;
		const int    _whichDeltaE;

		TVectorD    _startingValues;
		TMatrixDSym _covMat;

		antok::NeutralProblem* _myProblem;
		antok::KinematicFit*   _myFitter;

	};

}  // antok namespace

#endif  // NEUTRALFIT_H
