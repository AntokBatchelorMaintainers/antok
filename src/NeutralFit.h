#ifndef NEUTRALFIT_H
#define NEUTRALFIT_H

#include <vector>

#include "TMatrixDSym.h"

#include "NeutralProblem.h"

class TH1;

namespace antok {
	class NeutralFit {
		public:
			NeutralFit(const TVector3 &vertexPosition_,
			           const TVector3 &cluster1Position_,
			           const TVector3 &cluster2Position_,
			           const TVector3 &cluster1PositionError_,
			           const TVector3 &cluster2PositionError_,
			           double cluster1Energy_,
			           double cluster2Energy_,
			           double cluster1EnergyError_,
			           double cluster2EnergyError_,
			           double mass_,
			           double window_ = 0.1,
			           int whichDeltaE_ = 0);

			~NeutralFit() { delete myFitter; delete myProblem; }

			bool isInWindow() const;
			bool doFit();

			double getChi2() { return myFitter->getChi2(); }
			double getCL()   { return myFitter->getCL();   }

			const TLorentzVector &getLV1() { return lv1; }
			const TLorentzVector &getLV2() { return lv2; }
			const TLorentzVector &getLVSum() { return lvSum; }

			TH1 *gethPull(size_t i) { return hPulls[i]; }

		private:
			static bool first;
			static std::vector<TH1 *> hPulls;

			const TVector3 &vertexPosition;
			const TVector3 &cluster1Position;
			const TVector3 &cluster2Position;
			const TVector3 &cluster1PositionError;
			const TVector3 &cluster2PositionError;

			double cluster1Energy;
			double cluster2Energy;
			double cluster1EnergyError;
			double cluster2EnergyError;

			TLorentzVector lv1;
			TLorentzVector lv2;
			TLorentzVector lvSum;

			double mass;
			double window;

			int whichDeltaE;

			TVectorD startingValues;
			TMatrixDSym covMat;

			NeutralProblem *myProblem;
			KinematicFit *myFitter;

			void initPulls();

			TMatrixDSym covMatForCluster(const TVector3 &clusterPosition,
			                             const TVector3 &clusterPositionError,
			                             const double clusterEnergy,
			                             const double clusterEnergyError);
	};
}

#endif
