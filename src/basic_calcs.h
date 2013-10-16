#ifndef ANTOK_BASIC_CALCS_H
#define ANTOK_BASIC_CALCS_H

class TLorentzVector;
class TVector3;

namespace antok {

	TLorentzVector getBeamEnergy(TVector3 p3_beam, const TLorentzVector& LV_X);

	void getBoostToCenterOfMassSystem(const TLorentzVector& pBeam,
	                                  double& centerOfMassEnergy,
	                                  TVector3& boostVector);

	namespace utils {

		double getPositiveSolutionOfQuadraticEquation(const double& a,
		                                              const double& b,
		                                              const double& c);

	}

}

#endif
