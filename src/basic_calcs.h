#ifndef ANTOK_BASIC_CALCS_H
#define ANTOK_BASIC_CALCS_H


class TLorentzVector;
class TVector3;

namespace antok {

	TLorentzVector getBeamEnergy(TVector3 p3_beam, const TLorentzVector& LV_X, const double beam_mass, const double target_mass);
	/***
	 * Uses pion mass for the beam and proton mass for the target.
	 */
	TLorentzVector getBeamEnergy(TVector3 p3_beam, const TLorentzVector& LV_X);

	void getBoostToCenterOfMassSystem(const TLorentzVector& pBeam,
	                                  double& centerOfMassEnergy,
	                                  TVector3& boostVector);

	template<typename T> double sqr(T t) { return t*t; }

	namespace utils {

		double getPositiveSolutionOfQuadraticEquation(const double& a,
		                                              const double& b,
		                                              const double& c);

	}

}

#endif
