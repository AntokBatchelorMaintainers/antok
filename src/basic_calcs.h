#ifndef ANTOK_BASIC_CALCS_H
#define ANTOK_BASIC_CALCS_H

class TLorentzVector;
class TVector3;

namespace antok {

	TLorentzVector getBeamEnergy(TVector3 p3_beam, const TLorentzVector& LV_X);

	void getRPDDeltaPhiResProjection(const TLorentzVector& pBeam,
	                                 const TLorentzVector& pProton,
	                                 const TLorentzVector& pX,
	                                 double& delta_phi, double& res);

	void getRPDDeltaPhiResRotation(const TLorentzVector& pBeam,
	                               const TLorentzVector& pProton,
	                               const TLorentzVector& pX,
	                               double& delta_phi, double& res);

	void getBoostToCenterOfMassSystem(const TLorentzVector& pBeam,
	                                  double& centerOfMassEnergy,
	                                  TVector3& boostVector);

	void getRPDExpectedHitsParameters(const TLorentzVector& pBeam,
	                                  const TLorentzVector& pX,
	                                  const TVector3& vertex,
	                                  const double& xOffset,
	                                  const double& yOffset,
	                                  const double& xAngle,
	                                  const double& yAngle,
	                                  double& rpdPhi,
	                                  double& rpdZRingA,
	                                  double& rpdZRingB);

}

#endif

