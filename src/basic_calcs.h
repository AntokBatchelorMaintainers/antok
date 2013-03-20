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

}

#endif

