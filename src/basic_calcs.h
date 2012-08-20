#ifndef BASIC_CALCS_H
#define BASIC_CALCS_H

class TLorentzVector;
class TVector3;

namespace hlib {

	TLorentzVector get_beam_energy(TVector3 p3_beam, const TLorentzVector& LV_X);

	void get_RPD_delta_phi_res(const TLorentzVector& pBeam,
							   const TLorentzVector& pProton,
							   const TLorentzVector& pX,
							   double& delta_phi, double& res);

}

#endif
