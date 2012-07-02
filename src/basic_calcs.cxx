
#include<basic_calcs.h>
#include<constants.h>

#include<TLorentzVector.h>
#include<TVector3.h>

TLorentzVector hlib::get_beam_energy(TVector3 p3_beam, TLorentzVector pX) {

	using hlib::PION_MASS;
	using hlib::PROTON_MASS;

	TVector3 p3_Tot(pX.Vect());
	double theta = p3_Tot.Angle(p3_beam);
	double E_tot = pX.E();
	double p3_Tot_Mag = p3_Tot.Mag();
	double a0 = (PION_MASS * PION_MASS) * p3_Tot_Mag * TMath::Cos(theta);
	double a1 = (PROTON_MASS * E_tot) - 0.5 * (pX.M2() + (PION_MASS * PION_MASS));
	double a2 = PROTON_MASS - E_tot + (p3_Tot_Mag * TMath::Cos(theta));
	double E_beam = (a1/a2) * (1 + (0.5 * (std::sqrt(1 + ((2*a2*a0)/(a1*a1)))-1)));
	double p_beam = std::sqrt((E_beam * E_beam) - (PION_MASS * PION_MASS));
	p3_beam.SetMag(p_beam);
	TLorentzVector pBeam;
	pBeam.SetXYZM(p3_beam.X(), p3_beam.Y(), p3_beam.Z(), PROTON_MASS);
	return pBeam;

};

void hlib::get_RPD_delta_phi_res(TLorentzVector pBeam,
								 TLorentzVector pProton,
								 TLorentzVector pX,
								 double& delta_phi, double& res)
{

	double rpd_Phi = pProton.Phi() / TMath::TwoPi();
	if(rpd_Phi < 0.) rpd_Phi += 1.;
	int hit_slat = (int)((rpd_Phi*24.) + 0.5);
	if(hit_slat >= 24) hit_slat -= 24;
	if(hit_slat % 2) res = 3.75;
	else res = 7.5;
	res = (res/180.)*TMath::Pi();
	res = std::sqrt(res*res + 0.067260*0.067260);

	TVector3 proton3 = pProton.Vect();
	TVector3 pTot3 = pX.Vect();
	TVector3 pBeam3 = pBeam.Vect();
	TVector3 a = pTot3 - (pTot3.Dot(pBeam3.Unit()))*pBeam3.Unit();
	TVector3 b = proton3 - (proton3.Dot(pBeam3.Unit()))*pBeam3.Unit();
	TVector3 n = a.Cross(b);
	if(n.Dot(TVector3(0.,0.,1.)) > 0) {
		n *= -1.;
	}
	double y = ((n.Cross(a)).Dot(b)) / (((a.Cross(b)).Cross(a)).Mag()*b.Mag());
	double x = (a.Unit()).Dot(b.Unit());
	delta_phi = TMath::ATan2(y, x)-TMath::Pi();
	if(delta_phi < (-1.*TMath::Pi())) {
		delta_phi += TMath::TwoPi();
	}

};
