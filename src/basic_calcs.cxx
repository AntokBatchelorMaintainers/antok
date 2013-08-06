#include<basic_calcs.h>

#include<iostream>

#include<TLorentzVector.h>
#include<TLorentzRotation.h>
#include<TVector3.h>

#include<constants.h>

TLorentzVector antok::getBeamEnergy(TVector3 p3_beam, const TLorentzVector& pX) {

	const double& PION_MASS = antok::Constants::chargedPionMass();
	const double& PROTON_MASS = antok::Constants::protonMass();

	TVector3 p3_Tot(pX.Vect());
	double theta = p3_Tot.Angle(p3_beam);
	double E_tot = pX.E();
	double p3_Tot_Mag = p3_Tot.Mag();
	double a0 = (PION_MASS * PION_MASS) * p3_Tot_Mag * TMath::Cos(theta);
	double a1 = (PROTON_MASS * E_tot) - 0.5 * (pX.M2() + (PION_MASS * PION_MASS));
	double a2 = PROTON_MASS - E_tot + (p3_Tot_Mag * TMath::Cos(theta));
	double E_beam = (a1/(2*a2)) * (1 + std::sqrt(1 + ((2*a2*a0)/(a1*a1))));
	double p_beam = std::sqrt((E_beam * E_beam) - (PION_MASS * PION_MASS));
	p3_beam.SetMag(p_beam);
	TLorentzVector pBeam;
	pBeam.SetXYZM(p3_beam.X(), p3_beam.Y(), p3_beam.Z(), PION_MASS);
	return pBeam;

};

void antok::getRPDDeltaPhiResProjection(const TLorentzVector& pBeam,
                                        const TLorentzVector& pProton,
                                        const TLorentzVector& pX,
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

void antok::getRPDDeltaPhiResRotation(const TLorentzVector& pBeam,
                                      const TLorentzVector& pProton,
                                      const TLorentzVector& pX,
                                      double& delta_phi, double& res)
{

	double rpd_Phi = pProton.Phi();
	rpd_Phi /= TMath::TwoPi();

	if(rpd_Phi < 0.) rpd_Phi += 1.;
	int hit_slat = (int)((rpd_Phi*24.) + 0.5);
	if(hit_slat >= 24) hit_slat -= 24;
	if(hit_slat % 2) res = 3.75;
	else res = 7.5;
	res = (res/180.)*TMath::Pi();

	TVector3 z(0.,0.,1.);
	double angle = pBeam.Vect().Angle(z);
	TVector3 ra = pBeam.Vect().Cross(z);
	TRotation r;
	r.Rotate(angle, ra);
	TLorentzRotation Tr(r);

	TLorentzVector proton_lab(pProton);
	TLorentzVector pX_lab(pX);
	proton_lab.Transform(Tr);
	pX_lab.Transform(Tr);

	delta_phi = std::fabs(pX_lab.Phi() - proton_lab.Phi()) - TMath::Pi();
	res = std::sqrt(res*res + 0.067260*0.067260);

};

void antok::getBoostToCenterOfMassSystem(const TLorentzVector& pBeam,
                                         double& centerOfMassEnergy,
                                         TVector3& boostVector)
{

	const double& PROTON_MASS = antok::Constants::protonMass();

	double beamEnergy = pBeam.E();
	centerOfMassEnergy = std::sqrt(pBeam.M2() + PROTON_MASS*PROTON_MASS + 2*beamEnergy*PROTON_MASS);
	TLorentzVector lorentzBoost = pBeam;
	lorentzBoost.SetE(beamEnergy + PROTON_MASS);
	boostVector = lorentzBoost.BoostVector();

}

void antok::getRPDExpectedHitsParameters(const TLorentzVector& pBeam,
                                         const TLorentzVector& pX,
                                         const TVector3& vertex,
                                         const double& xOffset,
                                         const double& yOffset,
                                         const double& xAngle,
                                         const double& yAngle,
                                         double& rpdPhi,
                                         double& rpdZRingA,
                                         double& rpdZRingB)
{

	static const double INNER_RING_DIAMETER = 12.0; // cm, from 2008 spectro paper
	static const double OUTER_RING_DIAMETER = 75.0; // cm, from 2008 spectro paper

	const TVector3& beamVector = pBeam.Vect();
	const TVector3& XVector = pX.Vect();

	const TVector3 proton = beamVector - XVector;

	// TRANSFORM `proton` AND `vertex` HERE
	//--------------------------------------
	if(xOffset != 0. or yOffset != 0. or xAngle != 0. or yAngle != 0.) {
		std::cerr<<"No transformations implemented in getRPDExpectedHitsParameters! Aborting..."<<std::endl;
		throw;
	}
	// REMOVE `#include<iostream>` AFTER REMOVING THIS BLOCK!
	//--------------------------------------

	rpdPhi = proton.Phi();
	const double tanTheta = TMath::Tan(proton.Theta());
	rpdZRingA = vertex.Z() + INNER_RING_DIAMETER / tanTheta;
	rpdZRingB = vertex.Z() + OUTER_RING_DIAMETER / tanTheta;

}
