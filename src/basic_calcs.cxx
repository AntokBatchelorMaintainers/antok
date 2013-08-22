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
                                      double& delta_phi, double& res,
                                      double& phiProton, double& phiX)
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

	phiProton = proton_lab.Phi();
	phiX = pX_lab.Phi();

	delta_phi = std::fabs(phiX - phiProton) - TMath::Pi();
	res = std::sqrt(res*res + 0.067260*0.067260);

};

void antok::getRPDDeltaPhiResRotation(const TLorentzVector& pBeam,
                                      const TLorentzVector& pProton,
                                      const TLorentzVector& pX,
                                      double& delta_phi, double& res)
{
	double phiProton, phiX;
	getRPDDeltaPhiResRotation(pBeam, pProton, pX, delta_phi, res, phiProton, phiX);
}

void antok::getRPDDeltaPhiResPrediction(const TLorentzVector& pBeam,
                                        const TLorentzVector& pProton,
                                        const TLorentzVector& pX,
                                        double& delta_phi, double& res,
                                        double& phiProton, double& phiX)
{

	double rpd_Phi = pProton.Phi();
	rpd_Phi /= TMath::TwoPi();

	if(rpd_Phi < 0.) rpd_Phi += 1.;
	int hit_slat = (int)((rpd_Phi*24.) + 0.5);
	if(hit_slat >= 24) hit_slat -= 24;
	if(hit_slat % 2) res = 3.75;
	else res = 7.5;
	res = (res/180.)*TMath::Pi();

	const TVector3& beamVector = pBeam.Vect();
	const TVector3& XVector = pX.Vect();

	TVector3 proton = beamVector - XVector;

	phiX = proton.Phi();
	phiProton = pProton.Phi();

	delta_phi = phiX - phiProton;

	if(delta_phi < -TMath::Pi()) {
		delta_phi += TMath::TwoPi();
	}
	if(delta_phi > TMath::Pi()) {
		delta_phi -= TMath::TwoPi();
	}

	res = std::sqrt(res*res + 0.067260*0.067260);

};

void antok::getRPDDeltaPhiResPrediction(const TLorentzVector& pBeam,
                                      const TLorentzVector& pProton,
                                      const TLorentzVector& pX,
                                      double& delta_phi, double& res)
{
	double phiProton, phiX;
	getRPDDeltaPhiResPrediction(pBeam, pProton, pX, delta_phi, res, phiProton, phiX);
}

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

namespace {

	double __getPositiveSolutionOfQuadraticEquation(const double& a, const double& b, const double& c) {
		double result = (-b + std::sqrt(b*b - 4*a*c)) / (2*a);
		if(result < 0) {
			result = (-b - std::sqrt(b*b - 4*a*c)) / (2*a);
		}
		if(result < 0) {
			throw 42;
		}
		return result;
	}

}

void antok::getRPDExpectedHitsParameters(const TLorentzVector& pBeam,
                                         const TLorentzVector& pX,
                                         const TVector3& vertex,
                                         const double& xOffset,
                                         const double& yOffset,
                                         const double& xAngle,
                                         const double& yAngle,
                                         double& rpdPhiRingA,
                                         double& rpdPhiRingB,
                                         double& rpdZRingA,
                                         double& rpdZRingB)
{

	static const double INNER_RING_DIAMETER = 12.0; // cm, from 2008 spectro paper
	static const double OUTER_RING_DIAMETER = 75.0; // cm, from 2008 spectro paper
	static const double RPD_CENTER_Z_POSITION = -36.2; // cm, from COMGEANT's geom_hadron_2008.ffr

	const TVector3& beamVector = pBeam.Vect();
	const TVector3& XVector = pX.Vect();

	TVector3 proton = beamVector - XVector;

	TRotation rotation;
	rotation.SetToIdentity();
	rotation.RotateX(yAngle);
	rotation.RotateY(xAngle);
	rotation.Invert();

	TVector3 transformedVertex(vertex);
	const TVector3 translation(xOffset, yOffset, RPD_CENTER_Z_POSITION);
	transformedVertex -= translation;
	transformedVertex.Transform(rotation);

	proton.Transform(rotation);

	bool error = false;
	const double a = proton.X()*proton.X() + proton.Y()*proton.Y();
	const double b = 2 * (proton.X()*transformedVertex.X() + proton.Y()*transformedVertex.Y());
	double c = transformedVertex.X()*transformedVertex.X() + transformedVertex.Y()*transformedVertex.Y() - INNER_RING_DIAMETER*INNER_RING_DIAMETER;
	double n;
	try {
		n = __getPositiveSolutionOfQuadraticEquation(a, b, c);
	} catch (const int& number) {
		if(number != 42) {
			std::cerr<<"In antok::getRPDExpectedHitsParameters: Something went very "
			         <<"wrong when calculating the RPD hit coordinates for ring A. Aborting..."<<std::endl;
			throw;
		}
		error = true;
	}
	TVector3 coordinatesRingA = transformedVertex + n*proton;
	rpdZRingA = coordinatesRingA.Z();
	rpdPhiRingA = coordinatesRingA.Phi();
	if(error) {
		rpdZRingA = 0.;
		rpdPhiRingA = -10.;
		error = false;
	}
	c = transformedVertex.X()*transformedVertex.X() + transformedVertex.Y()*transformedVertex.Y() - OUTER_RING_DIAMETER*OUTER_RING_DIAMETER;
	try {
		n = __getPositiveSolutionOfQuadraticEquation(a, b, c);
	} catch (const int& number) {
		if(number != 42) {
			std::cerr<<"In antok::getRPDExpectedHitsParameters: Something went very "
			         <<"wrong when calculating the RPD hit coordinates for ring B. Aborting..."<<std::endl;
			throw;
		}
		error = true;
	}
	TVector3 coordinatesRingB = transformedVertex + n*proton;
	rpdZRingB = coordinatesRingB.Z();
	rpdPhiRingB = coordinatesRingB.Phi();
	if(error) {
		rpdZRingB = 0.;
		rpdPhiRingB = -10.;
	}

}
