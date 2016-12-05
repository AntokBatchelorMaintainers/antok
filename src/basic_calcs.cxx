#include<basic_calcs.h>

#include<iostream>

#include<TLorentzVector.h>
#include<TVector3.h>

#include<constants.h>

TLorentzVector antok::getBeamEnergy(TVector3 p3_beam, const TLorentzVector& pX, const double mass_beam) {

	const double& PROTON_MASS = antok::Constants::protonMass();

	const TVector3 p3_Tot(pX.Vect());
	const double E_X = pX.E();
	const double theta = p3_Tot.Angle(p3_beam);
	const double p3_Tot_Mag = p3_Tot.Mag();
	const double p3_cosTheta = p3_Tot_Mag * TMath::Cos(theta);
	const double a0 = (mass_beam * mass_beam) * p3_cosTheta;
	const double a1 = (PROTON_MASS * E_X) - 0.5 * (pX.M2() + (mass_beam * mass_beam));
	const double a2 = PROTON_MASS - E_X + (p3_cosTheta);
	const double E_beam = (a1/(2*a2)) * (1 + std::sqrt(1 + ((2*a2*a0)/(a1*a1))));
	const double p_beam = std::sqrt((E_beam * E_beam) - (mass_beam * mass_beam));
	p3_beam.SetMag(p_beam);
	TLorentzVector pBeam;
	pBeam.SetXYZM(p3_beam.X(), p3_beam.Y(), p3_beam.Z(), mass_beam);
	return pBeam;

};
TLorentzVector antok::getBeamEnergy(TVector3 p3_beam, const TLorentzVector& LV_X){ return getBeamEnergy(p3_beam, LV_X, antok::Constants::chargedPionMass()); }

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


double antok::utils::getPositiveSolutionOfQuadraticEquation(const double& a,
                                                            const double& b,
                                                            const double& c)
{
	double result = (-b + std::sqrt(b*b - 4*a*c)) / (2*a);
	if(result < 0) {
		result = (-b - std::sqrt(b*b - 4*a*c)) / (2*a);
	}
	if(result < 0) {
		throw 42;
	}
	return result;
}
