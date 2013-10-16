#include<basic_calcs.h>

#include<iostream>

#include<TLorentzVector.h>
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
