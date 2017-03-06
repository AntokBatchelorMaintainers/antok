#include<basic_calcs.h>

#include<iostream>

#include<TLorentzVector.h>
#include<TVector3.h>

#include<constants.h>

TLorentzVector antok::getBeamEnergy(TVector3 p3_beam, const TLorentzVector& pX, const double mass_beam, const double mass_target) {


// old implementation using the taylor expansion of sqrt(p^2+m^2)
//	const TVector3 p3_Tot(pX.Vect());
//	const double E_X = pX.E();
//	const double theta = p3_Tot.Angle(p3_beam);
//	const double p3_Tot_Mag = p3_Tot.Mag();
//	const double p3_cosTheta = p3_Tot_Mag * TMath::Cos(theta);
//	const double a0 = (mass_beam * mass_beam) * p3_cosTheta;
//	const double a1 = (mass_target * E_X) - 0.5 * (pX.M2() + (mass_beam * mass_beam));
//	const double a2 = mass_target- E_X + (p3_cosTheta);
//	const double E_beam = (a1/(2*a2)) * (1 + std::sqrt(1 + ((2*a2*a0)/(a1*a1))));

	const TVector3 p3_Tot(pX.Vect());
	const double cos_theta = std::cos(p3_Tot.Angle(p3_beam));
	const double E_x = pX.E();
	const double p_x = pX.P();
	const double mass_x = pX.M();

	const double alpha = 2.0*mass_target*E_x - mass_beam*mass_beam - mass_x*mass_x;
	const double beta = E_x - mass_target;
	const double c = alpha*alpha - 4.0*beta*beta*mass_beam*mass_beam;
	const double b = -4.0*alpha*p_x*cos_theta;
	const double a = 4.0*p_x*p_x*cos_theta*cos_theta - 4.0*beta*beta;


	const double p_beam = (-b + std::sqrt(b*b-4.0*a*c))/(2.0*a); // use positive solution





	p3_beam.SetMag(p_beam);
	TLorentzVector pBeam;
	pBeam.SetXYZM(p3_beam.X(), p3_beam.Y(), p3_beam.Z(), mass_beam);
	return pBeam;

};
TLorentzVector antok::getBeamEnergy(TVector3 p3_beam, const TLorentzVector& LV_X){
	return getBeamEnergy(p3_beam, LV_X, antok::Constants::chargedPionMass(), antok::Constants::protonMass());
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
