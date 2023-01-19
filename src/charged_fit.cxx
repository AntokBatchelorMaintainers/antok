// Based on antok function ´neutral_problem.cxx´  by T. Schlüter,
// modified to charged particle tracks by M. Wagner (HISKP Bonn),
// implemented by J. Beckers.
// Measured values that enter the fit (see COMPASS note 2007-10):
// energy E_i of particle i = {1, 2} and normalized momenta (x_i, y_i, z_i) of
// particle tracks [e.g. x_i = p_x_i/abs(p_i)].
//
// The set eta of variables that enters the fit is
// eta = (x_1, y_1, E_1, x_2, y_2, E_2)
// where x_i, y_i, z_i are the direction cosines.
// in other words, (x_i, y_i, z_i) is the unit vector pointing along
// the particle direction in the lab frame


#include <sstream>
#include <vector>

#include "TDirectory.h"
#include "TH1I.h"
#include "TLorentzVector.h"

#include "charged_fit.h"


// // initialize static member variables
bool              antok::ChargedFit::_first = true;
std::vector<TH1*> antok::ChargedFit::_pullsHists;


antok::ChargedFit::ChargedFit(const TVector3& particle1Momentum,
                              const TVector3& particle2Momentum,
                              const double&   particle1Energy,
                              const double&   particle2Energy,
                              const double    particle1Mass,
                              const double    particle2Mass,
                              const std::vector<double>& particle1MomentumCovariance,
                              const std::vector<double>& particle2MomentumCovariance,
                              const double    mass,
                              const double    precisionGoal)
	: _problem                    (mass, particle1Mass, particle2Mass, precisionGoal),
	  _kinFitter                  (nullptr),
	  _particle1Momentum          (particle1Momentum),
	  _particle2Momentum          (particle2Momentum),
	  _particle1Energy            (particle1Energy),
	  _particle2Energy            (particle2Energy),
	  _particle1Mass              (particle1Mass),
	  _particle2Mass              (particle2Mass),
	  _particle1MomentumCovariance(particle1MomentumCovariance),
	  _particle2MomentumCovariance(particle2MomentumCovariance),
	  _improvedLV_1               (0, 0, 0, 0),
	  _improvedLV_2               (0, 0, 0, 0),
	  _improvedLVSum              (0, 0, 0, 0),
	  _mass                       (mass),
	  _precisionGoal              (precisionGoal)
{
	if (_first) {
		_first = false;
		initPullHists();
	}

	// Construct covariance matrix for particles
	_covMat.ResizeTo(_nmbFittedVars, _nmbFittedVars);
	// Set submatrix for first particle
	TMatrixDSub(_covMat, 0, 2, 0, 2) = covForParticle(_particle1Momentum, _particle1Energy, _particle1MomentumCovariance);
	// Set submatrix for second particle
	TMatrixDSub(_covMat, 3, 5, 3, 5) = covForParticle(_particle2Momentum, _particle2Energy, _particle2MomentumCovariance);

	// Initial measured values for the fit
	const TVector3 particleDir_1 = _particle1Momentum.Unit();
	const TVector3 particleDir_2 = _particle2Momentum.Unit();
	_startValues.ResizeTo(_nmbFittedVars);
	_startValues[0] = particleDir_1.X();
	_startValues[1] = particleDir_1.Y();
	_startValues[2] = _particle1Energy;
	_startValues[3] = particleDir_2.X();
	_startValues[4] = particleDir_2.Y();
	_startValues[5] = _particle2Energy;

	_kinFitter = new KinematicFit(_problem, _startValues, _covMat);
	assert(_kinFitter);
}


void
antok::ChargedFit::initPullHists()
{
	const std::string varNames[_nmbFittedVars] = {"x1", "y1", "E1", "x2", "y2", "E2"};
	for (size_t i = 0; i < _nmbFittedVars; ++i) {
		std::stringstream histName, histTitle;
		histName << "hPulls" << i + 1;
		histTitle << "pulls for variable " << varNames[i];
		_pullsHists.push_back(new TH1I(histName.str().c_str(), histTitle.str().c_str(), 100, -3, 3));
	}
}


TMatrixDSym
antok::ChargedFit::covForParticle(const TVector3& particleMomentum,
	                              const double&   particleEnergy,
	                              const std::vector<double>& particleMomentumCovariance)
{

	const double cov[16] = {particleMomentumCovariance [0], particleMomentumCovariance [1], particleMomentumCovariance [2], particleMomentumCovariance [3],
	                        particleMomentumCovariance [4], particleMomentumCovariance [5], particleMomentumCovariance [6], particleMomentumCovariance [7],
	                        particleMomentumCovariance [8], particleMomentumCovariance [9], particleMomentumCovariance[10], particleMomentumCovariance[11],
	                        particleMomentumCovariance[12], particleMomentumCovariance[13], particleMomentumCovariance[14], particleMomentumCovariance[15]};
    
    TMatrixDSym Cov(4, cov);

	// Linear uncertainty propagation for variable transformation (X, Y, Z, E) -> (x, y, E),
	// where x and y are the direction cosines, i.e. x = X / R, y = Y / R
	const double X = particleMomentum.X();
	const double Y = particleMomentum.Y();
	const double Z = particleMomentum.Z();
	const double X2 = X * X;
	const double Y2 = Y * Y;
	const double Z2 = Z * Z;
	const double R = sqrt(X2 + Y2 + Z2);
	const double invR3 = 1 / (R * R * R);

	// Jacobian
	const double j[12] = {(Y2 + Z2) * invR3,  -X * Y * invR3,     -X * Z * invR3,  0,
	                      -X * Y * invR3,     (X2 + Z2) * invR3,  -Y * Z * invR3,  0,
	                      0,                  0,                  0,               1};
	const TMatrixD J(3, 4, j);

	// Finally, calculate the covariance matrix for (x, y, E).
	Cov.Similarity(J);

	return Cov;
}


void
antok::ChargedFit::fillPullHists(const TVectorD& improvedMeas)
{
	TMatrixDSym covImprovedMeas(_kinFitter->covImprovedMeasurements());
	_pullValues.resize(_nmbFittedVars);
	for (size_t i = 0; i < _nmbFittedVars; ++i) {
		_pullValues[i] = (_startValues[i] - improvedMeas[i]) / sqrt(_covMat[i][i] - covImprovedMeas[i][i]);
		_pullsHists[i]->Fill(_pullValues[i]);
	}
}


bool
antok::ChargedFit::doFit() {
	// Perform the fit, returning early if it fails
	if (not _kinFitter->doFit()) {
		return false;
	}

	// Get improved values
	TVectorD improvedMeas(_kinFitter->improvedMeasurements());
	fillPullHists(improvedMeas);
	// particle 1
	const double improved_x_1 = improvedMeas[0]; // px1/p
	const double improved_y_1 = improvedMeas[1];
	const double improved_z_1 = sqrt(1 - improved_x_1 * improved_x_1 - improved_y_1 * improved_y_1);
	const double improved_E_1 = improvedMeas[2];
    const double improved_p_1 = sqrt(improved_E_1*improved_E_1 - _particle1Mass*_particle1Mass);
	// particle 2
	const double improved_x_2 = improvedMeas[3];
	const double improved_y_2 = improvedMeas[4];
	const double improved_z_2 = sqrt(1 - improved_x_2 * improved_x_2 - improved_y_2 * improved_y_2);
	const double improved_E_2 = improvedMeas[5];
    const double improved_p_2 = sqrt(improved_E_2*improved_E_2 - _particle2Mass*_particle2Mass);
	// improved particle directions
	const TVector3 improvedDir_1(improved_x_1, improved_y_1, improved_z_1);
	const TVector3 improvedDir_2(improved_x_2, improved_y_2, improved_z_2);
	// improved particle 4-momenta
	_improvedLV_1  = TLorentzVector(improvedDir_1 * improved_p_1, improved_E_1);
	_improvedLV_2  = TLorentzVector(improvedDir_2 * improved_p_2, improved_E_2);
	_improvedLVSum = _improvedLV_1 + _improvedLV_2;

	double* covMatElems = _covMat.GetMatrixArray();
	_transfCov.assign(covMatElems,covMatElems+36);

	return true;
}
