// Measured values that enter the fit (see COMPASS note 2007-10):
// energy E_i of photon i = {1, 2} and position (X_i, Y_i, Z_i) of
// photon ECAL cluster relative to vertex that the photons originate
// from.
//
// The set eta of variables that enters the fit is
// eta = (x_1, y_1, E_1, x_2, y_2, E_2)
// where x_i, y_i, z_i are the direction cosines, i.e.
// x_i = X_i / R_i, y_i = Y_i / R_i, and z_i = Z_i / R_i = sqrt(1 - x_i^2 - y_i^2)
// with R_i^2 = X_i^2 + Y_i^2 + Z_i^2
// in other words, (x_i, y_i, z_i) is the unit vector pointing along
// the photon direction in the lab frame


#include <sstream>
#include <vector>

#include "TDirectory.h"
#include "TH1I.h"
#include "TLorentzVector.h"

#include "neutral_fit.h"


// // initialize static member variables
bool              antok::NeutralFit::_first = true;
std::vector<TH1*> antok::NeutralFit::_pullsHists;


antok::NeutralFit::NeutralFit(const TVector3& vertexPosition,
                              const TVector3& cluster1Position,
                              const TVector3& cluster2Position,
                              const TVector3& cluster1PositionVariance,
                              const TVector3& cluster2PositionVariance,
                              const double    cluster1Energy,
                              const double    cluster2Energy,
                              const double    cluster1EnergyVariance,
                              const double    cluster2EnergyVariance,
                              const double    mass,
                              const double    massLowerLimit,
                              const double    massUpperLimit,
                              const double    precisionGoal,
                              const int       whichEnergyVariance)
	: _problem                 (mass, precisionGoal),
	  _kinFitter               (nullptr),
	  _vertexPosition          (vertexPosition),
	  _cluster1Position        (cluster1Position),
	  _cluster2Position        (cluster2Position),
	  _cluster1PositionVariance(cluster1PositionVariance),
	  _cluster2PositionVariance(cluster2PositionVariance),
	  _cluster1Energy          (cluster1Energy),
	  _cluster2Energy          (cluster2Energy),
	  _cluster1EnergyVariance  (cluster1EnergyVariance),
	  _cluster2EnergyVariance  (cluster2EnergyVariance),
	  _improvedLV_1            (0, 0, 0, 0),
	  _improvedLV_2            (0, 0, 0, 0),
	  _improvedLVSum           (0, 0, 0, 0),
	  _mass                    (mass),
	  _massLowerLimit          (massLowerLimit),
	  _massUpperLimit          (massUpperLimit),
	  _precisionGoal           (precisionGoal),
	  _whichEnergyVariance     (whichEnergyVariance)
{
	if (_first) {
		_first = false;
		initPullHists();
	}

	// Construct covariance matrix for clusters
	_covMat.ResizeTo(_nmbFittedVars, _nmbFittedVars);
	// Set submatrix for first cluster
	TMatrixDSub(_covMat, 0, 2, 0, 2) = covForCluster(_cluster1Position, _cluster1PositionVariance, _cluster1Energy, _cluster1EnergyVariance);
	// Set submatrix for second cluster
	TMatrixDSub(_covMat, 3, 5, 3, 5) = covForCluster(_cluster2Position, _cluster2PositionVariance, _cluster2Energy, _cluster2EnergyVariance);

	// Initial measured values for the fit
	const double X_1 = _cluster1Position.X() - _vertexPosition.X();
	const double Y_1 = _cluster1Position.Y() - _vertexPosition.Y();
	const double Z_1 = _cluster1Position.Z() - _vertexPosition.Z();
	const double R_1 = sqrt(X_1 * X_1 + Y_1 * Y_1 + Z_1 * Z_1);
	const double X_2 = _cluster2Position.X() - _vertexPosition.X();
	const double Y_2 = _cluster2Position.Y() - _vertexPosition.Y();
	const double Z_2 = _cluster2Position.Z() - _vertexPosition.Z();
	const double R_2 = sqrt(X_2 * X_2 + Y_2 * Y_2 + Z_2 * Z_2);
	_startValues.ResizeTo(_nmbFittedVars);
	_startValues[0] = X_1 / R_1;
	_startValues[1] = Y_1 / R_1;
	_startValues[2] = _cluster1Energy;
	_startValues[3] = X_2 / R_2;
	_startValues[4] = Y_2 / R_2;
	_startValues[5] = _cluster2Energy;

	_kinFitter = new KinematicFit(_problem, _startValues, _covMat);
	assert(_kinFitter);
}


bool
antok::NeutralFit::massIsInWindow() const
{
	const TLorentzVector photonLV_1((_cluster1Position - _vertexPosition).Unit() * _cluster1Energy, _cluster1Energy);
	const TLorentzVector photonLV_2((_cluster2Position - _vertexPosition).Unit() * _cluster2Energy, _cluster2Energy);

	return ( ((photonLV_1 + photonLV_2).M() > _massLowerLimit) && ((photonLV_1 + photonLV_2).M() < _massUpperLimit) );
}


void
antok::NeutralFit::initPullHists()
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
antok::NeutralFit::covForCluster(const TVector3& clusterPosition,
	                               const TVector3& clusterPositionVariance,
	                               const double    clusterEnergy,
	                               const double    clusterEnergyVariance)
{
	double varE;
	switch (_whichEnergyVariance) {
		case 0:
			varE = clusterEnergyVariance;
			break;
		case 1:
			varE = clusterEnergy * clusterEnergy;
			break;
		case 2: {
			// see Sec. 6.3.3., p. 491 in [NIMA 577 (2007) 455]
			// !!!note: this is for ECAL2 only
			const double a = 0.055;
			const double b = 0.015;
			varE = a * a * clusterEnergy + b * b * clusterEnergy * clusterEnergy;
			break;
		}
		default:
			abort();
	}

	const double cov[16] = {clusterPositionVariance.X(),  0,                            0,                            0,
	                        0,                            clusterPositionVariance.Y(),  0,                            0,
	                        0,                            0,                            clusterPositionVariance.Z(),  0,
	                        0,                            0,                            0,                            varE};
  TMatrixDSym Cov(4, cov);

	// Linear uncertainty propagation for variable transformation (X, Y, Z, E) -> (x, y, E),
	// where x and y are the direction cosines, i.e. x = X / R, y = Y / R
	const double X = clusterPosition.X() - _vertexPosition.X();
	const double Y = clusterPosition.Y() - _vertexPosition.Y();
	const double Z = clusterPosition.Z() - _vertexPosition.Z();
	const double X2 = X * X;
	const double Y2 = Y * Y;
	const double Z2 = Z * Z;
	const double R = sqrt(X2 + Y2 + Z2);
	const double R3 = 1 / (R * R * R);
	// Jacobian
	const double j[12] = {(Y2 + Z2) * R3,  -X * Y * R3,     -X * Z * R3,  0,
	                      -X * Y * R3,     (X2 + Z2) * R3,  -Y * Z * R3,  0,
	                      0,               0,               0,            1};
	const TMatrixD J(3, 4, j);
	// Finally, calculate the covariance matrix for (x, y, E).
	Cov.Similarity(J);
	return Cov;
}


void
antok::NeutralFit::fillPullHists(const TVectorD& improvedMeas)
{
	TMatrixDSym covImprovedMeas(_kinFitter->covImprovedMeasurements());
	_pullValues.resize(_nmbFittedVars);
	for (size_t i = 0; i < _nmbFittedVars; ++i) {
		_pullValues[i] = (_startValues[i] - improvedMeas[i]) / sqrt(_covMat[i][i] - covImprovedMeas[i][i]);
		_pullsHists[i]->Fill(_pullValues[i]);
	}
}


bool
antok::NeutralFit::doFit() {
	// Perform the fit, returning early if it fails
	if (not _kinFitter->doFit()) {
		return false;
	}

	// Get improved values
	TVectorD improvedMeas(_kinFitter->improvedMeasurements());
	fillPullHists(improvedMeas);
	// cluster 1
	const double improved_x_1 = improvedMeas[0];
	const double improved_y_1 = improvedMeas[1];
	const double improved_z_1 = sqrt(1 - improved_x_1 * improved_x_1 - improved_y_1 * improved_y_1);
	const double improved_E_1 = improvedMeas[2];
	// cluster 2
	const double improved_x_2 = improvedMeas[3];
	const double improved_y_2 = improvedMeas[4];
	const double improved_z_2 = sqrt(1 - improved_x_2 * improved_x_2 - improved_y_2 * improved_y_2);
	const double improved_E_2 = improvedMeas[5];
	// improved photon directions
	const TVector3 improvedDir_1(improved_x_1, improved_y_1, improved_z_1);
	const TVector3 improvedDir_2(improved_x_2, improved_y_2, improved_z_2);
	// improved photon 4-momenta
	_improvedLV_1  = TLorentzVector(improvedDir_1 * improved_E_1, improved_E_1);
	_improvedLV_2  = TLorentzVector(improvedDir_2 * improved_E_2, improved_E_2);
	_improvedLVSum = _improvedLV_1 + _improvedLV_2;

	return true;
}
