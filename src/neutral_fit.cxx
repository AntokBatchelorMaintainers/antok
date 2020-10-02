#include <vector>

#include "TDirectory.h"
#include "TH1I.h"
#include "TLorentzVector.h"

#include "neutral_fit.h"

#define NUMVARS 6


bool              antok::NeutralFit::_first = true;
std::vector<TH1*> antok::NeutralFit::_hPulls;


antok::NeutralFit::NeutralFit(const TVector3& vertexPosition,
                              const TVector3& cluster1Position,
                              const TVector3& cluster2Position,
                              const TVector3& cluster1PositionError,
                              const TVector3& cluster2PositionError,
                              const double    cluster1Energy,
                              const double    cluster2Energy,
                              const double    cluster1EnergyError,
                              const double    cluster2EnergyError,
                              const double    mass,
                              const double    window,
                              const int       whichDeltaE)
	: _vertexPosition(vertexPosition),
	  _cluster1Position(cluster1Position),
	  _cluster2Position(cluster2Position),
	  _cluster1PositionError(cluster1PositionError),
	  _cluster2PositionError(cluster2PositionError),
	  _cluster1Energy(cluster1Energy),
	  _cluster2Energy(cluster2Energy),
	  _cluster1EnergyError(cluster1EnergyError),
	  _cluster2EnergyError(cluster2EnergyError),
	  _lv1(0, 0, 0, 0),
	  _lv2(0, 0, 0, 0),
	  _lvSum(0, 0, 0, 0),
	  _mass(mass),
	  _window(window),
	  _whichDeltaE(whichDeltaE)
{
	if (_first) {
		_first = false;
		initPulls();
	}

	_covMat.ResizeTo(NUMVARS, NUMVARS);
	TMatrixDSub(_covMat, 0, 2, 0, 2) = covMatForCluster(_cluster1Position, _cluster1PositionError, _cluster1Energy, _cluster1EnergyError);
	TMatrixDSub(_covMat, 3, 5, 3, 5) = covMatForCluster(_cluster2Position, _cluster2PositionError, _cluster2Energy, _cluster2EnergyError);

	// Initial values for the fit.
	const double x1 = _cluster1Position.X() - _vertexPosition.X();
	const double y1 = _cluster1Position.Y() - _vertexPosition.Y();
	const double z1 = _cluster1Position.Z() - _vertexPosition.Z();
	const double r1 = hypot(hypot(x1, y1), z1);

	const double x2 = _cluster2Position.X() - _vertexPosition.X();
	const double y2 = _cluster2Position.Y() - _vertexPosition.Y();
	const double z2 = _cluster2Position.Z() - _vertexPosition.Z();
	const double r2 = hypot(hypot(x2, y2), z2);

	_startingValues.ResizeTo(NUMVARS);
	_startingValues[0] = x1 / r1;
	_startingValues[1] = y1 / r1;
	_startingValues[2] = _cluster1Energy;
	_startingValues[3] = x2 / r2;
	_startingValues[4] = y2 / r2;
	_startingValues[5] = _cluster2Energy;

	// Setup the worker objects.
	_myProblem = new NeutralProblem(_mass);
	_myFitter  = new KinematicFit(*_myProblem, _startingValues, _covMat);
}


bool
antok::NeutralFit::isInWindow() const
{
	//TODO why are these copies needed? couldn't one use _vertexPosition etal. directly?
	TVector3 v(_vertexPosition.X(), _vertexPosition.Y(), _vertexPosition.Z());

	TVector3 a3(_cluster1Position.X(), _cluster1Position.Y(), _cluster1Position.Z());
	TLorentzVector a((a3 - v).Unit() * _cluster1Energy, _cluster1Energy);

	TVector3 b3(_cluster2Position.X(), _cluster2Position.Y(), _cluster2Position.Z());
	TLorentzVector b((b3 - v).Unit() * _cluster2Energy, _cluster2Energy);

	return fabs((a + b).M() - _mass) < _window * _mass;
}


void
antok::NeutralFit::initPulls()
{
	const char* varnames[NUMVARS] = {"x1", "y1", "E1",
	                                 "x2", "y2", "E2",};
	for (size_t i = 0; i < NUMVARS; ++i) {
		char name[20];
		char title[50];
		sprintf(name, "hPulls%lu", i + 1);
		sprintf(title, "pulls for variable %s", varnames[i]);
		_hPulls.push_back(new TH1I(name, title, 100, -3, 3));
	}
}


TMatrixDSym
antok::NeutralFit::covMatForCluster(const TVector3& clusterPosition,
	                                  const TVector3& clusterPositionError,
	                                  const double    clusterEnergy,
	                                  const double    clusterEnergyError)
{
	double sigE2;
	switch (_whichDeltaE) {
		case 0:
			sigE2 = clusterEnergyError;
			break;
		case 1:
			sigE2 = TMath::Power(clusterEnergy, 2);
			break;
		case 2:
			sigE2 = 0.055 * 0.055 * clusterEnergy + 0.015 * 0.015 * TMath::Power(clusterEnergy, 2);
			break;
		default:
			abort();
	}

	const double c[16] = {clusterPositionError.X(),  0,                         0,                         0,
	                      0,                         clusterPositionError.Y(),  0,                         0,
	                      0,                         0,                         clusterPositionError.Z(),  0,
	                      0,                         0,                         0,                         sigE2};
  TMatrixDSym C(4, c);

	// Linearized change of variables (X, Y, Z, E) -> (x, y, E).
	// x,y are the direction cosines, i.e. x = X / R, y = Y / R
	const double X = clusterPosition.X() - _vertexPosition.X();
	const double Y = clusterPosition.Y() - _vertexPosition.Y();
	const double Z = clusterPosition.Z() - _vertexPosition.Z();
	const double X2 = X * X, Y2 = Y * Y, Z2 = Z * Z;
	const double R = hypot(hypot(X, Y), Z);
	const double R3 = 1 / (R * R * R);

	const double a[12] = {(Y2 + Z2) * R3,  -X * Y * R3,     -X * Z * R3,  0,
	                      -X * Y * R3,     (X2 + Z2) * R3,  -Y * Z * R3,  0,
	                      0,               0,               0,            1};
	const TMatrixD A(3, 4, a);

	// Finally, calculate the covariance matrix for (x, y, E).
	C.Similarity(A);

	return C;
}


void
antok::NeutralFit::fillPulls(const TVectorD& enhanced)
{
	TMatrixDSym covEnh(_myFitter->getCovEnhanced());
	_pulls.resize(NUMVARS);
	for (size_t i = 0; i < NUMVARS; ++i) {
		_pulls[i] = (_startingValues[i] - enhanced[i]) / sqrt(_covMat[i][i] - covEnh[i][i]);
		_hPulls[i]->Fill(_pulls[i]);
	}
}

bool antok::NeutralFit::doFit() {
	// Do the fit, returning early if it fails.
	if (not _myFitter->doFit()) {
		return false;
	}

	TVectorD enhanced(_myFitter->getEnhanced());

	const double newx1 = enhanced[0];
	const double newy1 = enhanced[1];
	const double newz1 = sqrt(1 - newx1 * newx1 - newy1 * newy1);
	const double newE1 = enhanced[2];

	const double newx2 = enhanced[3];
	const double newy2 = enhanced[4];
	const double newz2 = sqrt(1 - newx2 * newx2 - newy2 * newy2);
	const double newE2 = enhanced[5];

	TVector3 p1(newx1, newy1, newz1);
	TVector3 p2(newx2, newy2, newz2);

	_lv1 = TLorentzVector(p1 * newE1, newE1);
	_lv2 = TLorentzVector(p2 * newE2, newE2);
	_lvSum = _lv1 + _lv2;

	fillPulls(enhanced);

	return true;
}
