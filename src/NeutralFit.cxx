#include <vector>

#include "TLorentzVector.h"
#include "TH1I.h"
#include "TDirectory.h"

#include "NeutralFit.h"

#define NUMVARS 6

bool antok::NeutralFit::first = true;
std::vector<TH1 *> antok::NeutralFit::hPulls;

antok::NeutralFit::NeutralFit(const TVector3 &vertexPosition_,
                              const TVector3 &cluster1Position_,
                              const TVector3 &cluster2Position_,
                              const TVector3 &cluster1PositionError_,
                              const TVector3 &cluster2PositionError_,
                              double cluster1Energy_,
                              double cluster2Energy_,
                              double cluster1EnergyError_,
                              double cluster2EnergyError_,
                              double mass_,
                              double window_,
                              int whichDeltaE_)
		: vertexPosition(vertexPosition_),
		  cluster1Position(cluster1Position_),
		  cluster2Position(cluster2Position_),
		  cluster1PositionError(cluster1PositionError_),
		  cluster2PositionError(cluster2PositionError_),
		  cluster1Energy(cluster1Energy_),
		  cluster2Energy(cluster2Energy_),
		  cluster1EnergyError(cluster1EnergyError_),
		  cluster2EnergyError(cluster2EnergyError_),
		  lv1(0, 0, 0, 0),
		  lv2(0, 0, 0, 0),
		  lvSum(0, 0, 0, 0),
		  mass(mass_),
		  window(window_),
		  whichDeltaE(whichDeltaE_) {
	if (first) {
		first = false;
		initPulls();
	}

	covMat.ResizeTo(NUMVARS, NUMVARS);
	TMatrixDSub(covMat, 0, 2, 0, 2) = covMatForCluster(cluster1Position, cluster1PositionError, cluster1Energy, cluster1EnergyError);
	TMatrixDSub(covMat, 3, 5, 3, 5) = covMatForCluster(cluster2Position, cluster2PositionError, cluster2Energy, cluster2EnergyError);

	// Initial values for the fit.
	const double x1 = cluster1Position.X() - vertexPosition.X();
	const double y1 = cluster1Position.Y() - vertexPosition.Y();
	const double z1 = cluster1Position.Z() - vertexPosition.Z();
	const double r1 = hypot(hypot(x1, y1), z1);

	const double x2 = cluster2Position.X() - vertexPosition.X();
	const double y2 = cluster2Position.Y() - vertexPosition.Y();
	const double z2 = cluster2Position.Z() - vertexPosition.Z();
	const double r2 = hypot(hypot(x2, y2), z2);

	startingValues.ResizeTo(NUMVARS);
	startingValues[0] = x1 / r1;
	startingValues[1] = y1 / r1;
	startingValues[2] = cluster1Energy;
	startingValues[3] = x2 / r2;
	startingValues[4] = y2 / r2;
	startingValues[5] = cluster1Energy;

	// Setup the worker objects.
	myProblem = new NeutralProblem(mass);
	myFitter = new KinematicFit(*myProblem, startingValues, covMat);
}


bool antok::NeutralFit::isInWindow() const {
	TVector3 v(vertexPosition.X(), vertexPosition.Y(), vertexPosition.Z());

	TVector3 a3(cluster1Position.X(), cluster1Position.Y(), cluster1Position.Z());
	TLorentzVector a((a3 - v).Unit() * cluster1Energy, cluster1Energy);

	TVector3 b3(cluster2Position.X(), cluster2Position.Y(), cluster2Position.Z());
	TLorentzVector b((b3 - v).Unit() * cluster2Energy, cluster2Energy);

	return fabs((a + b).M() - mass) < window * mass;
}


void antok::NeutralFit::initPulls() {
	const char *varname[NUMVARS] = {"x1", "y1", "E1",
	                                "x2", "y2", "E2",};
	for (int i = 0; i < NUMVARS; i++) {
		char name[20];
		char title[50];

		sprintf(name, "hPulls%d", i + 1);
		sprintf(title, "pulls for variable %s", varname[i]);
		hPulls.push_back(new TH1I(name, title, 100, -3, 3));
	}

}


TMatrixDSym antok::NeutralFit::covMatForCluster(const TVector3 &clusterPosition,
                                                const TVector3 &clusterPositionError,
                                                const double clusterEnergy,
                                                const double clusterEnergyError ) {
	double sigE2;
	switch (whichDeltaE) {
		case 0:
			sigE2 = clusterEnergyError;
			break;
		case 1:
			sigE2 = TMath::Power(clusterEnergy,2);
			break;
		case 2:
			sigE2 = 0.055 * 0.055 * clusterEnergy + 0.015 * 0.015 * TMath::Power(clusterEnergy,2);
			break;
		default:
			abort();
	}

	double c[16] = {clusterPositionError.X(),                        0,                        0,       0,
	                                       0, clusterPositionError.Y(),                        0,       0,
	                                       0,                        0, clusterPositionError.Z(),       0,
	                                       0,                        0,                        0,  sigE2 };
	TMatrixDSym C(4, c);

	/* Linearized change of variables (X, Y, Z, E) -> (x, y, E).
	   x,y are the direction cosines, i.e. x = X / R, y = Y / R  */
	double X = clusterPosition.X() - vertexPosition.X();
	double Y = clusterPosition.Y() - vertexPosition.Y();
	double Z = clusterPosition.Z() - vertexPosition.Z();
	double X2 = X * X, Y2 = Y * Y, Z2 = Z * Z;
	double R = hypot(hypot(X, Y), Z);
	double R3 = 1 / (R * R * R);

	double a[] =
			{(Y2 + Z2) * R3, -X * Y * R3, -X * Z * R3, 0,
			 -X * Y * R3, (X2 + Z2) * R3, -Y * Z * R3, 0,
			 0, 0, 0, 1,};
	TMatrixD A(3, 4, a);

	// Finally, calculate the covariance matrix for (x, y, E).
	C.Similarity(A);

	return C;
}

void antok::NeutralFit::fillPulls( TVectorD enhanced ) {
	TMatrixDSym covEnh(myFitter->getCovEnhanced());
	pulls.resize(NUMVARS);
	for (int i = 0; i < NUMVARS; i++) {
		pulls[i] = (startingValues[i] - enhanced[i])
		           / sqrt(covMat[i][i] - covEnh[i][i]);
		hPulls[i]->Fill(pulls[i]);
	}
}

bool antok::NeutralFit::doFit() {
	// Do the fit, returning early if it fails.
	if (!myFitter->doFit())
	{
		return false;
	}

	TVectorD enhanced(myFitter->getEnhanced());

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

	lv1 = TLorentzVector(p1 * newE1, newE1);
	lv2 = TLorentzVector(p2 * newE2, newE2);
	lvSum = lv1 + lv2;

	fillPulls(enhanced);

	return true;
}
