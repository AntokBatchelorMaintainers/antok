#include "assert.h"

#include "TVectorD.h"
#include "TDecompChol.h"

#include "kinematic_fit.h"


antok::KinematicFit::KinematicFit(problem &prob,
                                  const TVectorD &x_,
                                  const TVectorD &eta_,
                                  const TMatrixDSym &coveta_)
		: nSteps(0),
		  maxSteps(10),
		  x(x_),
		  eta(eta_),
		  coveta(coveta_),
		  myProblem(prob) {
	dx.ResizeTo(x);
	dx = 0;
	numParams = dx.GetNrows();

	deta.ResizeTo(eta);
	deta = 0;
}

antok::KinematicFit::KinematicFit(problem &prob,
                                  const TVectorD &eta_,
                                  const TMatrixDSym &coveta_)
		: nSteps(0),
		  maxSteps(10),
		  eta(eta_),
		  coveta(coveta_),
		  myProblem(prob) {
	x.ResizeTo(0);
	dx.ResizeTo(0);
	numParams = 0;

	deta.ResizeTo(eta);
	deta = 0;
}


void antok::KinematicFit::step() {
	const TVectorD *pc;
	if (numParams == 0) {
		pc = &myProblem.constraint(eta + deta);
	} else {
		pc = &myProblem.constraint(x + dx, eta + deta);
	}
	const TVectorD &c(*pc);

	TMatrixD A, B;
	if (numParams == 0) {
		myProblem.dConstraint(eta + deta, B);
	} else {
		myProblem.dConstraint(x + dx, eta + deta, A, B);
	}

	TMatrixDSym covB(coveta);
	covB.Similarity(B);

	TMatrixDSym GB(TDecompChol(covB).Invert());
	TVectorD delta(c);

	if (numParams != 0) {
		TVectorD xi(c);
		TMatrixDSym ATGBA(GB);
		ATGBA.SimilarityT(A);

		TMatrixDSym ATGBAinv(TDecompChol(ATGBA).Invert());

		xi *= TMatrixD(ATGBAinv, TMatrixD::kMult, TMatrixD(A, TMatrixD::kTransposeMult, GB));
		xi *= -1;

		dx += xi;
		TVectorD Axi(xi);
		Axi *= A;

		delta += Axi;
	}

	delta *= TMatrixD(coveta, TMatrixD::kMult, TMatrixD(B, TMatrixD::kTransposeMult, GB));
	deta -= delta;
}


bool antok::KinematicFit::doFit() {
	do {
		step();
	} while (!myProblem.converged() && ++nSteps < maxSteps);

	return myProblem.converged();
}


TMatrixDSym antok::KinematicFit::getCovParams() {
	assert (numParams > 0);

	TMatrixD A, B;
	myProblem.dConstraint(x + dx, eta + deta, A, B);

	TMatrixDSym covB(coveta);
	covB.Similarity(B);

	TMatrixDSym ATGBA(TDecompChol(covB).Invert());
	ATGBA.SimilarityT(A);

	return TMatrixDSym(TDecompChol(ATGBA).Invert());
}


TMatrixDSym antok::KinematicFit::getCovEnhanced() {
	TMatrixD A, B;
	if (numParams == 0) {
		myProblem.dConstraint(eta + deta, B);
	} else {
		myProblem.dConstraint(x + dx, eta + deta, A, B);
	}

	TMatrixDSym covB(coveta);
	covB.Similarity(B);
	TMatrixDSym GB(TDecompChol(covB).Invert());

	TMatrixDSym shift(GB);

	if (numParams != 0) {
		TMatrixDSym ATGBA(GB);
		ATGBA.SimilarityT(A);

		TMatrixDSym shiftA(TDecompChol(ATGBA).Invert());
		shiftA.Similarity(A);
		shiftA.Similarity(GB);

		shift -= shiftA;
	}

	shift.SimilarityT(B);
	shift.Similarity(coveta);

	TMatrixDSym ret(coveta);
	ret -= shift;

	return ret;
}


Double_t antok::KinematicFit::getChi2() {
	TMatrixDSym G(coveta);
	G.Invert();

	return G.Similarity(deta);
}
