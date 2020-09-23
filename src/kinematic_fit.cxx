#include "TDecompChol.h"

#include "kinematic_fit.h"


antok::KinematicFit::KinematicFit(problem&           prob,
                                  const TVectorD&    x,
                                  const TVectorD&    eta,
                                  const TMatrixDSym& coveta)
	: _nSteps(0),
	  _maxSteps(10),
	  _x(x),
	  _eta(eta),
	  _coveta(coveta),
	  _myProblem(prob)
{
	_dx.ResizeTo(_x);
	_dx = 0;
	_nParams = _dx.GetNrows();
	_deta.ResizeTo(_eta);
	_deta = 0;
}


antok::KinematicFit::KinematicFit(problem&           prob,
                                  const TVectorD&    eta,
                                  const TMatrixDSym& coveta)
	: _nSteps(0),
	  _maxSteps(10),
	  _eta(eta),
	  _coveta(coveta),
	  _myProblem(prob)
{
	_nParams = 0;
	_x.ResizeTo(0);
	_dx.ResizeTo(0);
	_deta.ResizeTo(_eta);
	_deta = 0;
}


void
antok::KinematicFit::step()
{
	//TODO this block could be simplified using conditional assignment
	const TVectorD* pc;
	if (_nParams == 0) {
		pc = &_myProblem.constraint(_eta + _deta);
	} else {
		pc = &_myProblem.constraint(_x + _dx, _eta + _deta);
	}
	const TVectorD& c(*pc);

	TMatrixD A, B;
	if (_nParams == 0) {
		_myProblem.dConstraint(_eta + _deta, B);
	} else {
		_myProblem.dConstraint(_x + _dx, _eta + _deta, A, B);
	}

	TMatrixDSym covB(_coveta);
	covB.Similarity(B);

	TMatrixDSym GB(TDecompChol(covB).Invert());
	TVectorD delta(c);

	if (_nParams != 0) {
		TVectorD xi(c);
		TMatrixDSym ATGBA(GB);
		ATGBA.SimilarityT(A);

		TMatrixDSym ATGBAinv(TDecompChol(ATGBA).Invert());

		xi *= TMatrixD(ATGBAinv, TMatrixD::kMult, TMatrixD(A, TMatrixD::kTransposeMult, GB));
		xi *= -1;

		_dx += xi;
		TVectorD Axi(xi);
		Axi *= A;

		delta += Axi;
	}

	delta *= TMatrixD(_coveta, TMatrixD::kMult, TMatrixD(B, TMatrixD::kTransposeMult, GB));
	_deta -= delta;
}


bool
antok::KinematicFit::doFit()
{
	do {
		step();
	} while ((not _myProblem.converged()) and (++_nSteps < _maxSteps));
	return _myProblem.converged();  //TODO isn't this always true?
}


TMatrixDSym
antok::KinematicFit::getCovParams()
{
	assert(_nParams > 0);

	TMatrixD A, B;
	_myProblem.dConstraint(_x + _dx, _eta + _deta, A, B);

	TMatrixDSym covB(_coveta);
	covB.Similarity(B);

	TMatrixDSym ATGBA(TDecompChol(covB).Invert());
	ATGBA.SimilarityT(A);

	return TMatrixDSym(TDecompChol(ATGBA).Invert());
}


TMatrixDSym
antok::KinematicFit::getCovEnhanced()
{
	TMatrixD A, B;
	if (_nParams == 0) {
		_myProblem.dConstraint(_eta + _deta, B);
	} else {
		_myProblem.dConstraint(_x + _dx, _eta + _deta, A, B);
	}

	TMatrixDSym covB(_coveta);
	covB.Similarity(B);
	TMatrixDSym GB(TDecompChol(covB).Invert());

	TMatrixDSym shift(GB);

	if (_nParams != 0) {
		TMatrixDSym ATGBA(GB);
		ATGBA.SimilarityT(A);

		TMatrixDSym shiftA(TDecompChol(ATGBA).Invert());
		shiftA.Similarity(A);
		shiftA.Similarity(GB);

		shift -= shiftA;
	}

	shift.SimilarityT(B);
	shift.Similarity(_coveta);

	TMatrixDSym ret(_coveta);
	ret -= shift;

	return ret;
}


Double_t
antok::KinematicFit::getChi2() const
{
	TMatrixDSym G(_coveta);
	G.Invert();

	return G.Similarity(_deta);
}
