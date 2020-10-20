// Equation numbers refer to Sec. 9.10 on p. 251ff in Siegmund Brandt,
// "Data Analysis", 4th Ed., Springer 2014, DOI
// 10.1007/978-3-319-03762-2


#include "TDecompChol.h"

#include "kinematic_fit.h"


antok::KinematicFit::KinematicFit(Problem&           prob,
                                  const TVectorD&    fitParsStart,
                                  const TVectorD&    etaStart,
                                  const TMatrixDSym& covEtaStart)
	: _problem     (prob),
	  _nmbSteps    (0),
	  _maxNmbSteps (10),
	  _fitParsStart(fitParsStart),
	  _etaStart    (etaStart),
	  _covEtaStart (covEtaStart)
{
	_deltaFitPars.ResizeTo(_fitParsStart);
	_deltaFitPars = 0;
	_nmbFitPars = _deltaFitPars.GetNrows();
	_deltaEta.ResizeTo(_etaStart);
	_deltaEta = 0;
}


antok::KinematicFit::KinematicFit(Problem&           prob,
                                  const TVectorD&    etaStart,
                                  const TMatrixDSym& covEtaStart)
	: _problem    (prob),
	  _nmbSteps   (0),
	  _maxNmbSteps(10),
	  _etaStart   (etaStart),
	  _covEtaStart(covEtaStart)
{
	_nmbFitPars = 0;
	_fitParsStart.ResizeTo(0);
	_deltaFitPars.ResizeTo(0);
	_deltaEta.ResizeTo(_etaStart);
	_deltaEta = 0;
}


//! calculates weight matrix G_B as defined in Eq. (9.10.12)
TMatrixDSym
antok::KinematicFit::calcGB(const TMatrixD& B) const
{
	TMatrixDSym C_y(_covEtaStart);           // covariance matrix C_y = G_y^{-1} for measured values
	TMatrixDSym& BCyBT = C_y.Similarity(B);  // B C_y B^T
	return BCyBT.Invert();
}


//! calculates term (A^T G_B A)^{-1}
TMatrixDSym
antok::KinematicFit::calcATGBAinv(const TMatrixD&    A,
                                  const TMatrixDSym& G_B) const
{
		TMatrixDSym ATGBA(G_B);
		ATGBA.SimilarityT(A);  // A^T G_B A
		return ATGBA.Invert();
}


void
antok::KinematicFit::step()
{
	// We follow the iterative procedure described on page 254f (unfortunately no equation numbers here)
	// at iteration number i:
	// fit parameters: xi_i = fitPars_i - fitPars_(i - 1); xi_i is calculated
	// hence, fitPars_i = _deltaFitPars_i + _fitParsStart with _deltaFitPars_i = sum_{l = 1}^i xi_i
	// improved mesurement values: delta_i = eta_i - eta_(i - 1); delta_i is calculated
	// hence, eta_i = _deltaEta_i + _etaStart with _deltaEta_i = sum_{l = 1}^i delta_i

	// Get values of contraint functions for
	//     fitPars_(i - 1) = _fitParsStart + _deltaFitPars_(i - 1)
	// and eta_(i - 1)     = _etaStart     + _deltaEta_(i - 1)
	// corresponds to vector c_i on p. 254 and to c in Eq. (9.10.5), respectively
	const TVectorD& constrVals = (_nmbFitPars == 0) ?
		_problem.constraintFuncs(_etaStart + _deltaEta) :
		_problem.constraintFuncs(_fitParsStart + _deltaFitPars, _etaStart + _deltaEta);
	// get Jacobian of contraint functions w.r.t. measured values and fit parameters
	TMatrixD JEta, JFitPars;  // correspond to A_i, B_i on p. 254 and A, B in Eqs. (9.10.3) and (9.10.4), respectively
	if (_nmbFitPars == 0) {
		_problem.jacobianConstraintFuncs(_etaStart + _deltaEta, JEta);
	} else {
		_problem.jacobianConstraintFuncs(_fitParsStart + _deltaFitPars, _etaStart + _deltaEta, JFitPars, JEta);
	}

	const TMatrixDSym G_B = calcGB(JEta);
	//TODO shouldn't we use formula at top of p. 255: delta_i = -G_y^{-1} B^T G_B(c - B _deltaEta_(i - 1) + A xi_i)?
	//     i.e. delta should be initialized with c_i - B_i _deltaEta_(i - 1)
	// calculate delta_i using Eq. (9.10.15) in the form delta_i = -G_y^{-1} B^T G_B(c + A xi_i)
	TVectorD delta(constrVals);  // c vector in Eq. (9.10.15)
	if (_nmbFitPars != 0) {
		//TODO shouldn't we use formula at bottom of p. 254: xi_i = -(A^T G_B A)^{-1} A^T G_B (c - B _deltaEta_(i - 1))?
		//     i.e. xi should be initialized with c_i - B_i _deltaEta_(i - 1)
		// calculate xi_i using Eq. (9.10.14), i.e. xi_i = -(A^T G_B A)^{-1} A^T G_B c
		const TMatrixDSym ATGBAinv    = calcATGBAinv(JFitPars, G_B);
		const TMatrixD    ATGB        = TMatrixD(JFitPars, TMatrixD::kTransposeMult, G_B);
		const TMatrixD    totalMatrix = (TMatrixD(ATGBAinv, TMatrixD::kMult, ATGB) *= -1);
		TVectorD xi(constrVals);  // c vector in Eq. (9.10.14)
		xi *= totalMatrix;
		_deltaFitPars += xi;  // accumulate correction for fit parameters, i.e. _deltaFitPars_i = _deltaFitPars_(i - 1) + xi

		TVectorD Axi(xi);
		Axi *= JFitPars;
		delta += Axi;  // (c + A xi_i)
	}
	const TMatrixD BTGB = TMatrixD(JEta, TMatrixD::kTransposeMult, G_B);
	delta *= (TMatrixD(_covEtaStart, TMatrixD::kMult, BTGB) *= -1);
	_deltaEta += delta;  // accumulate correction for improved measurements, i.e. _deltaEta_i = _deltaEta_(i - 1) + delta
}


bool
antok::KinematicFit::doFit()
{
	do {
		step();
	} while ((not _problem.isConverged()) and (++_nmbSteps < _maxNmbSteps));
	return _problem.isConverged();
}


TMatrixDSym
antok::KinematicFit::covFitPars()
{
	assert(_nmbFitPars > 0);
	// get Jacobian of contraint functions w.r.t. measured values and fit parameters
	TMatrixD JEta, JFitPars;  // A and B in Eqs. (9.10.3) and (9.10.4)
	_problem.jacobianConstraintFuncs(_fitParsStart + _deltaFitPars, _etaStart + _deltaEta, JFitPars, JEta);

	// calculate Eq. (9.10.19) C_xiTilde = (A^T G_B A)^{-1}
	const TMatrixDSym G_B = calcGB(JEta);
	return calcATGBAinv(JFitPars, G_B);
}


TMatrixDSym
antok::KinematicFit::covImprovedMeasurements()
{
	// get Jacobian of contraint functions w.r.t. measured values and fit parameters
	TMatrixD JEta, JFitPars;  // A and B in Eqs. (9.10.3) and (9.10.4)
	if (_nmbFitPars == 0) {
		_problem.jacobianConstraintFuncs(_etaStart + _deltaEta, JEta);
	} else {
		_problem.jacobianConstraintFuncs(_fitParsStart + _deltaFitPars, _etaStart + _deltaEta, JFitPars, JEta);
	}

	// calculate Eq. (9.10.20) in the form
	// C_etaTilde = C_y - Delta C_y
	// with Delta C_y = C_y B^T (G_B - matrix) B C_y
	// where matrix = G_B [A (A^T G_B A)^{-1} A^T] G_B is different from 0 only when fit parameters exist
	const TMatrixDSym G_B = calcGB(JEta);
	TMatrixDSym DeltaCy(G_B);
	if (_nmbFitPars != 0) {
		TMatrixDSym matrix(calcATGBAinv(JFitPars, G_B));
		matrix.Similarity(JFitPars);  // A (...) A^T
		matrix.Similarity(G_B);       // G_B (...) G_B^T; using G_B == G_B^T
		DeltaCy -= matrix;            // (G_B - matrix)
	}
	DeltaCy.SimilarityT(JEta);         // B^T (...) B
	DeltaCy.Similarity(_covEtaStart);  // C_y (...) C_y^T, using C_y == C_y^T
	TMatrixDSym covEtaTilde(_covEtaStart);
	covEtaTilde -= DeltaCy;            // C_y - Delta C_y
	return covEtaTilde;
}


double
antok::KinematicFit::chi2Value() const
{
	// get precision matrix for measured values
	TMatrixDSym precEtaStart(_covEtaStart);
	precEtaStart.Invert();

	return precEtaStart.Similarity(_deltaEta);
}
