// Based on code by Tobias Schlüter  <tobias.schlueter@physik.uni-muenchen.de>,
// available at
// https://gitlab.cern.ch/compass/hadron/hadrontools/-/tree/master/kinematicFit
// and
// http://wwwcompass.cern.ch/compass/notes/2007-10/source.tar.gz
// and documented in COMPASS note 2007-10 available at
// https://wwwcompass.cern.ch/compass/notes/2007-10/2007-10.pdf


#ifndef KINEMATICFIT_H
#define KINEMATICFIT_H

#include <cassert>

#include "TMath.h"
#include "TVectorD.h"

namespace antok {

	//! A class for kinematic fitting.
	/*! Implements the kinematic fitting algorithm given in Brandt, Ch.
	    9.10. It fits a vector eta of measured values with covariance
	    matrix covEta and a vector fitPars of derived parameters under
	    the constraint encoded in the user provided problem. */
	class KinematicFit {

	public:

		//! A class to define what the fitting algorithm needs to know
		class Problem {

		public:

			virtual ~Problem() { }

			//! Return values of constraint functions
			/*! The constraint functions must evaluate to 0, if constraints are fulfilled
			    \param[in] fitPars current estimate for the fit parameters
			    \param[in] eta     current estimate for the fitted measured values
			    \return the vector of values of the constraint functions */
			virtual const TVectorD& constraintFuncs(const TVectorD& fitPars, const TVectorD& eta) { assert(0); };
			//! Return values of constraint function
			/*! As above, but for the case without fit parameters; either
			    this function or the one above needs to be implemented. */
			virtual const TVectorD& constraintFuncs(const TVectorD& eta)                          { assert(0); }

			//! Return Jacobian matrix of the constraint funtions
			/*! \param[in]  fitPars current estimate for the fit parameters
			    \param[in]  eta     current estimate for the fitted measured values
			    \param[out] jacobianFitPars matrix of the first derivatives of the constraint functions w.r.t. fitPars
			    \param[out] jacobianEta     matrix of the first derivatives of the constraint functions w.r.t. eta */
			virtual void jacobianConstraintFuncs(const TVectorD& fitPars, const TVectorD& eta,
			                                     TMatrixD& jacobianFitPars, TMatrixD& jacobianEta) { assert(0); }
			//! Return Jacobian matrix of the constraint functions
			/*! As above, but for the case without fit parameters; either
			    this function or the one above needs to be implemented. */
			virtual void jacobianConstraintFuncs(const TVectorD& eta, TMatrixD& jacobianEta)       { assert(0); }

			//! Return whether convergence criterion was reached
			virtual bool   isConverged()    const = 0;
			//! Return number of constraints
			virtual size_t nmbConstraints() const = 0;

		};

		/*! General constructor
		  \param     prob         the user-defined fit problem to be solved
		  \param[in] fitParsStart start values for the fit parameters
		  \param[in] etaStart     the measured values
		  \param[in] covEtaStart  the covariance matrix of the measured values */
		KinematicFit(Problem& prob, const TVectorD& fitParsStart, const TVectorD& etaStart, const TMatrixDSym& covEtaStart);
		//! Constructor as above, but for the case with no fit parameters
		KinematicFit(Problem& prob, const TVectorD& etaStart, const TMatrixDSym& covEtaStart);
		~KinematicFit() { }

		//! Execute a single fitting step
		void step();
		//! Perform the full fit
		/*! Returns whether the convergence criterion was reached in the
		    maximum allowed number of iterations
		    \sa converged() */
		bool doFit();

		//! Return number of steps performed so far
		size_t nmbSteps() const { return _nmbSteps; }
		//! Return maximum allowed number of iterations
		size_t maxNmbSteps() const { return _maxNmbSteps; }
		//! Set maximum allowed number of iterations
		void setMaxNmbSteps(const size_t maxNmbSteps) { _maxNmbSteps = maxNmbSteps; }
		//! Return current values of the fit parameters
		TVectorD fitPars() const { return _fitParsStart + _deltaFitPars; }
		//! Return current values of the improved measured values
		TVectorD improvedMeasurements() const { return _etaStart + _deltaEta; }
		//! Returns covariance matrix for the fit parameters
		TMatrixDSym covFitPars();
		//! Returns covariance matrix for the improved measured values
		TMatrixDSym covImprovedMeasurements();
		//! Return chi^2 value of fit
		double chi2Value() const;
		//! Return P-value of fit
		double pValue() const {
			return TMath::Prob(chi2Value(), _problem.nmbConstraints());
		}

	private:

		// helper functions
		TMatrixDSym calcGB(const TMatrixD& B) const;
		TMatrixDSym calcATGBAinv(const TMatrixD& A, const TMatrixDSym& G_B) const;

		Problem&    _problem;       //!< the user-defined fit problem to be solved
		size_t      _nmbSteps;      //!< number of steps performed so far
		size_t      _maxNmbSteps;   //!< maximum allowed number of iterations
		int         _nmbFitPars;    //!< number of fit parameters; if == 0, code pertaining to fit parameters is skipped
		TVectorD    _fitParsStart;  //!< vector that holds start values for fit paramters
		TVectorD    _deltaFitPars;  //!< vector that holds current improvement for fit paramters
		TVectorD    _etaStart;      //!< vector that holds measured values
		TMatrixDSym _covEtaStart;   //!< covariance matrix for the measured values
		TVectorD    _deltaEta;      //!< vector that holds current improvement for measured values

	};

}

#endif  // KINEMATICFIT_H
