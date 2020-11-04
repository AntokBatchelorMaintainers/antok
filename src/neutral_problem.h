// Based on code by Tobias Schlüter  <tobias.schlueter@physik.uni-muenchen.de>,
// available at
// https://gitlab.cern.ch/compass/hadron/hadrontools/-/tree/master/kinematicFit
// and
// http://wwwcompass.cern.ch/compass/notes/2007-10/source.tar.gz
// and documented in COMPASS note 2007-10 available at
// https://wwwcompass.cern.ch/compass/notes/2007-10/2007-10.pdf


#ifndef NEUTRALPROBLEM_H
#define NEUTRALPROBLEM_H

#include "TVectorD.h"

#include "kinematic_fit.h"

namespace antok {

	/*! Class that defines constraint functions and other parameters for
	    fitting two neutral clusters to a given mass */
	class NeutralProblem : public KinematicFit::Problem {

	public:

		//! \param[in] mass is the desired mass for the reconstructed photon pair
		NeutralProblem(const double mass, const double precisionGoal);
		~NeutralProblem() {}

		//! Constraint function; evaluates to 0, if constraints are fulfilled
		/*! \param[in] eta current estimate for the fitted measured values
		    \return the value of the constraint function */
		const TVectorD& constraintFuncs(const TVectorD& eta);
		//! Jacobian matrix of the constraint funtion
		/*! \param[in] eta current estimate for the fitted measured values
		    \return    matrix of the first derivatives of the constraint function w.r.t. eta */
		void jacobianConstraintFuncs   (const TVectorD& eta, TMatrixD& jacobianEta);

		//! Convergence criterion
		bool   isConverged()    const { return fabs(_funcValue[0]) / _mass2 < _precisionGoal; }
		//! There is one constraint function
		size_t nmbConstraints() const { return 1; }

	private:

		const double _mass2;            //!< mass squared that the photon pair is constrained to
		TVectorD     _funcValue;        //!< holds value of constraint function
		const double _precisionGoal;    //!< sets the limit for which fit is accepted as converged

	};

}

#endif  // NEUTRALPROBLEM_H
