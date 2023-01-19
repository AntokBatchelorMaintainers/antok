// Based on antok function ´neutral_problem.h´,
// originally by Tobias Schlüter  <tobias.schlueter@physik.uni-muenchen.de>,
// available at
// https://gitlab.cern.ch/compass/hadron/hadrontools/-/tree/master/kinematicFit
// and
// http://wwwcompass.cern.ch/compass/notes/2007-10/source.tar.gz
// and documented in COMPASS note 2007-10 available at
// https://wwwcompass.cern.ch/compass/notes/2007-10/2007-10.pdf


#ifndef CHARGEDPROBLEM_H
#define CHARGEDPROBLEM_H

#include "TVectorD.h"

#include "kinematic_fit.h"

namespace antok {

	/*! Class that defines constraint functions and other parameters for
	    fitting two charged tracks to a given mass */
	class ChargedProblem : public KinematicFit::Problem {

	public:

		//! \param[in] mass is the desired mass for the particle pair
		//! \param[in] d1mass is the mass of the first daughter particle
		//! \param[in] d2mass is the mass of the second daughter particle
		ChargedProblem(const double mass, const double d1mass, const double d2mass, const double precisionGoal);
		~ChargedProblem() {}

		//! Constraint function; evaluates to 0, if constraints are fulfilled
		/*! \param[in] eta current estimate for the fitted measured values
		    \return the value of the constraint function */
		const TVectorD& constraintFuncs(const TVectorD& eta);
		//! Jacobian matrix of the constraint funtion
		/*! \param[in] eta current estimate for the fitted measured values
		    \return    matrix of the first derivatives of the constraint function w.r.t. eta */
		void jacobianConstraintFuncs(const TVectorD& eta, TMatrixD& jacobianEta);

		//! Convergence criterion
		bool   isConverged()    const { return fabs(_funcValue[0]) / _mass2 < _precisionGoal; }
		//! There is one constraint function
		size_t nmbConstraints() const { return 1; }

	private:

		const double _mass2;          //!< mass squared that the particle pair is constrained to
		const double _d1mass2;        //!< mass squared of the first daughter particle
		const double _d2mass2;        //!< mass squared of the second daughter particle
		TVectorD     _funcValue;      //!< holds value of constraint function
		const double _precisionGoal;  //!< sets the limit for which fit is accepted as converged

	};

}

#endif  // CHARGEDPROBLEM_H
