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


#include "charged_problem.h"
#include "neutral_problem.h"

using antok::QuantitiyHelper;


antok::ChargedProblem::ChargedProblem(const double mass,
                                      const double d1mass,
                                      const double d2mass,
                                      const double precisionGoal)
	: _mass2        (mass * mass),
	  _d1mass2      (d1mass * d1mass),
	  _d2mass2      (d2mass * d2mass),
	  _funcValue    (1),
	  _precisionGoal(precisionGoal)
{ }


// implements Eq. (5) in COMPASS note 2007-10, modified for charged particles
const TVectorD&
antok::ChargedProblem::constraintFuncs(const TVectorD& eta)
{
	const QuantitiyHelper q = QuantitiyHelper(eta);
	_funcValue[0] = 2 * ( q.E_1 * q.E_2 - std::sqrt(q.E_1*q.E_1 - _d1mass2) * std::sqrt(q.E_2*q.E_2 - _d2mass2) * ( q.x_1 * q.x_2 + q.y_1 * q.y_2 + q.z_1 * q.z_2 )) + _d1mass2 + _d2mass2 - _mass2;
	return _funcValue;
}


void
antok::ChargedProblem::jacobianConstraintFuncs(const TVectorD& eta,
                                               TMatrixD&       jacobianEta)
{
	const QuantitiyHelper q = QuantitiyHelper(eta);
	const double term = q.x_1 * q.x_2 + q.y_1 * q.y_2 + q.z_1 * q.z_2;
	const double p1 = std::sqrt(q.E_1*q.E_1 - _d1mass2);
	const double p2 = std::sqrt(q.E_2*q.E_2 - _d2mass2);
	const double J[6] = {
		2 * p1 * p2 * (-q.x_2 + q.z_2 * q.x_1 / q.z_1),  // d F / d x_1
		2 * p1 * p2 * (-q.y_2 + q.z_2 * q.y_1 / q.z_1),  // d F / d y_1
		2 * ( q.E_2 - p2 * q.E_1 / p1 * term ),          // d F / d E_1
		2 * p1 * p2 * (-q.x_1 + q.z_1 * q.x_2 / q.z_2),  // d F / d x_2
		2 * p1 * p2 * (-q.y_1 + q.z_1 * q.y_2 / q.z_2),  // d F / d y_2
		2 * ( q.E_1 - p1 * q.E_2 / p2 * term )           // d F / d E_2
	};
	jacobianEta.ResizeTo(1, 6);
	jacobianEta.SetMatrixArray(J);
}
