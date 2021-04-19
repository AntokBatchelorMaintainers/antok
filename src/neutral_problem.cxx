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


#include "neutral_problem.h"


antok::NeutralProblem::NeutralProblem(const double mass,
                                      const double precisionGoal)
	: _mass2        (mass * mass),
	  _funcValue    (1),
	  _precisionGoal(precisionGoal)
{ }


// implements Eq. (5) in COMPASS note 2007-10
const TVectorD&
antok::NeutralProblem::constraintFuncs(const TVectorD& eta)
{
	const QuantitiyHelper q = QuantitiyHelper(eta);
	_funcValue[0] = 2 * q.E_1 * q.E_2 * (1 - q.x_1 * q.x_2 - q.y_1 * q.y_2 - q.z_1 * q.z_2) - _mass2;
	return _funcValue;
}


void
antok::NeutralProblem::jacobianConstraintFuncs(const TVectorD& eta,
                                               TMatrixD&       jacobianEta)
{
	const QuantitiyHelper q = QuantitiyHelper(eta);
	const double term = 1 - q.x_1 * q.x_2 - q.y_1 * q.y_2 - q.z_1 * q.z_2;
	const double J[6] = {
		2 * q.E_1 * q.E_2 * (-q.x_2 + q.z_2 * q.x_1 / q.z_1),  // d F / d x_1
		2 * q.E_1 * q.E_2 * (-q.y_2 + q.z_2 * q.y_1 / q.z_1),  // d F / d y_1
		2 * q.E_2 * term,                                      // d F / d E_1
		2 * q.E_1 * q.E_2 * (-q.x_1 + q.z_1 * q.x_2 / q.z_2),  // d F / d x_2
		2 * q.E_1 * q.E_2 * (-q.y_1 + q.z_1 * q.y_2 / q.z_2),  // d F / d y_2
		2 * q.E_1 * term                                       // d F / d E_2
	};
	jacobianEta.ResizeTo(1, 6);
	jacobianEta.SetMatrixArray(J);
}
