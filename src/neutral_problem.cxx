#include "neutral_problem.h"


antok::NeutralProblem::NeutralProblem(const double mass)
	: _mass2(mass * mass)
{
	_c.ResizeTo(1);
}


const TVectorD&
antok::NeutralProblem::constraint(const TVectorD& eta)
{
	const double x1 = eta[0];
	const double y1 = eta[1];
	const double z1 = sqrt(1 - x1 * x1 - y1 * y1);
	const double E1 = eta[2];

	const double x2 = eta[3];
	const double y2 = eta[4];
	const double z2 = sqrt(1 - x2 * x2 - y2 * y2);
	const double E2 = eta[5];

	_c[0] = 2 * E1 * E2 * (1 - x1 * x2 - y1 * y2 - z1 * z2) - _mass2;
	return _c;
}


void
antok::NeutralProblem::dConstraint(const TVectorD& eta,
                                   TMatrixD&       B)
{
	//TODO reduce code repetition: same as in antok::NeutralProblem::constraint()
	const double x1 = eta[0];
	const double y1 = eta[1];
	const double z1 = sqrt(1 - x1 * x1 - y1 * y1);
	const double E1 = eta[2];

	const double x2 = eta[3];
	const double y2 = eta[4];
	const double z2 = sqrt(1 - x2 * x2 - y2 * y2);
	const double E2 = eta[5];

	const double b[] = {
		2 * E1 * E2 * (-x2 + z2 * x1 / z1),
		2 * E1 * E2 * (-y2 + z2 * y1 / z1),
		2 * E2 * (1 - x1 * x2 - y1 * y2 - z1 * z2),

		2 * E1 * E2 * (-x1 + z1 * x2 / z2),
		2 * E1 * E2 * (-y1 + z1 * y2 / z2),
		2 * E1 * (1 - x1 * x2 - y1 * y2 - z1 * z2),
	};

	B.ResizeTo(1, 6);
	B.SetMatrixArray(b);
}
