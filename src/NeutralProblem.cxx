#include "NeutralProblem.h"

antok::NeutralProblem::NeutralProblem(double mass_)
		: mass2(mass_*mass_) {
	c.ResizeTo(1);
}


const TVectorD& antok::NeutralProblem::constraint(const TVectorD& eta) {
	const double x1 = eta[0];
	const double y1 = eta[1];
	const double z1 = sqrt(1 - x1*x1 - y1*y1);
	const double E1 = eta[2];

	const double x2 = eta[3];
	const double y2 = eta[4];
	const double z2 = sqrt(1 - x2*x2 - y2*y2);
	const double E2 = eta[5];

	c[0] = 2*E1*E2*(1 - x1*x2 - y1*y2 - z1*z2) - mass2;
	return c;
}


void antok::NeutralProblem::dConstraint(const TVectorD& eta,TMatrixD& B) {
	const double x1 = eta[0];
	const double y1 = eta[1];
	const double z1 = sqrt(1 - x1*x1 - y1*y1);
	const double E1 = eta[2];

	const double x2 = eta[3];
	const double y2 = eta[4];
	const double z2 = sqrt(1 - x2*x2 - y2*y2);
	const double E2 = eta[5];

	const double b[] = {
			2*E1*E2*(-x2 + z2*x1/z1),
			2*E1*E2*(-y2 + z2*y1/z1),
			2*E2*(1 - x1*x2 - y1*y2 - z1*z2),

			2*E1*E2*(-x1 + z1*x2/z2),
			2*E1*E2*(-y1 + z1*y2/z2),
			2*E1*(1 - x1*x2 - y1*y2 - z1*z2),
	};

	B.ResizeTo(1, 6);
	B.SetMatrixArray(b);
}
