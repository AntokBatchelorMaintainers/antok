#ifndef NEUTRALPROBLEM_H
#define NEUTRALPROBLEM_H

#include "TVectorD.h"

#include "kinematic_fit.h"

namespace antok {

	class NeutralProblem : public KinematicFit::problem {

	public:

		NeutralProblem(const double mass);
		~NeutralProblem() {}

		const TVectorD& constraint(const TVectorD& eta);
		void dConstraint(const TVectorD& eta, TMatrixD& B);
		bool converged() { return _c.NormInf() / _mass2 < 1e-10; }

		size_t getNConstraints() { return 1; }

	private:

		const double _mass2;
		TVectorD     _c;

	};

}

#endif  // NEUTRALPROBLEM_H
