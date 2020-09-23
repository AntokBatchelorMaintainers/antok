#ifndef NEUTRALPROBLEM_H
#define NEUTRALPROBLEM_H

#include "TVectorD.h"
#include "kinematic_fit.h"

namespace antok {
	class NeutralProblem : public KinematicFit::problem {
		public:
			NeutralProblem(double mass_);

			~NeutralProblem() {}

			const TVectorD& constraint(const TVectorD& eta);
			void dConstraint(const TVectorD& eta, TMatrixD& B);
			bool converged() { return c.NormInf() / mass2 < 1e-10; }

			unsigned int getNConstraints() { return 1; }

		private:
			const double mass2;

			TVectorD c;
	};
}

#endif
