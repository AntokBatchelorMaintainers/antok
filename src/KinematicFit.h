#ifndef KINEMATICFIT_H
#define KINEMATICFIT_H

#include <assert.h>

#include "TVectorD.h"
#include "TMath.h"

namespace antok {
	class KinematicFit {
		public:
			class problem {
			public:
				virtual const TVectorD &constraint(const TVectorD &xp, const TVectorD &etap) { assert(0); };

				virtual const TVectorD &constraint(const TVectorD &etap) { assert(0); }

				virtual void dConstraint(const TVectorD &xp, const TVectorD &etap, TMatrixD &A, TMatrixD &B) { assert(0); }

				virtual void dConstraint(const TVectorD &etap, TMatrixD &B) { assert(0); }

				virtual bool converged() = 0;

				virtual unsigned int getNConstraints() = 0;

				virtual ~problem() {};
			};

			KinematicFit(problem &prob, const TVectorD &x_, const TVectorD &eta_, const TMatrixDSym &coveta_);

			KinematicFit(problem &prob, const TVectorD &eta_, const TMatrixDSym &coveta_);

			void step();

			bool doFit();

			unsigned int getNSteps() { return nSteps; };

			void setMaxSteps(unsigned int nMax) { maxSteps = nMax; }

			TVectorD getParameters() { return x + dx; }

			TVectorD getEnhanced() { return eta + deta; }

			TMatrixDSym getCovParams();

			TMatrixDSym getCovEnhanced();

			Double_t getChi2();

			Double_t getCL() {
				return TMath::Prob(getChi2(), myProblem.getNConstraints());
			}

		private:
			int numParams;

			unsigned int nSteps;
			unsigned int maxSteps;
			TVectorD x;
			TVectorD eta;
			TMatrixDSym coveta;
			TVectorD dx, deta;

			problem &myProblem;
	};
}

#endif
