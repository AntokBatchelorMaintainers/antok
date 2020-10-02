#ifndef KINEMATICFIT_H
#define KINEMATICFIT_H

#include <cassert>

#include "TMath.h"
#include "TVectorD.h"

namespace antok {

	class KinematicFit {

	public:

		class problem {

		public:

			virtual ~problem() { }

			virtual const TVectorD& constraint(const TVectorD& xp, const TVectorD& etap) { assert(0); };
			virtual const TVectorD& constraint(const TVectorD& etap) { assert(0); }
			virtual void dConstraint(const TVectorD& xp, const TVectorD& etap, TMatrixD& A, TMatrixD& B) { assert(0); }
			virtual void dConstraint(const TVectorD& etap, TMatrixD& B) { assert(0); }
			virtual bool converged() = 0;
			virtual size_t getNConstraints() = 0;
		};

		KinematicFit(problem& prob, const TVectorD& x, const TVectorD& eta, const TMatrixDSym& coveta);
		KinematicFit(problem& prob, const TVectorD& eta, const TMatrixDSym& coveta);
		~KinematicFit() { }

		void step();
		bool doFit();
		size_t getNSteps() const { return _nSteps; }
		void setMaxSteps(const size_t nMax) { _maxSteps = nMax; }
		TVectorD getParameters() const { return _x + _dx; }
		TVectorD getEnhanced() const { return _eta + _deta; }
		TMatrixDSym getCovParams();
		TMatrixDSym getCovEnhanced();
		Double_t getChi2() const;
		Double_t getCL() const {
			return TMath::Prob(getChi2(), _myProblem.getNConstraints());
		}

	private:

		int         _nParams;
		size_t      _nSteps;
		size_t      _maxSteps;
		TVectorD    _x;
		TVectorD    _eta;
		TMatrixDSym _coveta;
		TVectorD    _dx;
		TVectorD    _deta;
		problem&    _myProblem;

	};

}

#endif  // KINEMATICFIT_H
