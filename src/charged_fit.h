// Based on code by Tobias Schlüter  <tobias.schlueter@physik.uni-muenchen.de>,
// available at
// https://gitlab.cern.ch/compass/hadron/hadrontools/-/tree/master/kinematicFit
// and
// http://wwwcompass.cern.ch/compass/notes/2007-10/source.tar.gz
// and documented in COMPASS note 2007-10 available at
// https://wwwcompass.cern.ch/compass/notes/2007-10/2007-10.pdf


#ifndef CHARGEDFIT_H
#define CHARGEDFIT_H

#include <vector>
#include <iostream>

#include "TMatrixDSym.h"

#include "charged_problem.h"
#include "kinematic_fit.h"

class TH1;

namespace antok {

	/*! This class collects everything that is needed to fit two
	    secondary particles to a given invariant mass constraint */

	class ChargedFit {

	public:

		ChargedFit(const TVector3& particle1Momentum,
                   const TVector3& particle2Momentum,
                   const double&   particle1Energy,
                   const double&   particle2Energy,
                   const double    particle1Mass,
                   const double    particle2Mass,
                   const std::vector<double>& particle1MomentumCovariance,
                   const std::vector<double>& particle2MomentumCovariance,
                   const double    mass,
                   const double    massLowerLimit,
                   const double    massUpperLimit,
                   const double    precisionGoal);
		~ChargedFit() { delete _kinFitter; }

		//! Perform the fit
		bool doFit();

		//! Return whether the measured mass of the two-particle system is within the defined mass window
		bool   massIsInWindow() const;
		//! Return chi^2 value of fit
		double chi2Value()      const { return _kinFitter->chi2Value(); }
		//! Return P-value of fit
		double pValue()         const { return _kinFitter->pValue();    }
		//! Return number of performed iterations of fit
		int nmbIterations()     const { return _kinFitter->nmbSteps();  }

		//! Return improved Lorentz vector of first particle
		const TLorentzVector& getImprovedLV1()   const { return _improvedLV_1;  }
		//! Return improved Lorentz vector of second particle
		const TLorentzVector& getImprovedLV2()   const { return _improvedLV_2;  }
		//! Return sum of the improved Lorentz vectors of the particle pair,
		//! i.e. this Lorentz vector will fulfill the mass constraint
		//! getImprovedLVSum().M() == _mass
		const TLorentzVector& getImprovedLVSum() const { return _improvedLVSum; }
		const std::vector<double>& getTransfCov() const { return _transfCov; }

		//! Return vector of pull values
		const std::vector<double>& pullValues() const { return _pullValues; }
		//! Return pull histogram defined by index
		TH1* pullHist(const size_t index) { return _pullsHists[index]; }

	private:

		//! Initialize pull histograms
		void initPullHists();
		//! Fill pull histograms
		void fillPullHists(const TVectorD& improvedMeas);
		//! Calculate covariance matrix for fit variables
		TMatrixDSym covForParticle(const TVector3& particleMomentum,
	                               const double&   particleEnergy,
	                               const std::vector<double>& particleMomentumCovariance);

		antok::ChargedProblem _problem;    //!< the fit problem
		antok::KinematicFit*  _kinFitter;  //!< the kinematic fitter

		static const size_t      _nmbFittedVars = 6;  //!< dimension of vector of measured variables that enter the fit
		static bool              _first;              //!< flag that steers initialization of static variables
		static std::vector<TH1*> _pullsHists;         //!< vector with histograms for pulls
		std::vector<double>      _pullValues;         //!< pull values after fit

		const TVector3& _particle1Momentum;
		const TVector3& _particle2Momentum;
		const double&   _particle1Energy;
		const double&   _particle2Energy;
		const double    _particle1Mass;
		const double    _particle2Mass;
		const std::vector<double>& _particle1MomentumCovariance;
		const std::vector<double>& _particle2MomentumCovariance;

		TLorentzVector _improvedLV_1;   //!< Lorentz vector of first photon after fit
		TLorentzVector _improvedLV_2;   //!< Lorentz vector of second photon after fit
		TLorentzVector _improvedLVSum;  //!< sum of the Lorentz vectors of the photon pair after fit

		std::vector<double> _transfCov;

		const double _mass;                 //!< mass the photon pair is constrained to
		const double _massLowerLimit;       //!< lower mass limit of mass window used in massIsInWindow()
		const double _massUpperLimit;       //!< upper mass limit of mass window used in massIsInWindow()
		const double _precisionGoal;        //!< limit for determining if fit converged

		TVectorD    _startValues;  //!< start values for measured variables that enter the fit
		TMatrixDSym _covMat;       //!< covariance matrix of measured values

	};

}  // antok namespace

#endif  // CHARGEDFIT_H
