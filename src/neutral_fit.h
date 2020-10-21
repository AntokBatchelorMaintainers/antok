// Based on code by Tobias Schlüter  <tobias.schlueter@physik.uni-muenchen.de>,
// available at
// https://gitlab.cern.ch/compass/hadron/hadrontools/-/tree/master/kinematicFit
// and
// http://wwwcompass.cern.ch/compass/notes/2007-10/source.tar.gz
// and documented in COMPASS note 2007-10 available at
// https://wwwcompass.cern.ch/compass/notes/2007-10/2007-10.pdf


#ifndef NEUTRALFIT_H
#define NEUTRALFIT_H

#include <vector>

#include "TMatrixDSym.h"

#include "neutral_problem.h"

class TH1;

namespace antok {

	/*! This class collects everything that is needed to fit two
	    calorimeter clusters from photons to a given mass constraint */
	class NeutralFit {

	public:

		/*! \param[in] vertexPosition           position of the vertex from which the two photons originate
		    \param[in] cluster1Position         position of the first calorimeter cluster
		    \param[in] cluster2Position         position of the second calorimeter cluster
		    \param[in] cluster1PositionVariance variances of the position of the first calorimeter cluster
		    \param[in] cluster2PositionVariance variances of the position of the second calorimeter cluster
		    \param[in] cluster1Energy           energy of the first calorimeter cluster
		    \param[in] cluster2Energy           energy of the second calorimeter cluster
		    \param[in] cluster1EnergyVariance   variance of the energy of the first calorimeter cluster
		    \param[in] cluster2EnergyVariance   variance of the energy of the second calorimeter cluster
		    \param[in] mass                     mass the photon pair is constrained to
		    \param[in] massLowerLimit           lower mass limit used in massIsInWindow()
		    \param[in] massUpperLimit           upper mass limit used in massIsInWindow()
		    \param[in] whichEnergyVariance      defines how variance of cluster energy is calculated:
		                                        0 = variance as returned by PHAST
		                                        1 = square of cluster energy
		                                        2 = (ECAL2 resolution)^2 as given in COMPASS spectrometer paper [NIMA 577 (2007) 455] */
		NeutralFit(const TVector3& vertexPosition,
		           const TVector3& cluster1Position,
		           const TVector3& cluster2Position,
		           const TVector3& cluster1PositionVariance,
		           const TVector3& cluster2PositionVariance,
		           const double    cluster1Energy,
		           const double    cluster2Energy,
		           const double    cluster1EnergyVariance,
		           const double    cluster2EnergyVariance,
		           const double    mass,
		           const double    massLowerLimit,
		           const double    massUpperLimit,
		           const double    convergenceLimit = 1e-10,
		           const int       whichEnergyVariance = 0);
		~NeutralFit() { delete _kinFitter; }

		//! Perform the fit
		bool doFit();

		//! Return whether the measured mass of the two-photon system is within the defined mass window
		bool   massIsInWindow() const;
		//! Return chi^2 value of fit
		double chi2Value()      const { return _kinFitter->chi2Value(); }
		//! Return P-value of fit
		double pValue()         const { return _kinFitter->pValue();    }

		//! Return improved Lorentz vector of first photon
		const TLorentzVector& getImprovedLV1()   const { return _improvedLV_1;  }
		//! Return improved Lorentz vector of second photon
		const TLorentzVector& getImprovedLV2()   const { return _improvedLV_2;  }
		//! Return sum of the improved Lorentz vectors of the photon pair,
		//! i.e. this Lorentz vector will fulfill the mass constraint
		//! getImprovedLVSum().M() == _mass
		const TLorentzVector& getImprovedLVSum() const { return _improvedLVSum; }

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
		TMatrixDSym covForCluster(const TVector3& clusterPosition,
		                          const TVector3& clusterPositionVariance,
		                          const double    clusterEnergy,
		                          const double    clusterEnergyVariance);

		antok::NeutralProblem _problem;    //!< the fit problem
		antok::KinematicFit*  _kinFitter;  //!< the kinematic fitter

		static const size_t      _nmbFittedVars = 6;  //!< dimension of vector of measured variables that enter the fit
		static bool              _first;              //!< flag that steers initialization of static variables
		//TODO why do the histograms need to be static?
		static std::vector<TH1*> _pullsHists;         //!< vector with histograms for pulls
		std::vector<double>      _pullValues;         //!< pull values after fit

		const TVector3& _vertexPosition;            //!< position of the vertex from which the two photons originate
		const TVector3& _cluster1Position;          //!< position of the first calorimeter cluster
		const TVector3& _cluster2Position;          //!< position of the second calorimeter cluster
		const TVector3& _cluster1PositionVariance;  //!< variances of the position of the first calorimeter cluster
		const TVector3& _cluster2PositionVariance;  //!< variances of the position of the second calorimeter cluster
		const double    _cluster1Energy;            //!< energy of the first calorimeter cluster
		const double    _cluster2Energy;            //!< energy of the second calorimeter cluster
		const double    _cluster1EnergyVariance;    //!< variance of the energy of the first calorimeter cluster
		const double    _cluster2EnergyVariance;    //!< variance of the energy of the second calorimeter cluster

		TLorentzVector _improvedLV_1;   //!< Lorentz vector of first photon after fit
		TLorentzVector _improvedLV_2;   //!< Lorentz vector of second photon after fit
		TLorentzVector _improvedLVSum;  //!< sum of the Lorentz vectors of the photon pair after fit

		const double _mass;                 //!< mass the photon pair is constrained to
		const double _massLowerLimit;       //!< lower mass limit of mass window used in massIsInWindow()
		const double _massUpperLimit;       //!< upper mass limit of mass window used in massIsInWindow()
		const double _convergenceLimit;     //!< limit for determining if fit converged
		const int    _whichEnergyVariance;  //!< defines how variance of cluster energy is calculated:

		TVectorD    _startValues;  //!< start values for measured variables that enter the fit
		TMatrixDSym _covMat;       //!< covariance matrix of measured values

	};

}  // antok namespace

#endif  // NEUTRALFIT_H
