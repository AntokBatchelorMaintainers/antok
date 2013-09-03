#ifndef ANTOK_RPD_HELPER_HELPER_H
#define ANTOK_RPD_HELPER_HELPER_H

#include<vector>

#include<TVector2.h>

class TF1;

namespace antok {

	class FitParameters {

	  public:

		FitParameters() { }
		FitParameters(const std::vector<double>& parameters) { _parameters = parameters; }
		const double& getParameter(const unsigned int& i) const { return _parameters[i]; }
		void setParameter(const unsigned int& i, const double& value) { _parameters[i] = value; }

	  private:

		std::vector<double> _parameters;

	};

	class RpdHelperHelper {

	  public:

		static const RpdHelperHelper* getInstance();

		double getLikelihood(const double& deltaPhi,
		                     const double& rpdProtonPhi,
		                     const double& vertexX,
		                     const double& vertexY,
		                     const double& xMass,
		                     const double& predictedProtonMom) const;

		double getScalingFactor(const double& xMass,
		                        const double& predictedProtonMom) const;

		std::pair<double, double> getLimits(const double& rpdPhi,
		                                    const TVector2& vertexXY) const;

		std::pair<double, double> getSharpSigmas(const double& m,
		                                         const double& q) const;
		std::pair<double, double> getBroadSigmas(const double& m,
		                                         const double& q) const;

		double getBroadSigmasScalingFactor(const double& m,
		                                   const double& q) const;

		bool rpdProtonPhiValid(const double& rpdProtonPhi) const;

	  private:

		enum {
			LEFTSLAB,
			MIDDLESLAB,
			RIGHTSLAB
		};

		RpdHelperHelper();

		std::vector<double> getAllLimits() const;

		double getScalingFactor(const double& broadSigmasScalingFactor) const;

		double evaluatePlaneFunction(const double& x,
		                             const double& y,
		                             const FitParameters& fitParameters) const;

		double evaluateSharpSigmaFunction(const double& m,
		                                  const double& q,
		                                  const FitParameters& fitParameters) const;

		double evaluateBroadSigmaFunction(const double& m,
		                                  const double& q,
		                                  const FitParameters& fitParameters) const;

		std::pair<double, TVector2> rotateRpdAndVertex(const double& rpdPhi,
		                                               const TVector2 vertexXY) const;

		unsigned int getHitSlab(const double& correctedProtonPhi) const;

		static RpdHelperHelper* _instance;

		std::vector<double> RPD_SLAB_PHI_ANGLES;

		std::pair<antok::FitParameters, antok::FitParameters> _middleSlabLimitsParameters;
		std::pair<antok::FitParameters, antok::FitParameters> _leftSlabLimitsParameters;
		std::pair<antok::FitParameters, antok::FitParameters> _rightSlabLimitsParameters;

		std::pair<antok::FitParameters, antok::FitParameters> _sharpSigmaParameters;
		antok::FitParameters _broadSigmaScalingFactor;
		std::pair<double, double> _broadSigmaParameters;

		std::vector<double> _allowedProtonPhiValues;

		TF1* _likelihoodFunction;

		mutable double __correctedProtonPhi;
		mutable TVector2 __correctedVertexXY;

	};

}

#endif
