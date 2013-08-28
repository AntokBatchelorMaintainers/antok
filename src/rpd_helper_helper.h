#ifndef ANTOK_RPD_HELPER_HELPER_H
#define ANTOK_RPD_HELPER_HELPER_H

#include<vector>

class TVector2;

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

		static RpdHelperHelper* getInstance();

		std::pair<double, double> getLimits(const double& rpdPhi,
		                                    const TVector2& vertexXY);

		std::pair<double, double> getSharpSigmas(const double& m,
		                                         const double& q);

	  private:

		enum {
			LEFTSLAB,
			MIDDLESLAB,
			RIGHTSLAB
		};

		RpdHelperHelper();

		double evaluatePlaneFunction(const double& x,
		                             const double& y,
		                             const FitParameters& fitParameters);

		double evaluateSharpSigmaFunction(const double& m,
		                                  const double& q,
		                                  const FitParameters& fitParameters);

		std::pair<double, TVector2> rotateRpdAndVertex(const double& rpdPhi,
		                                               const TVector2 vertexXY);

		static RpdHelperHelper* _instance;

		std::vector<double> RPD_SLAB_PHI_ANGLES;

		std::pair<antok::FitParameters, antok::FitParameters> _middleSlabLimitsParameters;
		std::pair<antok::FitParameters, antok::FitParameters> _leftSlabLimitsParameters;
		std::pair<antok::FitParameters, antok::FitParameters> _rightSlabLimitsParameters;

		std::pair<antok::FitParameters, antok::FitParameters> _sharpSigmaParameters;

	};

}

#endif
