#ifndef ANTOK_USER_KBICKER_H
#define ANTOK_USER_KBICKER_H

#include<yaml-cpp/yaml.h>

class TLorentzVector;
class TVector3;

namespace antok {

	class Function;

	namespace user {

		namespace kbicker {

			antok::Function* getUserFunction(const YAML::Node& function,
			                                 const std::vector<std::string>& quantityNames,
			                                 int index);

			antok::Function* generateGetRpdExpectedHitsParameters(const YAML::Node& function, const std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetRpdPhi(const YAML::Node& function, const std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetCutOnExtraTracks(const YAML::Node& function, const std::vector<std::string>& quantityNames, int index);


			void getRPDDeltaPhiResProjection(const TLorentzVector& pBeam,
			                                 const TLorentzVector& pProton,
			                                 const TLorentzVector& pX,
			                                 double& delta_phi, double& res);

			void getRPDDeltaPhiResRotation(const TLorentzVector& pBeam,
			                               const TLorentzVector& pProton,
			                               const TLorentzVector& pX,
			                               double& delta_phi, double& res,
			                               double& phiProton, double& phiX);

			void getRPDDeltaPhiResRotation(const TLorentzVector& pBeam,
			                               const TLorentzVector& pProton,
			                               const TLorentzVector& pX,
			                               double& delta_phi, double& res);

			void getRPDDeltaPhiResPrediction(const TLorentzVector& pBeam,
			                                 const TLorentzVector& pProton,
			                                 const TLorentzVector& pX,
			                                 const TVector3& vertex,
			                                 double& delta_phi, double& likelihood,
			                                 double& phiProton, double& phiX);

			void getRPDDeltaPhiResPrediction(const TLorentzVector& pBeam,
			                                 const TLorentzVector& pProton,
			                                 const TLorentzVector& pX,
			                                 const TVector3& vertex,
			                                 double& delta_phi, double& likelihood);

			void getRPDExpectedHitsParameters(const TLorentzVector& pBeam,
			                                  const TLorentzVector& pX,
			                                  const TVector3& vertex,
			                                  const double& xOffset,
			                                  const double& yOffset,
			                                  const double& xAngle,
			                                  const double& yAngle,
			                                  double& rpdPhiRingA,
			                                  double& rpdPhiRingB,
			                                  double& rpdZRingA,
			                                  double& rpdZRingB);

			bool extraTracksCut(std::vector<double> trackTimes,
			                    std::vector<double> trackTimeSigmas,
			                    std::vector<double> trackNHits,
			                    std::vector<double> trackZFirst,
			                    std::vector<double> trackZLast,
			                    std::vector<double> trackQP);

		}

	}

}

#endif
