#ifndef ANTOK_USER_CDREIS_H
#define ANTOK_USER_CDREIS_H

#include <yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace user {

		namespace cdreis {

			antok::Function* getUserFunction                         (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);

			antok::Function* generateGetVector3VectorAttributes      (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetVectorLorentzVectorAttributes(const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetNominalMassDifferences       (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetRecoilLorentzVec             (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetECALCorrectedEnergy          (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetECALCorrectedTiming          (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetCleanedEcalClusters          (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetECALVariables                (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetPhotonLorentzVecs            (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetPhotonPairParticles          (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetPi0Pair                      (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetKinematicFittingMass         (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetOmega                        (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
			antok::Function* generateGetThreePionCombinationMass     (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		}

	}

}

#endif  // ANTOK_USER_CDREIS_H
