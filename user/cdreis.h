#ifndef ANTOK_USER_CDREIS_H
#define ANTOK_USER_CDREIS_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace user {

		namespace cdreis {

			antok::Function* getUserFunction(const YAML::Node& function,
			                                 std::vector<std::string>& quantityNames,
			                                 int index);

			antok::Function* generateGetRecoilLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetPhotonLorentzVecs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetVectorLorentzVectorAttributes(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetPi0s(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetPi0Pair(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetCleanedEcalClusters(const YAML::Node &function, std::vector<std::string> &quantityNames, int index);
			antok::Function* generateGetOmega(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetECALCorrectedTiming(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetECALCorrectedEnergy(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);

			const double MASS_PI0   = 0.1349766;
			const double MASS_OMEGA = 0.7826500;

			const double RESOLUTION_ECAL1        = 3 * 0.00884872;
			const double RESOLUTION_ECAL2        = 3 * 0.00408290;
			const double RESOLUTION_ECALCOMBINED = 3 * 0.00831821;

			const double RESOLUTION_OMEGA        = 3 * 0.0108543;

			const double THRESHOLD_ENERGY_ECAL1 = 0.6;
			const double THRESHOLD_ENERGY_ECAL2 = 1.2;

			const double THRESHOLD_TIMING_ECAL1 = 3.75;
			const double THRESHOLD_TIMING_ECAL2 = 3.00;

			const int POSITION_ECAL1 = 2500;

		}

	}

}

#endif
