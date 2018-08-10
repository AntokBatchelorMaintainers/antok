#ifndef ANTOK_USER_stefan_H
#define ANTOK_USER_stefan_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace user {

		namespace stefan {

			antok::Function* getUserFunction(const YAML::Node& function,
			                                 std::vector<std::string>& quantityNames,
			                                 int index);
		antok::Function* getCalcLTProjections(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* getCalcArmenterosAlpha(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* getCalcRICHPID(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* getDetermineKaonPionLV(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* getDetermineKaonPionLVLikelihood(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* getCalcCEDARPID(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* getCalcCEDARPIDMulitL(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* getCalcCEDARPIDOneL(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* getCalcAngles3P(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* getCalcRapidityXF(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);

		}

	}

}

#endif
