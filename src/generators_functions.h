#ifndef ANTOK_GENERATORS_FUNCTIONS_H
#define ANTOK_GENERATORS_FUNCTIONS_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace generators {

		antok::Function* generateAbs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateConvertIntToDouble(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateDiff(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateEnergy(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetBeamLorentzVector(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetGradXGradY(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetLorentzVectorAttributes(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetRpdPhi(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetTs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetVector3(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateMass(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateRadToDegree(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateSum(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateSum2(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);

	}

}

#endif

