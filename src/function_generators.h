#ifndef ANTOK_FUNCTION_GENERATORS_H
#define ANTOK_FUNCTION_GENERATORS_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace generators {

		antok::Function* registerAbs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* registerDiff(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* registerGetBeamLorentzVector(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* registerGetLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* registerGetRpdPhi(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* registerGetTs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* registerMass(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* registerSum(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* registerSum2(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);

		bool functionArgumentHandler(std::vector<std::pair<std::string, std::string> >& args,
		                             const YAML::Node& function,
		                             int index,
		                             bool argStringsAlreadyValues = false);

		std::string getFunctionArgumentHandlerErrorMsg(std::vector<std::string> quantityNames);

	}

}

#endif

