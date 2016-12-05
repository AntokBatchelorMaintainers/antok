#ifndef ANTOK_GENERATORS_FUNCTIONS_H
#define ANTOK_GENERATORS_FUNCTIONS_H

#include<yaml-cpp/yaml.h>


namespace antok {

	class Function;

	namespace generators {

		std::string mergeNameIndex( std::string const& name, int const index );
		bool functionArgumentHandler(std::vector<std::pair<std::string, std::string> >& args,
		                             const YAML::Node& function,
		                             int index,
		                             bool argStringsAlreadyValues = false);
		std::string getFunctionArgumentHandlerErrorMsg(std::vector<std::string> quantityNames);

		antok::Function* generateAbs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateConvertIntToDouble(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateDiff(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateQuotient(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateMul(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateEnergy(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetBeamLorentzVector(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetGradXGradY(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetLorentzVectorAttributes(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetTs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateGetVector3(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateMass(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateRadToDegree(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateSum(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateSum2(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);

	}

}

#endif

