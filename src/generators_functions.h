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
		/**
		* Sets the data pointers in the args vector to the address of the variable or to an constant if no variable name, but a number is given
		* @param args Vector of pairs where first: node/variable name, second: data pointer (will be set in this function)
		* @param function: Node of the function
		* @param index: Index of the function call (0 if this arguments have no index)
		* @return true if everything was ok
		*/
		template<typename T>
		bool functionrgumentHandlerPossibleConst(std::vector<std::pair<std::string, T*> >& args,
		                                         const YAML::Node& function,
		                                         int index);
		std::string getFunctionArgumentHandlerErrorMsg(std::vector<std::string> quantityNames);

		antok::Function* generateAbs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		antok::Function* generateLog(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
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

