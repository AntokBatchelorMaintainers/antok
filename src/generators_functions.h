#ifndef ANTOK_GENERATORS_FUNCTIONS_H
#define ANTOK_GENERATORS_FUNCTIONS_H

#include <yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace generators {

		std::string mergeNameIndex(const std::string& name, const int index);
		//TODO introduce typedef for args (at least in .cxx file)
		bool functionArgumentHandler(std::vector<std::pair<std::string, std::string> >& args,
		                             const YAML::Node&                                  function,
		                             const int                                          index,
		                             const bool                                         argStringsAlreadyValues = false);
		/**
		* Sets the data pointers in the args vector to the address of the variable or to a constant if no variable name, but a number is given
		* @param args Vector of pairs where first: node/variable name, second: data pointer (will be set in this function)
		* @param function: Node of the function
		* @param index: Index of the function call (0 if this arguments have no index)
		* @return true if everything was ok
		*/
		template <typename T>
		bool functionArgumentHandlerPossibleConst(std::vector<std::pair<std::string, T*> >& args,
		                                          const YAML::Node&                         function,
		                                          const int                                 index);
		std::string getFunctionArgumentHandlerErrorMsg(const std::vector<std::string>& quantityNames);

		//TODO since all functions have the same signature, it would maybe easier to use a macro to generate them
		antok::Function* generateAbs                       (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateLog                       (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateConvertIntToDouble        (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateDiff                      (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateQuotient                  (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateMul                       (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateEnergy                    (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateGetBeamLorentzVector      (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateGetGradXGradY             (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateGetLorentzVectorAttributes(const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateGetLorentzVec             (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateGetTs                     (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateGetVector3                (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateGetVectorEntry            (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateMass                      (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateRadToDegree               (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateSum                       (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);
		antok::Function* generateSum2                      (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index);

	}

}

#endif  // ANTOK_GENERATORS_FUNCTIONS_H
