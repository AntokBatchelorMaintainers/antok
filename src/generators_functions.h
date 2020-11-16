#ifndef ANTOK_GENERATORS_FUNCTIONS_H
#define ANTOK_GENERATORS_FUNCTIONS_H

#include <yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace generators {

		std::string mergeNameIndex(const std::string& name, const int index);
		bool nmbArgsIsExactly(const YAML::Node& function, const size_t& actualNmb, const size_t& requiredNmb);
		//TODO introduce typedef for args (at least in .cxx file) to avoid the long type name
		bool functionArgumentHandler(std::vector<std::pair<std::string, std::string>>& args,
		                             const YAML::Node&                                 function,
		                             const int                                         index,
		                             const bool                                        argStringsAlreadyValues = false);
		/**
		 * Sets values in given args vector to the constants found in YAML node
		 * @param args Vector of pairs where first: node name, second: node value (will be copied in this function)
		 * @param function: Node of the function
		 * @return true if everything was ok
		 */
		template <typename T>
		bool functionArgumentHandlerConst(std::map<std::string, T>& args,
		                                  const YAML::Node&         function);
		/**
		* Sets the data pointers in the args vector to the address of the variable or to a constant if no variable name, but a number is given
		* @param args Vector of pairs where first: node/variable name, second: data pointer (will be set in this function)
		* @param function: Node of the function
		* @param index: Index of the function call (0 if this arguments have no index)
		* @return true if everything was ok
		*/
		template <typename T>
		bool functionArgumentHandlerPossibleConst(std::vector<std::pair<std::string, T*>>& args,
		                                          const YAML::Node&                        function,
		                                          const int                                index);
		std::string getFunctionArgumentHandlerErrorMsg(const std::vector<std::string>& quantityNames);

		// macro for shorter prototyping of the antok functions
		#define FUNCTIONPROTOTYPE(fctName) antok::Function* fctName (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index)
		FUNCTIONPROTOTYPE(generateAbs);
		FUNCTIONPROTOTYPE(generateLog);
		FUNCTIONPROTOTYPE(generateConvertIntToDouble);
		FUNCTIONPROTOTYPE(generateDiff);
		FUNCTIONPROTOTYPE(generateQuotient);
		FUNCTIONPROTOTYPE(generateMul);
		FUNCTIONPROTOTYPE(generateEnergy);
		FUNCTIONPROTOTYPE(generateGetBeamLorentzVector);
		FUNCTIONPROTOTYPE(generateGetGradXGradY);
		FUNCTIONPROTOTYPE(generateGetLorentzVectorAttributes);
		FUNCTIONPROTOTYPE(generateGetLorentzVec);
		FUNCTIONPROTOTYPE(generateGetTs);
		FUNCTIONPROTOTYPE(generateGetVector3);
		FUNCTIONPROTOTYPE(generateGetVectorEntry);
		FUNCTIONPROTOTYPE(generateMass);
		FUNCTIONPROTOTYPE(generateMass2);
		FUNCTIONPROTOTYPE(generateRadToDegree);
		FUNCTIONPROTOTYPE(generateSum);
		FUNCTIONPROTOTYPE(generateSum2);
	}

}

#endif  // ANTOK_GENERATORS_FUNCTIONS_H
