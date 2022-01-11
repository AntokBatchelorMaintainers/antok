#ifndef ANTOK_GENERATORS_FUNCTIONS_H
#define ANTOK_GENERATORS_FUNCTIONS_H

#include <yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace generators {

		std::string mergeNameIndex(const std::string& name, const int index);
		bool nmbArgsIsExactly(const YAML::Node& function, const size_t& actualNmb, const size_t& requiredNmb);
		std::string getTypeOfArg(const YAML::Node&  function, const int index, const std::string& argName);
		bool functionArgumentHandler(std::vector<std::pair<std::string, std::string>>& args,
		                             const YAML::Node&                                 function,
		                             const int                                         index,
		                             const bool                                        argStringsAlreadyValues = false);
		/**
		 * Sets values in given args vector to the constants found in YAML node
		 * @param args Vector of pairs where first: node name, second: node value (will be copied in this function)
		 * @param function: Node of the function
		 * @param required: if set error message is printed if args do not exist
		 * @return true if everything was ok
		 */
		template <typename T>
		bool functionArgumentHandlerConst(std::map<std::string, T>& args,
		                                  const YAML::Node&         function,
		                                  const bool                required = true);
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
		#define FUNCTION_PROTOTYPE(fctName) antok::Function* fctName(const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index)
		FUNCTION_PROTOTYPE(generateAbs);
		FUNCTION_PROTOTYPE(generateLog);
		FUNCTION_PROTOTYPE(generateSqrt);
		FUNCTION_PROTOTYPE(generateConvertIntToDouble);
		FUNCTION_PROTOTYPE(generateDiff);
		FUNCTION_PROTOTYPE(generateQuotient);
		FUNCTION_PROTOTYPE(generateMul);
		FUNCTION_PROTOTYPE(generateEnergy);
		FUNCTION_PROTOTYPE(generateGetBeamLorentzVector);
		FUNCTION_PROTOTYPE(generateGetGradXGradY);
		FUNCTION_PROTOTYPE(generateGetLorentzVectorAttributes);
		FUNCTION_PROTOTYPE(generateGetLorentzVec);
		FUNCTION_PROTOTYPE(generateGetTs);
		FUNCTION_PROTOTYPE(generateGetVector2);
		FUNCTION_PROTOTYPE(generateGetVector3);
		FUNCTION_PROTOTYPE(generateGetVectorEntry);
		FUNCTION_PROTOTYPE(generateGetVectorSize);
		FUNCTION_PROTOTYPE(generateMass);
		FUNCTION_PROTOTYPE(generateMass2);
		FUNCTION_PROTOTYPE(generateRadToDegree);
		FUNCTION_PROTOTYPE(generateSum);
		FUNCTION_PROTOTYPE(generateSum2);

	}

}

#endif  // ANTOK_GENERATORS_FUNCTIONS_H
