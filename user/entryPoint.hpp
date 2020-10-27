#ifndef ANTOK_USER_ENTRYPOINT_H
#define ANTOK_USER_ENTRYPOINT_H

#include "functions.hpp"

// Users should only change things within the following block!
//
// To add your user functions, add the #include "username.h" line for your file
// and push back your getUserFunction() into the userFunctions std::vector<>
// below.
#include "kbicker.h"
#include "hubers.h"
#include "cdreis.h"
#include "swallner.h"

namespace {

	std::vector<antok::Function* (*)(const YAML::Node&, const std::vector<std::string>&, int)>
	__getUserFunctions(const YAML::Node&               function,
	                   const std::vector<std::string>& quantityNames,
	                   int                             index)
	{
		std::vector<antok::Function* (*)(const YAML::Node&, const std::vector<std::string>&, int)> userFunctions;
		// User defined function have to be added here

		userFunctions.push_back(&antok::user::kbicker::getUserFunction);
		userFunctions.push_back(&antok::user::hubers::getUserFunction);
		userFunctions.push_back(&antok::user::cdreis::getUserFunction);
		userFunctions.push_back(&antok::user::stefan::getUserFunction);

		// End of user defined functions
		return userFunctions;
	}

}
// End of user changable block, no changes by users after this line!

namespace antok {

	namespace user {

		antok::Function*
		getUserFunction(const YAML::Node&               function,
		                const std::vector<std::string>& quantityNames,
		                int                             index)
		{
			std::vector<antok::Function* (*)(const YAML::Node&, const std::vector<std::string>&, int)>
				userFunctions = __getUserFunctions(function, quantityNames, index);
			antok::Function* retval = nullptr;
			for (size_t i = 0; i < userFunctions.size(); ++i) {
				antok::Function* userFunction = (*userFunctions[i])(function, quantityNames, index);
				if (userFunction == nullptr) {
					continue;
				}
				if (retval != nullptr) {
					std::cerr << "More than one user defined function responded to function name '"
					          << function["Name"] << "' (duplicate names?). Aborting..." << std::endl;
					return nullptr;
				}
				retval = userFunction;
			}
			return retval;
		}

	}

}

#endif  // ANTOK_USER_ENTRYPOINT_H
