#ifndef ANTOK_USER_ENTRYPOINT_H
#define ANTOK_USER_ENTRYPOINT_H

#include<functions.hpp>

// Users should only change things within this block!
//
// To add your user functions, add the #include<username.h> line to your file
// and push back your getUserFunction() into the userFunctions std::vector<>
// below.
#include<kbicker.h>
namespace {
	std::vector<antok::Function* (*)(const YAML::Node&, std::vector<std::string>&, int)> __getUserFunctions(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {
		std::vector<antok::Function* (*)(const YAML::Node&, std::vector<std::string>&, int)> userFunctions;
		// User defined function have to be added here

		userFunctions.push_back(&antok::user::kbicker::getUserFunction);

		// End of user defined functions
		return userFunctions;
	}
}
// End of user changable block, no changes by users after this line!

namespace antok {

	namespace user {

		antok::Function* getUserFunction(const YAML::Node& function, std::vector<std::string>& quantityNames, int index) {

			std::vector<antok::Function* (*)(const YAML::Node&, std::vector<std::string>&, int)> userFunctions;
			userFunctions = __getUserFunctions(function, quantityNames, index);

			antok::Function* retval = 0;
			bool alreadyThere = false;
			for(unsigned int i = 0; i < userFunctions.size(); ++i) {
				if(alreadyThere) {
					std::cerr<<"More than one user defined function responded (dublicate names?). Aborting..."<<std::endl;
					return 0;
				}
				retval = (*userFunctions[i])(function, quantityNames, index);
				if(retval) {
					alreadyThere = true;
				}
			}
			return retval;

		}

	}

}

#endif
