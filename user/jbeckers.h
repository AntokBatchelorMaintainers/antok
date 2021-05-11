#ifndef ANTOK_USER_jbeckers_H
#define ANTOK_USER_jbeckers_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace user {

		namespace jbeckers {

			antok::Function* getUserFunction(const YAML::Node& function, const std::vector<std::string>& quantityNames, int index);

		    antok::Function* generateScale(const YAML::Node& function, const std::vector<std::string>& quantityNames, int index);

		}

	}

}

#endif
