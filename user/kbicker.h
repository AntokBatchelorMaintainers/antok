#ifndef ANTOK_USER_KBICKER_H
#define ANTOK_USER_KBICKER_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace user {

		namespace kbicker {

			antok::Function* getUserFunction(const YAML::Node& function,
			                                 std::vector<std::string>& quantityNames,
			                                 int index);

			antok::Function* generateGetRpdExpectedHitsParameters(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetRpdPhi(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);

		}

	}

}

#endif
