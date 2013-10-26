#ifndef ANTOK_USER_HUBERS_H
#define ANTOK_USER_HUBERS_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Function;

	namespace user {

		namespace hubers {

			antok::Function* getUserFunction(const YAML::Node& function,
			                                 std::vector<std::string>& quantityNames,
			                                 int index);

			antok::Function* generateGetNeuronalBeamEnergy(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);

			void getNeuronalBeamEnergy(const double& X,
																 const double& Y, 
																 const double& dX, 
																 const double& dY, 
																 double& E); 
 

		}

	}

}

#endif
