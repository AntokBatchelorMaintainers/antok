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

			antok::Function* generateSqrt(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateFrac(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateThetaRICH(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetPt(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateEnforceEConservation(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetNeuronalBeam(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetTheta(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetThetaZCut(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetBadSpill(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetShifted(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetScaledCluster(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetCleanedClusters(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetMaximumCluster(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetNeutralLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetFormFactor(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetBgTrackCut(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetClosestPi0(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetClosestPi0Pi0(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateRunSpill(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateGetClusterPosCor(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
			antok::Function* generateExtrapNeutral(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);

			void getNeuronalBeamEnergy(const double& X, const double& Y, const double& dX, const double& dY, double& E);

			double PEDepGammaCorrection(double Egamma, double cx, double cy);
			double LinearGammaCorrection(double Egamma);
			double IntraCellX(double cx);
			double IntraCellY(double cy);

		}

	}

}

#endif
