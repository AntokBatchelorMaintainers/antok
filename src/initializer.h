#ifndef ANTOK_INITIALIZER_H
#define ANTOK_INITIALIZER_H

#include<string>

#include<yaml-cpp/yaml.h>

class TFile;
class TTree;

namespace antok {

	class Event;
	class Cutter;
	class ObjectManager;
	class Plotter;

	class Initializer {

	  public:

		static Initializer* instance();

		bool readConfigFile(const std::string& filename);

		bool initAll() {
			return (
				initializeData() and
				initializeEvent() and
				initializeCutter() and
				initializePlotter()
			);
		}

		bool initializeCutter();
		bool initializeData();
		bool initializeEvent();
		bool initializePlotter();

	  private:

		bool registerAbs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		bool registerDiff(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		bool registerGetBeamLorentzVector(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		bool registerGetLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		bool registerGetRpdPhi(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		bool registerGetTs(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		bool registerMass(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		bool registerSum(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);
		bool registerSum2(const YAML::Node& function, std::vector<std::string>& quantityNames, int index);

		bool argumentHandler(std::vector<std::pair<std::string, std::string> >& args,
		                     const YAML::Node& function,
		                     int index,
		                     bool argStringsAlreadyValues = false);

		Initializer();

		static Initializer* _initializer;

		YAML::Node* _config;

		std::string getYAMLStringSafe(const YAML::Node& node);

		std::string getArgumentHandlerErrorMsg(std::vector<std::string> quantityNames);
		std::string getVariableInsertionErrorMsg(std::vector<std::string> quantityNames, std::string variableName = "");
		std::string getVariableInsertionErrorMsg(std::string variableName);

	};

}

#endif

