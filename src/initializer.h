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

		bool registerGetBeamLorentzVector(const YAML::Node& function, std::string& quantityName, int index);
		bool registerGetLorentzVec(const YAML::Node& function, std::string& quantityName, int index);
		bool registerGetTs(const YAML::Node& function, std::vector<std::string>& quantityName, int index);
		bool registerMass(const YAML::Node& function, std::string& quantityName, int index);
		bool registerSum(const YAML::Node& function, std::string& quantityName, int index);
		bool registerSum2(const YAML::Node& function, std::string& quantityName, int index);

		bool argumentHandler(std::vector<std::string>& argNames, int index);

		Initializer();

		static Initializer* _initializer;

		YAML::Node* _config;

	};

}

#endif

