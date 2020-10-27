#ifndef ANTOK_INITIALIZER_H
#define ANTOK_INITIALIZER_H

#include <string>

#include <yaml-cpp/yaml.h>
#include "generators_functions.h"

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

		bool initAll()
		{
			return (    initializeData()
			        and initializeInput()
			        and initializeEvent()
			        and initializeCutter()
			        and initializePlotter());
		}

		bool initializeCutter();
		bool initializeData();
		bool initializeInput();
		bool initializeEvent();
		bool initializePlotter();
		bool updateInput();

	private:

		Initializer();

		static Initializer* _initializer;
		YAML::Node*         _config;

		antok::Function* getFunction(const YAML::Node function, const std::vector<std::string> quantityNames, const int index);

	};

}

#endif  // ANTOK_INITIALIZER_H
