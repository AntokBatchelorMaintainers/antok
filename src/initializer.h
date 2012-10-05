#ifndef ANTOK_INITIALIZER_H
#define ANTOK_INITIALIZER_H

#include<string>

#include<yaml-cpp/yaml.h>

class TFile;
class TTree;

namespace antok {

	class Data;
	class Event;
	class Cutter;
	class Plotter;

	class Initializer {

	  public:

		static Initializer* instance();

		bool readConfigFile(const std::string& filename);

		bool init();

		antok::Cutter& get_cutter();
		antok::Data& get_data(TFile* infile, TTree*& intree);
		antok::Event& get_event();
		antok::Plotter& get_plotter();

	  private:

		Initializer();

		static Initializer* _initializer;

		YAML::Node* _config;

		antok::Cutter* _cutter;
		antok::Data* _data;
		antok::Event* _event;
		antok::Plotter* _plotter;

	};

}

#endif

