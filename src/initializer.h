#ifndef INITIALIZER_H
#define INITIALIZER_H

#include<string>

namespace antok {

	class Event;
	class Cutter;
	class Plotter;

	class Initializer {

	  public:

		static Initializer* instance();

		bool readConfigFile(const std::string& filename);

		bool init();

		antok::Cutter& get_cutter();
		antok::Event& get_event();
		antok::Plotter& get_plotter();

	  private:

		Initializer();

		static Initializer* _initializer;

		antok::Cutter* _cutter;
		antok::Event* _event;
		antok::Plotter* _plotter;

	};

}

#endif

