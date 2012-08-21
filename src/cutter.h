#ifndef HLIB_CUTTER_H
#define HLIB_CUTTER_H

#include<string>
#include<vector>

#include<cut.hpp>
#include<event.h>

namespace hlib {

	class Cutter {
	
	  public:

		static Cutter* instance();

		int get_cutmask(const hlib::Event& event);

		std::string get_abbreviations(int bitmask);

		unsigned int get_no_cuts() { return _cuts.size(); };

	  private:

		Cutter();
		~Cutter();

		static Cutter* _cutter;

		std::vector<hlib::Cut*> _cuts;

	};

}

#endif

