#ifndef HLIB_CUTTER_H
#define HLIB_CUTTER_H

#include<string>
#include<vector>

#include<cut.hpp>
#include<event.h>

namespace hlib {

	class Cutter {
	
	  public:
		Cutter();

		int get_cutmask(const hlib::Event& event);

		std::string get_abbreviations(int bitmask);

	  private:
		std::vector<hlib::Cut*> _cuts;

	};

}

#endif

