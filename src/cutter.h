#ifndef HLIB_CUTTER_H
#define HLIB_CUTTER_H

#include<string>
#include<vector>

namespace hlib {

	class Cutter {
	
	  public:
		Cutter() { };

		virtual int get_cutmask(hlib::Event event);

		std::string get_abbreviations(int bitmask);

	  private:
		std::vector<hlib::Cut*> _cuts;

	};

}

#endif

