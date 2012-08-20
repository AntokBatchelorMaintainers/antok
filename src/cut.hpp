#ifndef HLIB_CUT_H
#define HLIB_CUT_H

namespace hlib {

	class Cut {
	
	  public:

	    Cut() { };

		virtual bool apply() const = 0;

		const std::string shortname = "ERROR";
		const std::string longname = "ERROR";
		const std::string abbreviation = "ERROR";
	
	}

	class TrigMask : Cut {

	  public:
		
		bool apply() const { return !(data.TrigMask&0x1); };



	}



}

#endif

