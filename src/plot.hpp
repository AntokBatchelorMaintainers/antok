#ifndef ANTOK_PLOT_HPP
#define ANTOK_PLOT_HPP

namespace antok {

	class Plot {

	  public:

		virtual ~Plot() { };
		virtual void fill(antok::bitmask cutmask) = 0;

	};

}

#endif

