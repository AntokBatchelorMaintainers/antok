#ifndef ANTOK_PLOTTER_H
#define ANTOK_PLOTTER_H

#include<vector>

namespace antok {

	class Plot;

	class Plotter {

		friend class Initializer;

	  public:

		static Plotter* instance();

		void fill(long cutPattern);

	  private:

		Plotter();

		static Plotter* _plotter;

		std::vector<antok::Plot*> _plots;

	};

}

#endif
