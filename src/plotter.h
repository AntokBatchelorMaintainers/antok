#ifndef HLIB_PLOTTER_H
#define HLIB_PLOTTER_H

#include<vector>

#include<TH1D.h>
#include<TH2D.h>

#include<event.h>
#include<plot.h>

namespace hlib {

	class Plotter {

	  public:
		Plotter();

		void fill(const hlib::Event& event);


	  private:

		std::vector<hlib::Plot> _plots;

	};

}

#endif
