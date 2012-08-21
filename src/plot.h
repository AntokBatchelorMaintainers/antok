#ifndef HLIB_PLOT_H
#define HLIB_PLOT_H

#include<vector>

#include<TH1.h>

#include<event.h>

namespace hlib {

	class Plot {

		Plot(std::vector<int> cutmasks, TH1* hist_template, int dim = 1);

		void fill(int cutmask, double data1, double data2 = -999999999.);

		std::vector<TH1*> get_histograms();

	  private:
		std::vector<int> _cutmasks;
		std::vector<TH1*> _histograms;
		int _dim;

	};

}

#endif
