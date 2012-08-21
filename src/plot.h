#ifndef HLIB_PLOT_H
#define HLIB_PLOT_H

#include<vector>

#include<TH1.h>

#include<event.h>

namespace hlib {

	class Plot {

	  public:

		Plot(std::vector<int> cutmasks, TH1* hist_template, double* data1, double* data2 = NULL);

		void fill(int cutmask);

		std::vector<TH1*> get_histograms() { return _histograms; };

		const TH1* get_template() const { return _hist_template; };

	  private:
		std::vector<int> _cutmasks;
		TH1* _hist_template;
		std::vector<TH1*> _histograms;

		double* _data1;
		double* _data2;

	};

}

#endif
