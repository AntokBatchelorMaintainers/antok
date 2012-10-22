#ifndef ANTOK_PLOT_H
#define ANTOK_PLOT_H

#include<map>
#include<vector>

#include<TH1.h>

#include<event.h>

namespace antok {

	class Plot {

	  public:

		Plot(std::map<std::string, std::vector<long> >& cutmasks, TH1* hist_template, double* data1, double* data2 = 0);

		void fill(long cutmask);

		int getNoCutsForCutTrain(std::string cutTrainName);

	  private:
		TH1* _histTemplate;
		std::vector<std::pair<TH1*, long> > _histograms;

		double* _data1;
		double* _data2;

	};

}

#endif
