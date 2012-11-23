#ifndef ANTOK_PLOT_HPP
#define ANTOK_PLOT_HPP

namespace antok {

	class Plot {

	  public:

		Plot() { };
		virtual ~Plot() { };
		virtual void fill(long cutmask) = 0;
		std::map<std::string, std::map<unsigned int, TH1*> > getHistogramOrder() const { return _histogramOrder; };

	  protected:
		std::map<std::string, std::map<unsigned int, TH1*> > _histogramOrder;

	};

}

#endif

