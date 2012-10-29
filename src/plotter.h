#ifndef ANTOK_PLOTTER_H
#define ANTOK_PLOTTER_H

#include<map>
#include<string>
#include<vector>

namespace YAML {
	class Node;
}

class TH1;

namespace antok {

	class Plot;

	namespace plotUtils {

		struct GlobalPlotOptions {

			GlobalPlotOptions(const YAML::Node& optionNode);
			bool plotsForSequentialCuts;
			bool plotsWithSingleCutsOn;
			bool plotsWithSingleCutsOff;
			std::string statisticsHistInName;
			std::string statisticsHistOutName;
			std::map<std::string, std::vector<long> > cutMasks;

		  private:

			bool handleOnOffOption(std::string optionName, const YAML::Node& option, std::string location) const;

		};

		struct waterfallHistogramContainer {

			waterfallHistogramContainer(TH1* hist,
			                            std::vector<std::pair<const char*, const bool*> > cuts_)
				: histogram(hist),
				  cuts(cuts_) { };

			TH1* histogram;
			std::vector<std::pair<const char*, const bool*> > cuts;

		};

	}


	class Plotter {

		friend class Initializer;

	  public:

		static Plotter* instance();

		void fill(long cutPattern);

		static bool handleAdditionalCuts(const YAML::Node& cuts, std::map<std::string, std::vector<long> >& map);

	  private:

		Plotter();

		static Plotter* _plotter;

		std::vector<antok::Plot*> _plots;

		std::vector<antok::plotUtils::waterfallHistogramContainer> _waterfallHistograms;

	};

}

#endif

