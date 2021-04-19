#ifndef ANTOK_PLOTTER_H
#define ANTOK_PLOTTER_H

#include<map>
#include<string>
#include<vector>

#include<TH1.h>

namespace YAML {
	class Node;
}

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
			std::map<std::string, std::vector<long>> cutMasks;

		  private:

		};

		struct waterfallHistogramContainer {

			waterfallHistogramContainer(TH1*                                                    histogram,
			                            const std::vector<std::pair<std::string, const bool*>>& cuts)
				: _histogram(histogram),
				  _cuts(cuts)
			{
				TAxis* xAxis = _histogram->GetXaxis();
				int startBin = 1;
				for (int i = 1; i <= _histogram->GetNbinsX(); ++i) {
					if (std::string(xAxis->GetBinLabel(i)) == "") {
						startBin = i;
						break;
					}
				}
				const int nBinsNeeded = startBin + _cuts.size() - 1;
				if (nBinsNeeded > _histogram->GetNbinsX()) {
					_histogram->SetBins(nBinsNeeded,
															xAxis->GetBinLowEdge(1),
															xAxis->GetBinUpEdge(startBin - 1) + xAxis->GetBinWidth(startBin - 1) * (nBinsNeeded - (startBin - 1)));
				}
				for (size_t i = 0; i < _cuts.size(); ++i) {
					xAxis->SetBinLabel(startBin + i, _cuts[i].first.c_str());
				}
			};

			TH1*                                             _histogram;
			std::vector<std::pair<std::string, const bool*>> _cuts;

		};

	}

	class Plotter {

		friend class Initializer;

	  public:

		static Plotter* instance();

		void fill(const long& cutPattern);
		void addInputfileToWaterfallHistograms(const TH1D* waterfall);

		static bool handleAdditionalCuts(const YAML::Node& cuts, std::map<std::string, std::vector<long>>& map);

	  private:

		Plotter() { };

		static Plotter* _plotter;

		std::vector<antok::Plot*> _plots;

		std::vector<antok::plotUtils::waterfallHistogramContainer> _waterfallHistograms;

	};

}

#endif

