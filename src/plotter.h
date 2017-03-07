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
			std::map<std::string, std::vector<long> > cutMasks;

		  private:

		};

		struct waterfallHistogramContainer {

			waterfallHistogramContainer(TH1* hist,
			                            std::vector<std::pair<const char*, const bool*> > cuts_)
				: histogram(hist),
				  cuts(cuts_) {
						int startBin = 1;
						for(int i = 1; i <= hist->GetNbinsX(); ++i) {
							if(std::string(hist->GetXaxis()->GetBinLabel(i)) == "") {
								startBin = i;
								break;
							}
						}
						const int nBinsNeeded = startBin + cuts.size()-1;
						if(nBinsNeeded > hist->GetNbinsX()){
							histogram->SetBins( nBinsNeeded, hist->GetXaxis()->GetBinLowEdge(1),
							                    hist->GetXaxis()->GetBinUpEdge(startBin - 1) + hist->GetXaxis()->GetBinWidth(startBin - 1) * (nBinsNeeded-(startBin-1)));
						}
						for(unsigned int i = 0; i < cuts.size(); ++i) {
							histogram->GetXaxis()->SetBinLabel(startBin + i, cuts[i].first);
						}
					};

			TH1* histogram;
			std::vector<std::pair<const char*, const bool*> > cuts;

		};

	}

	class Plotter {

		friend class Initializer;

	  public:

		static Plotter* instance();

		void fill(const long& cutPattern);
		void addInputfileToWaterfallHistograms(const TH1D* waterfall);

		static bool handleAdditionalCuts(const YAML::Node& cuts, std::map<std::string, std::vector<long> >& map);

	  private:

		Plotter() { };

		static Plotter* _plotter;

		std::vector<antok::Plot*> _plots;

		std::vector<antok::plotUtils::waterfallHistogramContainer> _waterfallHistograms;

	};

}

#endif

