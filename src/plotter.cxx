#include<plotter.h>

#include<yaml-cpp/yaml.h>

#include<TH1.h>
#include<stdexcept>

#include<object_manager.h>
#include<plot.hpp>
#include<yaml_utils.hpp>

using antok::YAMLUtils::handleOnOffOption;

namespace {

	// Inserts a new bin with the given label at bin i_bos to the given histogram.
	// All bins, starting from binIndex, will be moved one bin to the right
	void insertLabeledBinInHist(TH1* h, const char* label, const int binIndex);

}


antok::Plotter* antok::Plotter::_plotter = nullptr;

antok::Plotter* antok::Plotter::instance() {
	if(_plotter == nullptr) {
		_plotter = new antok::Plotter();
	}
	return _plotter;
}


void
antok::Plotter::fill(const antok::bitmask& cutPattern)
{
	for (size_t i = 0; i < _plots.size(); ++i) {
		_plots[i]->fill(cutPattern);
	}
	for (size_t iCutTrain = 0; iCutTrain < _waterfallHistograms.size(); ++iCutTrain) {
		TH1* hist = _waterfallHistograms[iCutTrain]._histogram;
		const std::vector<std::pair<std::string, const bool*>>& cuts = _waterfallHistograms[iCutTrain]._cuts;
		for (size_t iCut = 0; iCut < cuts.size(); ++iCut) {
			if (*(cuts[iCut].second)) {
				hist->Fill(cuts[iCut].first.c_str(), 1);
			} else {
				break;
			}
		}
	}
}


void
antok::Plotter::addInputfileToWaterfallHistograms(const TH1D* inFileWaterfall)
{
	for (std::vector<antok::plotUtils::waterfallHistogramContainer>::iterator w = _waterfallHistograms.begin();
	     w != _waterfallHistograms.end(); ++w) {
	  TH1* existingWaterfall = w->_histogram;
		for (int iBin = 1; iBin <= inFileWaterfall->GetNbinsX(); ++iBin) {
			const char* label = inFileWaterfall->GetXaxis()->GetBinLabel(iBin);
			if (std::string(label) != "") {
#if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
				int iBinLabel = existingWaterfall->GetXaxis()->FindFixBin(label);
#else
				existingWaterfall->SetBit(TH1::kCanRebin, false);
				int iBinLabel = existingWaterfall->GetXaxis()->FindBin(label);
#endif
				if (iBinLabel < 1) {  // cannot find label from new histogram in existing waterfall histogram
					insertLabeledBinInHist(existingWaterfall, label, iBin);
					iBinLabel = iBin;
				}
				if (iBinLabel != iBin) {
					std::cerr << iBinLabel << " != " << iBin;
					throw std::logic_error("Order of bins in waterfall histogram changes from file to file.");
				}
				existingWaterfall->SetBinContent(iBinLabel, existingWaterfall->GetBinContent(iBinLabel) + inFileWaterfall->GetBinContent(iBin));
			}
		}
	}
}


namespace {

	std::string __getCutnamesOffNode(const YAML::Node& withCut, const std::string& cutTrainName) {

		std::string cutName = antok::YAMLUtils::getString(withCut["ShortName"]);
		if(cutName == "") {
			std::cerr<<"Warning: One of the \"WithCuts\" or \"WithoutCuts\" appears to have a missing or invalid ";
			std::cout<<"\"ShortName\" for \"CutTrain\" \""<<cutTrainName<<"\", skipping entry."<<std::endl;
			return "";
		}
		antok::Cutter& cutter = antok::ObjectManager::instance()->getCutter();
		if(cutter.cutInCutTrain(cutName, cutTrainName)) {
			return cutName;
		} else {
			std::cerr<<"Warning: \"WithCuts\" or \"WithoutCuts\" entry \""<<cutName;
			std::cout<<"\" is not in \"CutTrain\" \""<<cutTrainName<<"\", skipping entry."<<std::endl;
			return "";
		}

	}

	antok::bitmask __handleCutList(const YAML::Node& cutNode, const std::string& cutTrainName, bool invertSelection) {

		std::vector<std::string> cutNames;
		if(cutNode.IsMap()) {
			std::string cutName = __getCutnamesOffNode(cutNode, cutTrainName);
			if(cutName == "") {
				return -1;
			}
			cutNames.push_back(cutName);
		} else if (cutNode.IsSequence()) {
			for(YAML::const_iterator innerCuts_it = cutNode.begin(); innerCuts_it != cutNode.end(); ++innerCuts_it) {
				std::string cutName = __getCutnamesOffNode(*innerCuts_it, cutTrainName);
				if(cutName == "") {
					cutNames.clear();
					break;
				}
				cutNames.push_back(cutName);
			}
			if(cutNames.empty()) {
				return -1;
			}
		} else {
			std::cerr<<"Warning: One of the \"WithCuts\" or \"WithoutCuts\" appears to have invalid format for \"CutTrain\" \""<<cutTrainName<<"\", skipping entry."<<std::endl;
			return -1;
		}
		antok::Cutter& cutter = antok::ObjectManager::instance()->getCutter();
		antok::bitmask cutmask = cutter.getCutmaskForNames(cutNames);
		if(invertSelection) {
			antok::bitmask allCuts = cutter.getAllCutsCutmaskForCutTrain(cutTrainName);
			cutmask = (~cutmask)&allCuts;
		}
		return cutmask;
	}

}

bool antok::Plotter::handleAdditionalCuts(const YAML::Node& trainList, std::map<std::string, std::vector<antok::bitmask> >& map) {

	using antok::YAMLUtils::hasNodeKey;

	bool error = false;
	for(YAML::const_iterator trainList_it = trainList.begin(); trainList_it != trainList.end(); ++trainList_it) {
		const YAML::Node& entry = *trainList_it;
		if((not (hasNodeKey(entry, "CutTrain"))) or (not hasNodeKey(entry["CutTrain"], "Name"))) {
			std::cerr<<"Warning: \"CutTrain\" missing or invalid for an entry, skipping it."<<std::endl;
			continue;
		}
		std::string cutTrainName = antok::YAMLUtils::getString(entry["CutTrain"]["Name"]);
		if(cutTrainName == "") {
			std::cerr<<"Warning: Could not convert \"CutTrain\" to std::string for an entry, skipping it."<<std::endl;
			continue;
		}
		if(not (hasNodeKey(entry, "WithCuts") and hasNodeKey(entry, "WithoutCuts"))) {
			std::cerr<<"Warning: Either \"WithCuts\" or \"WithoutCuts\" missing for \"CutTrain\" \""<<cutTrainName<<"\", skipping entry."<<std::endl;
			continue;
		}
		if(not (entry["WithCuts"].IsSequence() and entry["WithoutCuts"].IsSequence())) {
			std::cerr<<"Warning: Either \"WithCuts\" or \"WithoutCuts\" is not a sequence for \"CutTrain\" \""<<cutTrainName<<"\", skipping entry."<<std::endl;
			continue;
		}
		bool innerError = false;
		std::vector<antok::bitmask> cutMasks;
		for(YAML::const_iterator withCuts_it = entry["WithCuts"].begin(); withCuts_it != entry["WithCuts"].end(); ++withCuts_it) {
			const YAML::Node& withCut = *withCuts_it;
			antok::bitmask cutmask = __handleCutList(withCut, cutTrainName, false);
			cutMasks.push_back(cutmask);
		}
		if(innerError) {
			continue;
		}
		for(YAML::const_iterator withoutCuts_it = entry["WithoutCuts"].begin(); withoutCuts_it != entry["WithoutCuts"].end(); ++withoutCuts_it) {
			const YAML::Node& withoutCut = *withoutCuts_it;
			antok::bitmask cutmask = __handleCutList(withoutCut, cutTrainName, true);
			cutMasks.push_back(cutmask);
		}
		if(innerError) {
			continue;
		}
		map[cutTrainName] = cutMasks;
	}
	return (not error);

}

antok::plotUtils::GlobalPlotOptions::GlobalPlotOptions(const YAML::Node& optionNode) {

	using antok::YAMLUtils::hasNodeKey;

	plotsForSequentialCuts = false;
	plotsWithSingleCutsOn = false;
	plotsWithSingleCutsOff = false;
	statisticsHistInName = "";
	statisticsHistOutName = "";

	if(not optionNode) {
		std::cerr<<"Warning: \"GlobalPlotOptions\" not found in configuration file."<<std::endl;
		return;
	}

	plotsForSequentialCuts = handleOnOffOption("PlotsForSequentialCuts", optionNode, "GobalPlotOptions");
	plotsWithSingleCutsOn = handleOnOffOption("PlotsWithSingleCutsOn", optionNode, "GobalPlotOptions");
	plotsWithSingleCutsOff = handleOnOffOption("PlotsWithSingleCutsOff", optionNode, "GobalPlotOptions");

	if(not hasNodeKey(optionNode, "StatisticsHistogram")) {
		std::cerr<<"Warning: \"StatisticsHistogram\" not found in \"GobalPlotOptions\", switching it off"<<std::endl;
	} else {
		const YAML::Node& statsHistOpt = optionNode["StatisticsHistogram"];
		bool state = handleOnOffOption("State", statsHistOpt, "StatisticsHistogram");
		if(state) {
			if(not (hasNodeKey(statsHistOpt, "InputName") and hasNodeKey(statsHistOpt, "OutputName"))) {
				std::cerr<<"Warning: \"InputName\" or \"OutputName\" not found in \"GobalPlotOptions\"'s \"StatisticsHistogram\", switching histogram off."<<std::endl;
			} else {
				statisticsHistInName = antok::YAMLUtils::getString(statsHistOpt["InputName"]);
				statisticsHistOutName = antok::YAMLUtils::getString(statsHistOpt["OutputName"]);
				if(statisticsHistInName == "" or statisticsHistOutName == "") {
					std::cerr<<"Warning: Could not convert either \"InputName\" or \"OutputName\" to std::string, switching histogram off."<<std::endl;
					statisticsHistInName = "";
					statisticsHistOutName = "";
				}
			}
		}
	}

	if(hasNodeKey(optionNode, "GlobalCuts")) {
		if(not antok::Plotter::handleAdditionalCuts(optionNode["GlobalCuts"], cutMasks)) {
			std::cerr<<"Warning: There was a problem when processing the \"GlobalCuts\" in \"GobalPlotOptions\"."<<std::endl;
		}
	}

	if(hasNodeKey(optionNode, "HistogramNameAppendix")) {
		std::string histNameAppendix = antok::YAMLUtils::getString(optionNode["HistogramNameAppendix"]);
		if(histNameAppendix == "") {
			std::cerr<<"Warning: \"HistogramNameAppendix\" in \"GlobalPlotOptions\" is empty or not a valid string."<<std::endl;
		}
		antok::ObjectManager::instance()->registerHistogramNameAppendix(histNameAppendix);
	}

}


namespace {

	// Inserts a new bin with the given label at bin binIndex to the given histogram.
	// All bins, starting from binIndex, will be moved one bin to the right
	void
	insertLabeledBinInHist(TH1*        h,
	                       const char* label,
	                       const int   binIndex)
	{
		int nBins = h->GetNbinsX();
		TAxis* xAxis = h->GetXaxis();

		// last bin is not free --> increase number of bins by one
		if (std::string(xAxis->GetBinLabel(nBins)) != "") {
			h->SetBins(nBins + 1,
			           h->GetXaxis()->GetBinLowEdge(1),
			           h->GetXaxis()->GetBinUpEdge(nBins) + h->GetXaxis()->GetBinWidth(nBins));
			nBins += 1;
		}

		for (int i = nBins - 1 ; i >= binIndex; --i) {
			const char* binLabel = h->GetXaxis()->GetBinLabel(i);
			if (std::string(binLabel) != "") {
				xAxis->SetBinLabel(i + 1, binLabel);
				xAxis->SetBinLabel(i, "");
				h->SetBinContent(i + 1, h->GetBinContent(i));
				h->SetBinContent(i, 0);
			}
		}
		xAxis->SetBinLabel(binIndex, label);
	}

}
