#include<plotter.h>

#include<yaml-cpp/yaml.h>

#include<cutter.h>
#include<object_manager.h>
#include<plot.hpp>
#include<yaml_utils.hpp>

antok::Plotter* antok::Plotter::_plotter = 0;

antok::Plotter* antok::Plotter::instance() {
	if(_plotter == 0) {
		_plotter = new antok::Plotter();
	}
	return _plotter;
}

antok::Plotter::Plotter() {

/*
	std::vector<int> cutmasks;
	std::vector<int> standard_cutmasks;

	standard_cutmasks.push_back(0);
	standard_cutmasks.push_back(256);
	standard_cutmasks.push_back(384);
	standard_cutmasks.push_back(448);
	standard_cutmasks.push_back(480);
	standard_cutmasks.push_back(496);
	standard_cutmasks.push_back(504);
	standard_cutmasks.push_back(508);
	standard_cutmasks.push_back(510);
	standard_cutmasks.push_back(511);

	standard_cutmasks.push_back(48);
	standard_cutmasks.push_back(56);
	standard_cutmasks.push_back(112);
	standard_cutmasks.push_back(120);
	standard_cutmasks.push_back(176);
	standard_cutmasks.push_back(184);
	standard_cutmasks.push_back(240);
	standard_cutmasks.push_back(248);

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	_plots.push_back(antok::Plot(cutmasks, new TH1D("mass_5pi", "mass_5Pi", 500, 0, 7), &XMass));
	antok::Data& data = antok::ObjectManager::instance()->getData();
	_plots.push_back(antok::Plot(cutmasks, new TH1D("mass_5pi_newWay", "mass_5Pi_newWay", 500, 0, 7), data.getDoubleAddr("XMass")));
	cutmasks.clear();



	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(128);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("mom_5pi", "mom_5Pi", 500, 0, 250), &XMom));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(128);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("calc_beam_E", "calc_beam_E", 500, 0, 250), &CalcBeamE));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(8);
	cutmasks.push_back(264);
	cutmasks.push_back(136);
	cutmasks.push_back(392);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("rpd_mult", "rpd_mult", 10, 0, 10), &RPDMult));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(16);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("rpd_P_mass", "rpd_P_mass", 1000, 0, 10), &ProtonMass));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(4);
	_plots.push_back(antok::Plot(cutmasks, new TH2D("vtx_pos", "vtx_pos", 1000, -5, 5, 1000, -5, 5), &PrimVX, &PrimVY));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(2);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("vtx_z", "vtx_z", 2000, -200, 200), &PrimVZ));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(64);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("t_prime", "t_prime", 1000, -5, 5), &TPrim));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(128);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi", "delta_phi", 500, -7, 7), &RPDDeltaPhi));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(1);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("trigger_mask", "trigger_mask", 15, 0, 15), &TrigMask));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	_plots.push_back(antok::Plot(cutmasks, new TH1D("beam_time", "beam_time", 100, -10, 10), &beam_time));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	_plots.push_back(antok::Plot(cutmasks, new TH2D("cedar_theta", "cedar_theta", 3000, -300, 300, 3000, -300, 300), &cedarTheta_X, &cedarTheta_Y));
	cutmasks.clear();

	cutmasks.push_back(384);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi", "delta_phi", 500, -7, 7), &RPDDeltaPhi));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi_fhaas", "delta_phi", 500, -7, 7), &RPDDeltaPhi_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH2D("delta_phi_comparison", "delta_phi_comparison", 1000, -7, 7, 1000, -7, 7), &RPDDeltaPhi, &RPDDeltaPhi_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi_abs", "delta_abs_phi", 500, -7, 7), &RPDDeltaPhiAbs));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi_abs_fhaas", "delta_abs_phi", 500, -7, 7), &RPDDeltaPhiAbs_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH2D("delta_phi_abs_comparison", "delta_phi_abs_comparison", 1000, -7, 7, 1000, -7, 7), &RPDDeltaPhiAbs, &RPDDeltaPhiAbs_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("phi_res", "phi_res", 500, -7, 7), &RPDPhiRes));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("phi_res_fhaas", "phi_res_fhaas", 500, -7, 7), &RPDPhiRes_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH2D("phi_res_comparison", "phi_res_comparison", 500, -7, 7, 500, -7, 7), &RPDPhiRes, &RPDPhiRes_fhaas));
*/
};

void antok::Plotter::fill(long cutPattern) {

	for(unsigned int i = 0; i < _plots.size(); ++i) {
		_plots[i]->fill(cutPattern);
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

	long __handleCutList(const YAML::Node& cutNode, const std::string& cutTrainName, bool invertSelection) {

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
			if(cutNames.size() == 0) {
				return -1;
			}
		} else {
			std::cerr<<"Warning: One of the \"WithCuts\" or \"WithoutCuts\" appears to have invalid format for \"CutTrain\" \""<<cutTrainName<<"\", skipping entry."<<std::endl;
			return -1;
		}
		antok::Cutter& cutter = antok::ObjectManager::instance()->getCutter();
		long cutmask = cutter.getCutmaskForNames(cutNames);
		if(invertSelection) {
			long allCuts = cutter.getAllCutsCutmaskForCutTrain(cutTrainName);
			cutmask = (~cutmask)&allCuts;
		}
		return cutmask;
	}

}

bool antok::Plotter::handleAdditionalCuts(const YAML::Node& trainList, std::map<std::string, std::vector<long> >& map) {

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
		std::vector<long> cutMasks;
		for(YAML::const_iterator withCuts_it = entry["WithCuts"].begin(); withCuts_it != entry["WithCuts"].end(); ++withCuts_it) {
			const YAML::Node& withCut = *withCuts_it;
			long cutmask = __handleCutList(withCut, cutTrainName, false);
			if(cutmask < 0) {
				innerError = true;
				error = true;
				break;
			}
			cutMasks.push_back(cutmask);
		}
		if(innerError) {
			continue;
		}
		for(YAML::const_iterator withoutCuts_it = entry["WithoutCuts"].begin(); withoutCuts_it != entry["WithoutCuts"].end(); ++withoutCuts_it) {
			const YAML::Node& withoutCut = *withoutCuts_it;
			long cutmask = __handleCutList(withoutCut, cutTrainName, true);
			if(cutmask < 0) {
				innerError = true;
				error = true;
				break;
			}
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

	if(not hasNodeKey(optionNode, "StaticsticsHistogram")) {
		std::cerr<<"Warning: \"StaticsticsHistogram\" not found in \"GobalPlotOptions\", switching it off"<<std::endl;
	} else {
		const YAML::Node& statsHistOpt = optionNode["StaticsticsHistogram"];
		bool state = handleOnOffOption("State", statsHistOpt, "StaticsticsHistogram");
		if(state) {
			if(not (hasNodeKey(statsHistOpt, "InputName") and hasNodeKey(statsHistOpt, "OutputName"))) {
				std::cerr<<"Warning: \"InputName\" or \"OutputName\" not found in \"GobalPlotOptions\"'s \"StaticsticsHistogram\", switching histogram off."<<std::endl;
			} else {
				statisticsHistInName = antok::YAMLUtils::getString(statsHistOpt["InputName"]);
				statisticsHistOutName = antok::YAMLUtils::getString(statsHistOpt["OutputName"]);
				if(statisticsHistInName == "" or statisticsHistOutName == "") {
					std::cerr<<"Warning: Could not convert either \"InputName\" or \"OutputName\" to std::string, switching histogram off."<<std::endl;
				}
			}
		}
	}

	if(not hasNodeKey(optionNode, "GlobalCuts")) {
		std::cerr<<"Warning: \"GlobalCuts\" not found in \"GobalPlotOptions\", not adding additional cuts."<<std::endl;
	} else {
		if(not antok::Plotter::handleAdditionalCuts(optionNode["GlobalCuts"], cutMasks)) {
			std::cerr<<"Warning: There was a problem when processing the \"GlobalCuts\" in \"GobalPlotOptions\"."<<std::endl;
		}
	}

}

bool antok::plotUtils::GlobalPlotOptions::handleOnOffOption(std::string optionName, const YAML::Node& option, std::string location) const {

	using antok::YAMLUtils::hasNodeKey;

	if(not hasNodeKey(option, optionName)) {
		std::cerr<<"Warning: \""<<optionName<<"\" not found in \""<<location<<"\", switching it off"<<std::endl;
	} else {
		std::string optionValue = antok::YAMLUtils::getString(option[optionName]);
		if(optionValue == "On") {
			return true;
		} else if (optionValue == "Off") {
			// returning at end of function
		} else {
			std::cerr<<"Warning: \""<<location<<"\"'s \""<<optionName<<"\" is \""<<optionValue<<"\" instead of \"On\" or \"Off\", switching it off"<<std::endl;
		}
	}
	return false;

}

