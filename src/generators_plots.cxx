#include<generators_plots.h>

#include<TH1D.h>
#include<TH2D.h>

#include<cutter.h>
#include<object_manager.h>
#include<plotter.h>
#include<template_plot.hpp>
#include<yaml_utils.hpp>

namespace {

	void __cleanDuplicatesFromMap(std::map<std::string, std::vector<long> >& map) {

		for(std::map<std::string, std::vector<long> >::iterator it = map.begin(); it != map.end(); ++it) {
			std::vector<long>& vec = it->second;
			for(unsigned int i = 0; i < vec.size(); ++i) {
				for(unsigned int j = (i+1); j < vec.size(); ++j) {
					if(vec[i] == vec[j]) {
						vec.erase(vec.begin() + j);
					}
				}
			}
		}

	}

	std::map<std::string, std::vector<long> > __mergeMaps(const std::map<std::string, std::vector<long> >& map1, const std::map<std::string, std::vector<long> >& map2) {
		std::map<std::string, std::vector<long> > returnMap = map1;
		for(std::map<std::string, std::vector<long> >::const_iterator it = map2.begin(); it != map2.end(); ++it) {
			std::map<std::string, std::vector<long> >::const_iterator finder = returnMap.find(it->first);
			if(finder == returnMap.end()) {
				returnMap[it->first] = it->second;
			} else {
				const std::vector<long>& vec = it->second;
				returnMap[it->first].insert(returnMap[it->first].end(), vec.begin(), vec.end());
			}
		}
		return returnMap;
	}

	std::map<std::string, std::vector<long> > __getCutmasks(antok::plotUtils::GlobalPlotOptions plotOptions) {
		antok::Cutter& cutter = antok::ObjectManager::instance()->getCutter();
		std::map<std::string, std::vector<long> > returnCutmasks;
		returnCutmasks = __mergeMaps(plotOptions.cutMasks, returnCutmasks);
		if(plotOptions.plotsWithSingleCutsOff) {
			returnCutmasks = __mergeMaps(cutter.getCutmasksAllCutsOffSeparately(), returnCutmasks);
		}
		if(plotOptions.plotsWithSingleCutsOn) {
			returnCutmasks = __mergeMaps(cutter.getCutmasksAllCutsOnSeparately(), returnCutmasks);
		}
		if(plotOptions.plotsForSequentialCuts) {
			returnCutmasks = __mergeMaps(cutter.getWaterfallCutmasks(), returnCutmasks);
		}
		__cleanDuplicatesFromMap(returnCutmasks);
		return returnCutmasks;
	}

}

antok::Plot* antok::generators::generate1DPlot(const YAML::Node& plot, antok::plotUtils::GlobalPlotOptions plotOptions) {

	std::string plotName = antok::YAMLUtils::getString(plot["Name"]);

	double lowerBound = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["LowerBound"], &lowerBound)) {
		std::cerr<<"Could not get \"LowerBound\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	double upperBound = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["UpperBound"], &upperBound)) {
		std::cerr<<"Could not get \"UpperBound\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	int nBins = 0;
	if(not antok::YAMLUtils::getValue<int>(plot["NBins"], &nBins)) {
		std::cerr<<"Could not get \"NBins\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

	std::string variableName = antok::YAMLUtils::getString(plot["Variable"]);
	if(variableName == "") {
		std::cerr<<"\"Variable\" entry invalid for \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

	antok::Data& data = ObjectManager::instance()->getData();
	std::map<std::string, std::vector<long> > cutmasks = __getCutmasks(plotOptions);
	std::string variableType = data.getType(variableName);
	antok::Plot* antokPlot = 0;

	if(variableType == "double") {
		antokPlot = new antok::TemplatePlot<double>(cutmasks, new TH1D(plotName.c_str(), plotName.c_str(), nBins, lowerBound, upperBound), data.getAddr<double>(variableName));
	} else if (variableType == "int") {
		antokPlot = new antok::TemplatePlot<int>(cutmasks, new TH1D(plotName.c_str(), plotName.c_str(), nBins, lowerBound, upperBound), data.getAddr<int>(variableName));
	} else if(variableType == "") {
		std::cerr<<"Could not find \"Variable\" \""<<variableName<<"\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	} else {
		std::cerr<<"\"Variable\"'s type \""<<variableType<<"\" not supported by \"Plot\" (in \""<<plotName<<"\")."<<std::endl;
		return 0;
	}

	return antokPlot;

}

antok::Plot* antok::generators::generate2DPlot(const YAML::Node& plot, antok::plotUtils::GlobalPlotOptions plotOptions) {

	std::string plotName = antok::YAMLUtils::getString(plot["Name"]);

	double lowerBound1 = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["LowerBounds"][0], &lowerBound1)) {
		std::cerr<<"Could not get first of the \"LowerBounds\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	double lowerBound2 = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["LowerBounds"][1], &lowerBound2)) {
		std::cerr<<"Could not get second of the \"LowerBounds\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

	double upperBound1 = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["UpperBounds"][0], &upperBound1)) {
		std::cerr<<"Could not get first of the \"UpperBounds\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	double upperBound2 = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["UpperBounds"][1], &upperBound2)) {
		std::cerr<<"Could not get second of the \"UpperBounds\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

	int nBins1 = 0;
	if(not antok::YAMLUtils::getValue<int>(plot["NBins"][0], &nBins1)) {
		std::cerr<<"Could not get first of the \"NBins\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	int nBins2 = 0;
	if(not antok::YAMLUtils::getValue<int>(plot["NBins"][1], &nBins2)) {
		std::cerr<<"Could not get first of the \"NBins\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

	std::string variable1Name = antok::YAMLUtils::getString(plot["Variables"][0]);
	if(variable1Name == "") {
		std::cerr<<"First of the \"Variables\" entries invalid for \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	std::string variable2Name = antok::YAMLUtils::getString(plot["Variables"][1]);
	if(variable2Name == "") {
		std::cerr<<"Second of the \"Variables\" entries invalid for \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

	antok::Data& data = ObjectManager::instance()->getData();
	std::map<std::string, std::vector<long> > cutmasks = __getCutmasks(plotOptions);

	std::string variableType = data.getType(variable1Name);
	std::string variable2Type = data.getType(variable2Name);
	if(variableType != variable2Type) {
		std::cerr<<"Cannot plot 2D \"Plot\" \""<<plotName<<"\" with \"Variables\" of different type (\""<<variableType<<"\"<>\""<<variable2Type<<"\")."<<std::endl;
		return 0;
	}
	antok::Plot* antokPlot = 0;
	if(variableType == "double") {
		antokPlot = new antok::TemplatePlot<double>(cutmasks,
		                                            new TH2D(plotName.c_str(), plotName.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
		                                            data.getAddr<double>(variable1Name),
		                                            data.getAddr<double>(variable2Name));
	} else if (variableType == "int") {
		antokPlot = new antok::TemplatePlot<int>(cutmasks,
		                                         new TH2D(plotName.c_str(), plotName.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
		                                         data.getAddr<int>(variable1Name),
		                                         data.getAddr<int>(variable2Name));
	} else if(variableType == "") {
		std::cerr<<"Could not find \"Variable\" \""<<variable1Name<<"\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	} else {
		std::cerr<<"\"Variable\"'s type \""<<variableType<<"\" not supported by \"Plot\" (in \""<<plotName<<"\")."<<std::endl;
		return 0;
	}

	return antokPlot;

}

