#include<generators_plots.h>

#include<TH1D.h>
#include<TH2D.h>

#include<cutter.h>
#include<object_manager.h>
#include<template_plot.hpp>
#include<yaml_utils.hpp>

antok::Plot* antok::generators::generate1DPlot(const YAML::Node& plot) {

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
	antok::Cutter& cutter = ObjectManager::instance()->getCutter();
	std::map<std::string, std::vector<long> > cutmasks = cutter.getWaterfallCutmasks();
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

antok::Plot* antok::generators::generate2DPlot(const YAML::Node& plot) {

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
	antok::Cutter& cutter = ObjectManager::instance()->getCutter();
	std::map<std::string, std::vector<long> > cutmasks = cutter.getWaterfallCutmasks();

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

