#include<generators_plots.h>

#include<TH1D.h>
#include<TH2D.h>
#include<TH3D.h>

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

	std::map<std::string, std::vector<long> > __mergeMaps(const std::map<std::string, std::vector<long> >& map1,
	                                                      const std::map<std::string, std::vector<long> >& map2)
	{
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

	void __getCutmasks(const antok::plotUtils::GlobalPlotOptions& plotOptions,
	                   const YAML::Node& plot,
	                   std::map<std::string, std::vector<long> >& cutmasks)
	{
		antok::Cutter& cutter = antok::ObjectManager::instance()->getCutter();
		cutmasks = __mergeMaps(plotOptions.cutMasks, cutmasks);

		if (antok::YAMLUtils::hasNodeKey(plot, "PlotsWithSingleCutsOff")) {
			const bool plotsWithSingleCutsOff = antok::YAMLUtils::handleOnOffOption("PlotsWithSingleCutsOff", plot, plot["Name"].as<std::string>());
			if (plotsWithSingleCutsOff) {
				cutmasks = __mergeMaps(cutter.getCutmasksAllCutsOffSeparately(), cutmasks);
			}
		} else if (plotOptions.plotsWithSingleCutsOff) {
			cutmasks = __mergeMaps(cutter.getCutmasksAllCutsOffSeparately(), cutmasks);
		}

		if (antok::YAMLUtils::hasNodeKey(plot, "PlotsWithSingleCutsOn")) {
			const bool plotsWithSingleCutsOn = antok::YAMLUtils::handleOnOffOption("PlotsWithSingleCutsOn", plot, plot["Name"].as<std::string>());
			if (plotsWithSingleCutsOn) {
				cutmasks = __mergeMaps(cutter.getCutmasksAllCutsOnSeparately(), cutmasks);
			}
		} else if (plotOptions.plotsWithSingleCutsOn) {
			cutmasks = __mergeMaps(cutter.getCutmasksAllCutsOnSeparately(), cutmasks);
		}

		if (antok::YAMLUtils::hasNodeKey(plot, "PlotsForSequentialCuts")) {
			const bool plotsForSequentialCuts = antok::YAMLUtils::handleOnOffOption("PlotsForSequentialCuts", plot, plot["Name"].as<std::string>());
			if (plotsForSequentialCuts) {
				cutmasks = __mergeMaps(cutter.getWaterfallCutmasks(), cutmasks);
			}

		} else if (plotOptions.plotsForSequentialCuts) {
			cutmasks = __mergeMaps(cutter.getWaterfallCutmasks(), cutmasks);
		}

		__cleanDuplicatesFromMap(cutmasks);
	}

	std::vector<int> __getIndices(const YAML::Node& plot, std::string plotName) {
		std::vector<int> indices;
		try {
			indices = plot["Indices"].as<std::vector<int> >();
		} catch (const YAML::TypedBadConversion<std::vector<int> >& e) {
			std::cerr<<"Could not convert YAML sequence to std::vector<int> when parsing \"Indices\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		} catch (const YAML::TypedBadConversion<int>& e) {
			std::cerr<<"Could not convert entries in YAML sequence to int when parsing \"Indices\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		}
		return indices;
	}

	template<typename T>
	std::vector<T*>* __getDataVector(const YAML::Node& plot, std::string plotName, std::string variableName, const std::vector<int>& indices) {

		std::vector<T*>* vecData = new std::vector<T*>();
		antok::Data& data = antok::ObjectManager::instance()->getData();
		for(unsigned int i = 0; i < indices.size(); ++i) {
			std::stringstream strStr;
			strStr<<variableName<<indices[i];
			T* addr = data.getAddr<T>(strStr.str());
			if(addr == 0) {
				std::cerr<<"All variables for \"Plot\" with \"Indices\" need to have the same type (in \""<<plotName<<"\")."<<std::endl;
				return 0;
			}
			vecData->push_back(addr);
		}
		if(vecData->empty()) {
			std::cerr<<"Got no indices when trying to compile data for \"Plot\" \""<<plotName<<"\"."<<std::endl;
			return 0;
		}
		return(vecData);

	}

	std::string __getPlotNameWithAxisLabels(const YAML::Node& plot) {

		std::string plotName = antok::YAMLUtils::getString(plot["Name"]);
		std::string plotNameWithAxis = plotName;
		if(antok::YAMLUtils::hasNodeKey(plot, "AxisLabels")) {
			std::string xAxisLabel = "";
			std::string yAxisLabel = "";
			if(not antok::YAMLUtils::getValue<std::string>(plot["AxisLabels"][0], &xAxisLabel)) {
				std::cerr<<"Could not get first of the \"AxisLabels\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
				return "";
			}
			if(not antok::YAMLUtils::getValue<std::string>(plot["AxisLabels"][1], &yAxisLabel)) {
				std::cerr<<"Could not get second of the \"AxisLabels\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
				return "";
			}
			std::stringstream strStr;
			strStr<<plotNameWithAxis<<';'<<xAxisLabel<<';'<<yAxisLabel;
			plotNameWithAxis = strStr.str();
		}

		return plotNameWithAxis;

	}

}

antok::Plot* antok::generators::generate1DPlot(const YAML::Node& plot, const antok::plotUtils::GlobalPlotOptions& plotOptions) {

	using antok::YAMLUtils::hasNodeKey;

	std::string plotName = antok::YAMLUtils::getString(plot["Name"]);
	std::string plotNameWithAxisLables = __getPlotNameWithAxisLabels(plot);
	if(plotNameWithAxisLables == "") {
		std::cerr<<"Could not get plot name with axis labels in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

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
	if(lowerBound >= upperBound){
		std::cerr<<"\"LowerBound\" >= \"UpperBound\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
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
	std::map<std::string, std::vector<long> > cutmasks;
	if(hasNodeKey(plot, "CustomCuts")) {
		if(not antok::Plotter::handleAdditionalCuts(plot["CustomCuts"], cutmasks)) {
			std::cerr<<"Warning: There was a problem when processing the \"CustomCuts\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		}
	}
	__getCutmasks(plotOptions, plot, cutmasks);

	antok::Plot* antokPlot = 0;

	if(not hasNodeKey(plot, "Indices")) {

		std::string variableType = data.getType(variableName);
		if(variableType == "double") {
			antokPlot = new antok::TemplatePlot<double>(cutmasks,
			                                            new TH1D(plotName.c_str(),
			                                                     plotNameWithAxisLables.c_str(),
			                                                     nBins,
			                                                     lowerBound,
			                                                     upperBound),
			                                            data.getAddr<double>(variableName));
		} else if (variableType == "int") {
			antokPlot = new antok::TemplatePlot<int>(cutmasks,
			                                         new TH1D(plotName.c_str(),
			                                                  plotNameWithAxisLables.c_str(),
			                                                  nBins,
			                                                  lowerBound,
			                                                  upperBound),
			                                         data.getAddr<int>(variableName));
		} else if (variableType == "std::vector<int>") {
			antokPlot = new antok::TemplatePlot<int>(cutmasks,
			                                            new TH1D(plotName.c_str(),
			                                                     plotNameWithAxisLables.c_str(),
			                                                     nBins,
			                                                     lowerBound,
			                                                     upperBound),
			                                            data.getAddr<std::vector<int> >(variableName));
		} else if (variableType == "std::vector<double>") {
			antokPlot = new antok::TemplatePlot<double>(cutmasks,
			                                            new TH1D(plotName.c_str(),
			                                                     plotNameWithAxisLables.c_str(),
			                                                     nBins,
			                                                     lowerBound,
			                                                     upperBound),
			                                            data.getAddr<std::vector<double> >(variableName));
		} else if(variableType == "") {
			std::cerr<<"Could not find \"Variable\" \""<<variableName<<"\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
			return 0;
		} else {
			std::cerr<<"\"Variable\"'s type \""<<variableType<<"\" not supported by \"Plot\" (in \""<<plotName<<"\")."<<std::endl;
			return 0;
		}

	} else {

		std::vector<int> indices = __getIndices(plot, plotName);
		if(indices.empty()) {
			return 0;
		}
		std::stringstream strStr;
		strStr<<variableName<<indices[0];
		std::string variableType = data.getType(strStr.str());
		if(variableType == "double") {
			std::vector<double*>* vecData = __getDataVector<double>(plot, plotName, variableName, indices);
			if(not vecData) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double>(cutmasks, new TH1D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins, lowerBound, upperBound), vecData);
		} else if(variableType == "int") {
			std::vector<int*>* vecData = __getDataVector<int>(plot, plotName, variableName, indices);
			if(not vecData) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int>(cutmasks, new TH1D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins, lowerBound, upperBound), vecData);
		} else if(variableType == "std::vector<int>") {
			std::vector<std::vector<int>*>* vecData = __getDataVector<std::vector<int> >(plot, plotName, variableName, indices);
			if(not vecData) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int>(cutmasks, new TH1D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins, lowerBound, upperBound), vecData);
		} else if(variableType == "std::vector<double>") {
			std::vector<std::vector<double>*>* vecData = __getDataVector<std::vector<double> >(plot, plotName, variableName, indices);
			if(not vecData) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double>(cutmasks, new TH1D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins, lowerBound, upperBound), vecData);
		} else if(variableType == "") {
			std::cerr<<"Could not find \"Variable\" \""<<variableName<<"\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
			return 0;
		} else {
			std::cerr<<"\"Variable\"'s type \""<<variableType<<"\" not supported by \"Plot\" (in \""<<plotName<<"\")."<<std::endl;
			return 0;
		}

	}

	return antokPlot;

}

antok::Plot* antok::generators::generate2DPlot(const YAML::Node& plot, const antok::plotUtils::GlobalPlotOptions& plotOptions) {

	using antok::YAMLUtils::hasNodeKey;

	std::string plotName = antok::YAMLUtils::getString(plot["Name"]);
	std::string plotNameWithAxisLables = __getPlotNameWithAxisLabels(plot);
	if(plotNameWithAxisLables == "") {
		std::cerr<<"Could not get plot name with axis labels in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

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
	if(lowerBound1 >= upperBound1){
		std::cerr<<"\"LowerBoundX\" >= \"UpperBoundX\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	double upperBound2 = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["UpperBounds"][1], &upperBound2)) {
		std::cerr<<"Could not get second of the \"UpperBounds\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	if(lowerBound2 >= upperBound2){
		std::cerr<<"\"LowerBoundY\" >= \"UpperBoundY\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
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
	std::map<std::string, std::vector<long> > cutmasks;
	if(hasNodeKey(plot, "CustomCuts")) {
		if(not antok::Plotter::handleAdditionalCuts(plot["CustomCuts"], cutmasks)) {
			std::cerr<<"Warning: There was a problem when processing the \"CustomCuts\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		}
	}
	__getCutmasks(plotOptions, plot, cutmasks);

	antok::Plot* antokPlot = 0;

	if(not hasNodeKey(plot, "Indices")) {

		std::string variableType = data.getType(variable1Name);
		std::string variable2Type = data.getType(variable2Name);

		if(variableType == "double" && variable2Type == "double") {
			antokPlot = new antok::TemplatePlot<double,double>(cutmasks,
			                                                          new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                          data.getAddr<double>(variable1Name),
			                                                          data.getAddr<double>(variable2Name));
		} else if (variableType == "int" && variable2Type == "int") {
			antokPlot = new antok::TemplatePlot<int,int>(cutmasks,
			                                                 new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                 data.getAddr<int>(variable1Name),
			                                                 data.getAddr<int>(variable2Name));
		} else if (variableType == "std::vector<int>" && variable2Type == "std::vector<int>") {
			antokPlot = new antok::TemplatePlot<int,int>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             data.getAddr<std::vector<int> >(variable1Name),
			                                             data.getAddr<std::vector<int> >(variable2Name));
		} else if (variableType == "std::vector<double>" && variable2Type == "std::vector<double>") {
			antokPlot = new antok::TemplatePlot<double,double>(cutmasks,
			                                                          new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                          data.getAddr<std::vector<double> >(variable1Name),
			                                                          data.getAddr<std::vector<double> >(variable2Name));
		} else if (variableType == "double" && variable2Type == "int") {
			antokPlot = new antok::TemplatePlot<double,int>(cutmasks,
			                                                    new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                    data.getAddr<double >(variable1Name),
			                                                    data.getAddr<int>(variable2Name));
		} else if (variableType == "int" && variable2Type == "double") {
			antokPlot = new antok::TemplatePlot<int,double>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                data.getAddr<int >(variable1Name),
			                                                data.getAddr<double>(variable2Name));
		} else if (variableType == "int" && variable2Type == "std::vector<int>") {
			antokPlot = new antok::TemplatePlot<int,int>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                data.getAddr<int>(variable1Name),
			                                                data.getAddr<std::vector<int>>(variable2Name));
		} else if (variableType == "int" && variable2Type == "std::vector<double>") {
			antokPlot = new antok::TemplatePlot<int,double>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             data.getAddr<int>(variable1Name),
			                                             data.getAddr<std::vector<double>>(variable2Name));
		} else if (variableType == "std::vector<int>" && variable2Type == "int") {
			antokPlot = new antok::TemplatePlot<int,int>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             data.getAddr<std::vector<int>>(variable1Name),
			                                             data.getAddr<int>(variable2Name));
		} else if (variableType == "std::vector<double>" && variable2Type == "int") {
			antokPlot = new antok::TemplatePlot<double,int>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                data.getAddr<std::vector<double>>(variable1Name),
			                                                data.getAddr<int>(variable2Name));
		} else if (variableType == "double" && variable2Type == "std::vector<int>") {
			antokPlot = new antok::TemplatePlot<int,int>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             data.getAddr<int>(variable1Name),
			                                             data.getAddr<std::vector<int>>(variable2Name));
		} else if (variableType == "double" && variable2Type == "std::vector<double>") {
			antokPlot = new antok::TemplatePlot<double,double>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                data.getAddr<double>(variable1Name),
			                                                data.getAddr<std::vector<double>>(variable2Name));
		} else if (variableType == "std::vector<int>" && variable2Type == "double") {
			antokPlot = new antok::TemplatePlot<int,double>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             data.getAddr<std::vector<int>>(variable1Name),
			                                             data.getAddr<double>(variable2Name));
		} else if (variableType == "std::vector<double>" && variable2Type == "double") {
			antokPlot = new antok::TemplatePlot<double,double>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                data.getAddr<std::vector<double>>(variable1Name),
			                                                data.getAddr<double>(variable2Name));
		} else if (variableType == "std::vector<double>" && variable2Type == "std::vector<int>") {
			antokPlot = new antok::TemplatePlot<double,int>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                data.getAddr<std::vector<double>>(variable1Name),
			                                                data.getAddr<std::vector<int>>(variable2Name));
		} else if (variableType == "std::vector<int>" && variable2Type == "std::vector<double>") {
			antokPlot = new antok::TemplatePlot<int,double>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                data.getAddr<std::vector<int>>(variable1Name),
			                                                data.getAddr<std::vector<double>>(variable2Name));
		} else if(variableType == "") {
			std::cerr<<"Could not find \"Variable\" \""<<variable1Name<<"\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
			return 0;
		} else if(variable2Type == "") {
			std::cerr<<"Could not find \"Variable\" \""<<variable2Name<<"\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
			return 0;
		} else {
			std::cerr<<"\"Variable\" types \""<<variableType<<"\" and \""<<variable2Type<<"\" not supported by \"Plot\" (in \""<<plotName<<"\")."<<std::endl;
			return 0;
		}

	} else {

		std::vector<int> indices = __getIndices(plot, plotName);
		if(indices.empty()) {
			return 0;
		}
		std::stringstream strStr;
		strStr<<variable1Name<<indices[0];
		std::string variableType = data.getType(strStr.str());
		strStr.str("");
		strStr<<variable2Name<<indices[0];
		std::string variable2Type =  data.getType(strStr.str());

		if(variableType == "double" && variable2Type == "double") {
			std::vector<double*>* vec1Data = __getDataVector<double>(plot, plotName, variable1Name, indices);
			std::vector<double*>* vec2Data = __getDataVector<double>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double,double>(cutmasks,
			                                                          new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                          vec1Data,
			                                                          vec2Data);
		} else if (variableType == "int" && variable2Type == "int") {
			std::vector<int*>* vec1Data = __getDataVector<int>(plot, plotName, variable1Name, indices);
			std::vector<int*>* vec2Data = __getDataVector<int>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int,int>(cutmasks,
			                                                 new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                 vec1Data,
			                                                 vec2Data);
		} else if (variableType == "std::vector<double>" && variable2Type == "std::vector<double>") {
			std::vector<std::vector<double>*>* vec1Data = __getDataVector<std::vector<double> >(plot, plotName, variable1Name, indices);
			std::vector<std::vector<double>*>* vec2Data = __getDataVector<std::vector<double> >(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double,double>(cutmasks,
			                                                          new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                          vec1Data,
			                                                          vec2Data);
		} else if (variableType == "std::vector<int>" && variable2Type == "std::vector<int>") {
			std::vector<std::vector<int>*>* vec1Data = __getDataVector<std::vector<int> >(plot, plotName, variable1Name, indices);
			std::vector<std::vector<int>*>* vec2Data = __getDataVector<std::vector<int> >(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int,int>(cutmasks,
			                                                   new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                   vec1Data,
			                                                   vec2Data);
		} else if (variableType == "std::vector<int>" && variable2Type == "std::vector<double>") {
			std::vector<std::vector<int>*>* vec1Data = __getDataVector<std::vector<int> >(plot, plotName, variable1Name, indices);
			std::vector<std::vector<double>*>* vec2Data = __getDataVector<std::vector<double> >(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int,double>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             vec1Data,
			                                             vec2Data);
		} else if (variableType == "std::vector<double>" && variable2Type == "std::vector<int>") {
			std::vector<std::vector<double>*>* vec1Data = __getDataVector<std::vector<double> >(plot, plotName, variable1Name, indices);
			std::vector<std::vector<int>*>* vec2Data = __getDataVector<std::vector<int> >(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double,int>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                vec1Data,
			                                                vec2Data);
		} else if (variableType == "double" && variable2Type == "int") {
			std::vector<double*>* vec1Data = __getDataVector<double>(plot, plotName, variable1Name, indices);
			std::vector<int*>* vec2Data = __getDataVector<int>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double,int>(cutmasks,
			                                                    new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                    vec1Data,
			                                                    vec2Data);
		} else if (variableType == "int" && variable2Type == "double") {
			std::vector<int*>* vec1Data = __getDataVector<int>(plot, plotName, variable1Name, indices);
			std::vector<double*>* vec2Data = __getDataVector<double>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int,double>(cutmasks,
			                                                    new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                    vec1Data,
			                                                    vec2Data);
		} else if (variableType == "int" && variable2Type == "std::vector<int>") {
			std::vector<int*>* vec1Data = __getDataVector<int>(plot, plotName, variable1Name, indices);
			std::vector<std::vector<int>*>* vec2Data = __getDataVector<std::vector<int>>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int,int>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                vec1Data,
			                                                vec2Data);
		} else if (variableType == "std::vector<int>" && variable2Type == "int") {
			std::vector<std::vector<int>*>* vec1Data = __getDataVector<std::vector<int>>(plot, plotName, variable1Name, indices);
			std::vector<int*>* vec2Data = __getDataVector<int>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int,int>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             vec1Data,
			                                             vec2Data);
		} else if (variableType == "int" && variable2Type == "std::vector<double>") {
			std::vector<int*>* vec1Data = __getDataVector<int>(plot, plotName, variable1Name, indices);
			std::vector<std::vector<double>*>* vec2Data = __getDataVector<std::vector<double>>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int,double>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             vec1Data,
			                                             vec2Data);
		} else if (variableType == "std::vector<double>" && variable2Type == "int") {
			std::vector<std::vector<double>*>* vec1Data = __getDataVector<std::vector<double>>(plot, plotName, variable1Name, indices);
			std::vector<int*>* vec2Data = __getDataVector<int>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double,int>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             vec1Data,
			                                             vec2Data);
		} else if (variableType == "double" && variable2Type == "std::vector<int>") {
			std::vector<double*>* vec1Data = __getDataVector<double>(plot, plotName, variable1Name, indices);
			std::vector<std::vector<int>*>* vec2Data = __getDataVector<std::vector<int>>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double,int>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             vec1Data,
			                                             vec2Data);
		} else if (variableType == "std::vector<int>" && variable2Type == "double") {
			std::vector<std::vector<int>*>* vec1Data = __getDataVector<std::vector<int>>(plot, plotName, variable1Name, indices);
			std::vector<double*>* vec2Data = __getDataVector<double>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<int,double>(cutmasks,
			                                             new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                             vec1Data,
			                                             vec2Data);
		} else if (variableType == "double" && variable2Type == "std::vector<double>") {
			std::vector<double*>* vec1Data = __getDataVector<double>(plot, plotName, variable1Name, indices);
			std::vector<std::vector<double>*>* vec2Data = __getDataVector<std::vector<double>>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double,double>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                vec1Data,
			                                                vec2Data);
		} else if (variableType == "std::vector<double>" && variable2Type == "double") {
			std::vector<std::vector<double>*>* vec1Data = __getDataVector<std::vector<double>>(plot, plotName, variable1Name, indices);
			std::vector<double*>* vec2Data = __getDataVector<double>(plot, plotName, variable2Name, indices);
			if((not vec1Data) or (not vec2Data)) {
				return 0;
			}
			antokPlot = new antok::TemplatePlot<double,double>(cutmasks,
			                                                new TH2D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1, nBins2, lowerBound2, upperBound2),
			                                                vec1Data,
			                                                vec2Data);

		} else if(variableType == "") {
			std::cerr<<"Could not find \"Variable\" \""<<variable1Name<<"\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
			return 0;
		} else if(variable2Type == "") {
			std::cerr<<"Could not find \"Variable\" \""<<variable2Name<<"\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
			return 0;
		} else {
			std::cerr<<"\"Variable\" types \""<<variableType<<"\" and \""<<variable2Type<<"\" not supported by \"Plot\" (in \""<<plotName<<"\")."<<std::endl;
			return 0;
		}

	}

	return antokPlot;

}


antok::Plot* antok::generators::generate3DPlot(const YAML::Node& plot, const antok::plotUtils::GlobalPlotOptions& plotOptions) {

	using antok::YAMLUtils::hasNodeKey;

	std::string plotName = antok::YAMLUtils::getString(plot["Name"]);
	std::string plotNameWithAxisLables = __getPlotNameWithAxisLabels(plot);
	if(plotNameWithAxisLables == "") {
		std::cerr<<"Could not get plot name with axis labels in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

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
	double lowerBound3 = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["LowerBounds"][2], &lowerBound3)) {
		std::cerr<<"Could not get third of the \"LowerBounds\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

	double upperBound1 = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["UpperBounds"][0], &upperBound1)) {
		std::cerr<<"Could not get first of the \"UpperBounds\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	if(lowerBound1 >= upperBound1){
		std::cerr<<"\"LowerBoundX\" >= \"UpperBoundX\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	double upperBound2 = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["UpperBounds"][1], &upperBound2)) {
		std::cerr<<"Could not get second of the \"UpperBounds\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	double upperBound3 = 0.;
	if(not antok::YAMLUtils::getValue<double>(plot["UpperBounds"][2], &upperBound3)) {
		std::cerr<<"Could not get third of the \"UpperBounds\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	if(lowerBound2 >= upperBound2){
		std::cerr<<"\"LowerBoundY\" >= \"UpperBoundY\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	if(lowerBound3 >= upperBound3){
		std::cerr<<"\"LowerBoundZ\" >= \"UpperBoundZ\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

	int nBins1 = 0;
	if(not antok::YAMLUtils::getValue<int>(plot["NBins"][0], &nBins1)) {
		std::cerr<<"Could not get first of the \"NBins\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	int nBins2 = 0;
	if(not antok::YAMLUtils::getValue<int>(plot["NBins"][1], &nBins2)) {
		std::cerr<<"Could not get second of the \"NBins\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}
	int nBins3 = 0;
	if(not antok::YAMLUtils::getValue<int>(plot["NBins"][2], &nBins3)) {
		std::cerr<<"Could not get third of the \"NBins\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
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
	std::string variable3Name = antok::YAMLUtils::getString(plot["Variables"][2]);
	if(variable3Name == "") {
		std::cerr<<"Third of the \"Variables\" entries invalid for \"Plot\" \""<<plotName<<"\"."<<std::endl;
		return 0;
	}

	antok::Data& data = ObjectManager::instance()->getData();
	std::map<std::string, std::vector<long> > cutmasks;
	if(hasNodeKey(plot, "CustomCuts")) {
		if(not antok::Plotter::handleAdditionalCuts(plot["CustomCuts"], cutmasks)) {
			std::cerr<<"Warning: There was a problem when processing the \"CustomCuts\" in \"Plot\" \""<<plotName<<"\"."<<std::endl;
		}
	}
	__getCutmasks(plotOptions, plot, cutmasks);

	antok::Plot* antokPlot = 0;

	if(not hasNodeKey(plot, "Indices")) {

		std::string variableType = data.getType(variable1Name);
		std::string variable2Type = data.getType(variable2Name);
		std::string variable3Type = data.getType(variable3Name);

		if(variableType == "double" && variable2Type == "double" && variable3Type == "double") {
			antokPlot = new antok::TemplatePlot<double,double,double>(cutmasks,
			                                            new TH3D(plotName.c_str(), plotNameWithAxisLables.c_str(), nBins1, lowerBound1, upperBound1,
			                                                                                                       nBins2, lowerBound2, upperBound2,
			                                                                                                       nBins3, lowerBound3, upperBound3),
			                                            data.getAddr<double>(variable1Name),
			                                            data.getAddr<double>(variable2Name),
			                                            data.getAddr<double>(variable3Name));
		} else {
			std::cerr<<"\"Variable\" types \""<<variableType<<"\",\""<<variable2Type<<"\", and \""<<variable3Type<<"\" not supported by \"Plot\" (in \""<<plotName<<"\")."<<std::endl;
			return 0;
		}

	} else {
		std::cerr<<"\"Variables with indices not (yet) supported for 3D plots by \"Plot\" (in \""<<plotName<<"\")."<<std::endl;
		return 0;

	}


	return antokPlot;

}
