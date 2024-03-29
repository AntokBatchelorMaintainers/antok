#include<generators_cuts.h>

#include<iostream>
#include<sstream>
#include<fstream>

#include<TLorentzVector.h>

#include<cut.hpp>
#include<data.h>
#include<initializer.h>
#include<object_manager.h>
#include<yaml_utils.hpp>

namespace {

	antok::Cut* __generateCut(const YAML::Node& cut,
	                          const std::string& shortName,
	                          const std::string& longName,
	                          const std::string& abbreviation,
	                          bool*& result);

	template<typename T>
	antok::Cut* __getEqualityCut(const YAML::Node& cut,
	                             const std::string& shortName,
	                             const std::string& longName,
	                             const std::string& abbreviation,
	                             bool* const result,
	                             antok::cuts::equalityMethod mode)
	{

		T* value = antok::YAMLUtils::getAddress<T>(cut["Value"]);
		if(value == nullptr) {
			std::cerr<<"Problem processing \"Value\" entry in \"Equality\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		T* variable = antok::YAMLUtils::getAddress<T>(cut["Variable"]);
		if(variable == nullptr) {
			std::cerr<<"Problem processing \"Variable\" entry in \"Equality\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		return (new antok::cuts::EqualityCut<T>(shortName, longName, abbreviation, result, variable, value, mode));

	}

	antok::Cut* __getEllipticCut(const YAML::Node& cut,
	                             const std::string& shortName,
	                             const std::string& longName,
	                             const std::string& abbreviation,
	                             bool* const result,
	                             antok::cuts::ellipticMethod mode)
	{


		double* meanX = antok::YAMLUtils::getAddress<double>(cut["meanX"]);
		if(meanX == nullptr) {
			std::cerr<<"Problem processing \"meanX\" entry in \"Elliptic\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		double* meanY = antok::YAMLUtils::getAddress<double>(cut["meanY"]);
		if(meanY == nullptr) {
			std::cerr<<"Problem processing \"meanY\" entry in \"Elliptic\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		double* cutX = antok::YAMLUtils::getAddress<double>(cut["cutX"]);
		if(cutX == nullptr) {
			std::cerr<<"Problem processing \"cutX\" entry in \"Elliptic\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		double* cutY = antok::YAMLUtils::getAddress<double>(cut["cutY"]);
		if(cutY == nullptr) {
			std::cerr<<"Problem processing \"cutY\" entry in \"Elliptic\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		double* X = antok::YAMLUtils::getAddress<double>(cut["X"]);
		if(X == nullptr) {
			std::cerr<<"Problem processing \"X\" entry in \"Elliptic\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		double* Y = antok::YAMLUtils::getAddress<double>(cut["Y"]);
		if(Y == nullptr) {
			std::cerr<<"Problem processing \"Y\" entry in \"Elliptic\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		return (new antok::cuts::EllipticCut(shortName, longName, abbreviation, result,
		                                        meanX, meanY, cutX, cutY, X, Y, mode));

	}


	template<typename T>
	antok::Cut* __getIsNotNANCut(const YAML::Node& cut,
	                             const std::string& shortName,
	                             const std::string& longName,
	                             const std::string& abbreviation,
	                             bool* const result)
	{

		T* variable = antok::YAMLUtils::getAddress<T>(cut["Variable"]);
		if(variable == nullptr) {
			std::cerr<<"Problem processing \"Variable\" entry in \"IsNotNAN\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		return (new antok::cuts::IsNotNANCut<T>(shortName, longName, abbreviation, result, variable));
	}

	template< typename T>
	antok::Cut* __generateRangeCutT(const YAML::Node& cut,
	                                const std::string& shortName,
	                                const std::string& longName,
	                                const std::string& abbreviation,
	                                bool* const result)
	{

		using antok::YAMLUtils::hasNodeKey;

		if(not (hasNodeKey(cut, "Type") and hasNodeKey(cut, "Variable"))) {
			std::cerr<<"A required entry is missing for \"Range\" cut \""<<shortName<<"\" (either \"Type\" or \"Variable\")"<<std::endl;
			return 0;
		}
		int method = antok::cuts::rangeFail;
		T* lowerBoundOrCentralValue = nullptr;
		T* upperBoundOrRange = nullptr;

		//=======================================
		// Range means central_value +- range
		//=======================================
		if( (hasNodeKey(cut, "CentralValue") or hasNodeKey(cut, "Range") ) and (hasNodeKey(cut, "LowerBound") or hasNodeKey(cut, "UpperBound") )) {
			std::cerr<<"Entries \"LowerBound\"/\"UpperBound\" invalid in combination with \"CentralValue\"/\"Range\" in \"Range\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		if( hasNodeKey(cut, "CentralValue") xor hasNodeKey(cut, "Range") ){
			std::cerr<<"Entries \"CentralValue\" and \"Range\" necessary in \"Range\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}




		if(hasNodeKey(cut, "LowerBound") and hasNodeKey(cut, "UpperBound")) {
			lowerBoundOrCentralValue = antok::YAMLUtils::getAddress<T>(cut["LowerBound"]);
			upperBoundOrRange = antok::YAMLUtils::getAddress<T>(cut["UpperBound"]);
			if(lowerBoundOrCentralValue == nullptr or upperBoundOrRange == nullptr) {
				std::cerr<<"Entries \"LowerBound\"/\"UpperBound\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type T."<<std::endl;
				return 0;
			}
			method = antok::cuts::rangeExcl;
		} else if (hasNodeKey(cut, "CentralValue") and hasNodeKey(cut, "Range") ){
			lowerBoundOrCentralValue = antok::YAMLUtils::getAddress<T>(cut["CentralValue"]);
			upperBoundOrRange = antok::YAMLUtils::getAddress<T>(cut["Range"]);
			if(lowerBoundOrCentralValue == nullptr or upperBoundOrRange == nullptr) {
				std::cerr<<"Entries \"CentralValue\"/\"Range\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type double."<<std::endl;
				return 0;
			}
			method = antok::cuts::centralValExcl;
		} else if(hasNodeKey(cut, "LowerBound")) {
			lowerBoundOrCentralValue = antok::YAMLUtils::getAddress<T>(cut["LowerBound"]);
			if(lowerBoundOrCentralValue == nullptr) {
				std::cerr<<"Entry \"LowerBound\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type T."<<std::endl;
				return 0;
			}
			method = antok::cuts::rangeOpenHighExcl;
		} else if(hasNodeKey(cut, "UpperBound")) {
			upperBoundOrRange = antok::YAMLUtils::getAddress<T>(cut["UpperBound"]);
			if(upperBoundOrRange == nullptr) {
				std::cerr<<"Entries \"LowerBound\"/\"UpperBound\" invalid in \"Range\" cut \""<<shortName<<"\", has to be either a variable name or of type double."<<std::endl;
				return 0;
			}
			method = antok::cuts::rangeOpenLowExcl;
		} else {
			std::cerr<<"Either \"LowerBound\" or \"UpperBound\" has to be present in \"Range\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		std::string type = antok::YAMLUtils::getString(cut["Type"]);
		if(type == "Inclusive") {
			method = method + 1; // change method to inclusive version
		} else if(type == "Exclusive") {
			// nothing to do
		} else if(type == "") {
			std::cerr<<"Could not convert \"Type\" to std::string for \"Range\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		} else {
			std::cerr<<"\"Type\" entry in \"Range\" cut \""<<shortName<<"\" has either to be \"Exclusive\" or \"Inclusive\""<<std::endl;
			return 0;
		}
		std::string varName = antok::YAMLUtils::getString(cut["Variable"]);
		if(varName == "") {
			std::cerr<<"Could not convert \"Range\" cut \""<<shortName<<"\"'s \"Variable\" entry to std::string."<<std::endl;
			return 0;
		}
		antok::Data& data = antok::ObjectManager::instance()->getData();

		T* variable = data.getAddr<T>(varName);
		if(variable == nullptr) {
			std::cerr<<"Could not find \"Range\" cut \""<<shortName<<"\"'s \"Variable\" entry \""<<varName<<"\" in Data."<<std::endl;
			return 0;
		}
		return (new antok::cuts::RangeCut<T>(shortName, longName, abbreviation, result, lowerBoundOrCentralValue, upperBoundOrRange, variable, (antok::cuts::rangeMethod) method));

	};
	antok::Cut* __generateRangeCut(const YAML::Node& cut,
	                                const std::string& shortName,
	                                const std::string& longName,
	                                const std::string& abbreviation,
	                                bool* const result)
	{

		std::string varName = antok::YAMLUtils::getString(cut["Variable"]);
		antok::Data& data = antok::ObjectManager::instance()->getData();
		const std::string varType = data.getType(varName);
		if( varType == "double" ) return __generateRangeCutT<double>( cut, shortName, longName, abbreviation, result );
		if( varType == "int" )    return __generateRangeCutT<int>   ( cut, shortName, longName, abbreviation, result );
		return nullptr;
	}

	antok::Cut* __generateEqualityCut(const YAML::Node& cut,
	                                  const std::string& shortName,
	                                  const std::string& longName,
	                                  const std::string& abbreviation,
	                                  bool* const result)
	{

		using antok::YAMLUtils::hasNodeKey;

		if(not (hasNodeKey(cut, "Type") and hasNodeKey(cut, "Value") and hasNodeKey(cut, "Variable"))) {
			std::cerr<<"One of the required arguments (\"Type\", \"Value\" and \"Variable\") for \"Equality\" cut \""<<shortName<<"\" is missing."<<std::endl;
			return 0;
		}

		std::string type = antok::YAMLUtils::getString(cut["Type"]);
		if(type == "") {
			std::cerr<<"Could not convert \"Type\" to std::string in \"Equality\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		antok::cuts::equalityMethod mode = antok::cuts::eqFail;
		if(type == "==") {
			mode = antok::cuts::equality;
		} else if (type == "!=") {
			mode = antok::cuts::inequality;
		} else {
			std::cerr<<"\"Type\" \""<<type<<"\" not supported by \"Equality\" cut."<<std::endl;
			return 0;
		}

		std::string variableName = antok::YAMLUtils::getString(cut["Variable"]);
		if(variableName == "") {
			std::cerr<<"Could not convert \"Equality\" cut \""<<shortName<<"\"'s \"Variable\" to std::string."<<std::endl;
			return 0;
		}

		antok::Data& data = antok::ObjectManager::instance()->getData();
		std::string typeName = data.getType(variableName);

		antok::Cut* antokCut = nullptr;
		if(typeName == "double") {
			antokCut = __getEqualityCut<double>(cut, shortName, longName, abbreviation, result, mode);
		} else if (typeName == "int") {
			antokCut = __getEqualityCut<int>(cut, shortName, longName, abbreviation, result, mode);
		} else if (typeName == "Long64_t") {
			antokCut = __getEqualityCut<Long64_t>(cut, shortName, longName, abbreviation, result, mode);
		} else if (typeName == "TLorentzVector") {
			antokCut = __getEqualityCut<TLorentzVector>(cut, shortName, longName, abbreviation, result, mode);
		} else {
			std::cerr<<"Variable type \""<<typeName<<"\" of variable \"" << variableName << "\" not supported in \"Equality\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		return antokCut;

	};

	antok::Cut* __generateEllipticCut(const YAML::Node& cut,
	                                  const std::string& shortName,
	                                  const std::string& longName,
	                                  const std::string& abbreviation,
	                                  bool* const result)
	{

		using antok::YAMLUtils::hasNodeKey;

		if(not (hasNodeKey(cut, "Type") and hasNodeKey(cut, "meanX") and hasNodeKey(cut, "meanY")
		   and hasNodeKey(cut, "cutX") and hasNodeKey(cut, "cutY") and hasNodeKey(cut, "X")
		   and hasNodeKey(cut, "Y"))) {
			std::cerr<<"One of the required arguments (\"Type\", \"meanX\", \"meanY\", \"cutX\", \"cutY\", \"X\", and \"Y\") for \"Equality\" cut \""<<shortName<<"\" is missing."<<std::endl;
			return 0;
		}

		std::string type = antok::YAMLUtils::getString(cut["Type"]);
		if(type == "") {
			std::cerr<<"Could not convert \"Type\" to std::string in \"Elliptic\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		antok::cuts::ellipticMethod mode = antok::cuts::ellipticFail;
		if(type == "Inclusive") {
			mode = antok::cuts::ellipticInclusive;
		} else if (type == "Exclusive") {
			mode = antok::cuts::ellipticExclusive;
		} else {
			std::cerr<<"\"Type\" \""<<type<<"\" not supported by \"Elliptic\" cut."<<std::endl;
			return 0;
		}

		std::string variableNameX = antok::YAMLUtils::getString(cut["X"]);
		if(variableNameX == "") {
			std::cerr<<"Could not convert \"Elliptic\" cut \""<<shortName<<"\"'s \"X\" to std::string."<<std::endl;
			return 0;
		}
		std::string variableNameY = antok::YAMLUtils::getString(cut["Y"]);
		if(variableNameY == "") {
			std::cerr<<"Could not convert \"Elliptic\" cut \""<<shortName<<"\"'s \"Y\" to std::string."<<std::endl;
			return 0;
		}

		antok::Data& data = antok::ObjectManager::instance()->getData();
		std::string typeName = data.getType(variableNameX);

		antok::Cut* antokCut = nullptr;
		if(typeName == "double") {
			antokCut = __getEllipticCut(cut, shortName, longName, abbreviation, result, mode);
		} else {
			std::cerr<<"Type \""<<typeName<<"\" not supported in \"Equality\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		return antokCut;

	};

	antok::Cut* __generateTriggerMaskCut(const YAML::Node& cut,
	                                     const std::string& shortName,
	                                     const std::string& longName,
	                                     const std::string& abbreviation,
	                                     bool* const result)
	{

		using antok::YAMLUtils::hasNodeKey;

		if(not (hasNodeKey(cut, "Type") and hasNodeKey(cut, "Mask") and hasNodeKey(cut, "Variable"))) {
			std::cerr<<"One of the required arguments (\"Type\", \"Mask\" and \"Variable\") for \"TriggerMask\" cut \""<<shortName<<"\" is missing."<<std::endl;
			return 0;
		}

		std::string type = antok::YAMLUtils::getString(cut["Type"]);
		if(type == "") {
			std::cerr<<"Could not convert \"Type\" to std::string in \"TriggerMask\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		int mode = 0;
		if(type != "Inclusive") {
			std::cerr<<"Only \"Inclusive\" supported as \"Type\" for TriggerMask cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		int* maskAddr = antok::YAMLUtils::getAddress<int>(cut["Mask"]);
		if(maskAddr == nullptr) {
			std::cerr<<"\"Mask\" entry invalid in \"TriggerMask\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		int* variable = antok::YAMLUtils::getAddress<int>(cut["Variable"]);
		if(variable == nullptr) {
			std::cerr<<"\"Variable\" entry invalid in \"TriggerMask\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		return (new antok::cuts::TriggerMaskCut(shortName, longName, abbreviation, result, maskAddr, variable, mode));

	};

	antok::Cut* __generateGroupCut(const YAML::Node& cut,
	                               const std::string& shortName,
	                               const std::string& longName,
	                               const std::string& abbreviation,
	                               bool* const result)
	{

		using antok::YAMLUtils::hasNodeKey;

		if(not hasNodeKey(cut, "Cuts")) {
			std::cerr<<"The required argument \"Cuts\" is missing for \"GroupCut\" \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		if(not hasNodeKey(cut, "Type")) {
			std::cerr<<"The required argument \"Type\" is missing for \"GroupCut\" \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		std::string type = antok::YAMLUtils::getString(cut["Type"]);
		antok::cuts::groupMethod mode = antok::cuts::groupFail;
		if(type == "And") {
			mode = antok::cuts::groupAnd;
		} else if (type == "Or") {
			mode = antok::cuts::groupOr;
		} else if (type == "Nand") {
			mode = antok::cuts::groupNand;
		} else if (type == "") {
			std::cerr<<"Could not convert \"GroupCut\" \""<<shortName<<"\"'s \"Type\" to std::string."<<std::endl;
			return 0;
		} else {
			std::cerr<<"\"Type\" in \"GroupCut\" \""<<shortName<<"\" has to be either \"And\" or \"Or\"."<<std::endl;
			return 0;
		}

		if(not cut["Cuts"].IsSequence()) {
			std::cerr<<"The \"Cuts\" entry for a \"GroupCut\" has to be a YAML sequence (in cut \""<<shortName<<"\")."<<std::endl;
			return 0;
		}

		std::vector<antok::Cut*> cuts;
		std::vector<bool*> results;

		unsigned int index = 0;
		for(YAML::const_iterator cuts_it = cut["Cuts"].begin(); cuts_it != cut["Cuts"].end(); ++cuts_it) {

			std::string innerShortName = shortName;
			std::string innerLongName = longName;
			std::string innerAbbreviation = abbreviation;

			std::stringstream strStr;
			strStr<<innerShortName<<index;
			innerShortName = strStr.str();
			strStr.str("");
			strStr<<innerLongName<<index;
			innerLongName = strStr.str();
			strStr.str("");
			strStr<<innerAbbreviation<<index;
			innerAbbreviation = strStr.str();

			const YAML::Node& cutEntry = (*cuts_it);

			if(not hasNodeKey(cutEntry, "Cut")) {
				std::cerr<<"Cut \""<<shortName<<"\" does not have required entry \"Cut\"."<<std::endl;
				return 0;
			}

			const YAML::Node& cut = cutEntry["Cut"];
			bool* innerResult = nullptr;

			antok::Cut* antokCut = __generateCut(cut, innerShortName, innerLongName, innerAbbreviation, innerResult);

			if(antokCut == nullptr) {
				std::cerr<<"Could not generate cut \""<<shortName<<"\" in \"Group\" cut \""<<shortName<<"\"."<<std::endl;
				return 0;
			}

			cuts.push_back(antokCut);
			results.push_back(innerResult);

			++index;

		}

		return (new antok::cuts::CutGroup(shortName, longName, abbreviation, result, cuts, results, mode));

	};

	antok::Cut* __generateListCut(const YAML::Node& cut,
	                              const std::string& shortName,
	                              const std::string& longName,
	                              const std::string& abbreviation,
	                              bool* const result)
	{

		using antok::YAMLUtils::hasNodeKey;
		antok::cuts::listMethod mode = antok::cuts::listFail;

		if(not (hasNodeKey(cut, "Type") and hasNodeKey(cut, "ListAddress") and hasNodeKey(cut, "Variable"))) {
			std::cerr<<"A required entry is missing for \"ListCut\" cut \""<<shortName<<"\" (either \"Type\" or \"ListAddress\" or \"Variable\")"<<std::endl;
			return 0;
		}

		std::string type = antok::YAMLUtils::getString(cut["Type"]);
		if(type == "Inclusive") {
			mode = antok::cuts::listInclusive;
		} else if(type == "Exclusive") {
			mode = antok::cuts::listExclusive;
		} else if(type == "") {
			std::cerr<<"Could not convert \"Type\" to std::string for \"ListCut\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		} else {
			std::cerr<<"\"Type\" entry in \"ListCut\" cut \""<<shortName<<"\" has either to be \"Exclusive\" or \"Inclusive\""<<std::endl;
			return 0;
		}

		std::string listAddress = antok::YAMLUtils::getString(cut["ListAddress"]);
		std::vector<int> list;
		{
			std::ifstream listFile;
			listFile.open(listAddress);
			if (not listFile) {
				std::cerr<<"Could not open file named \"" << listAddress << "\" for \"ListCut\" cut \""<<shortName<<"\"."<<std::endl;
				return 0;
			}
			int    variable;
			while (listFile >> variable) {
				list.push_back(variable);
			}
			if (not listFile.eof()) {
				std::cerr << "ERROR: Invalid run numbers at end of file '" << listAddress << "'; "
						<< "last good entry for run " << variable << "." << std::endl;
				return nullptr;
			}
			listFile.close();
		}

		std::string variableStr = antok::YAMLUtils::getString(cut["Variable"]);
		
		antok::Data& data = antok::ObjectManager::instance()->getData();
		std::string variableType = data.getType(variableStr);

		if(variableType != "int") {
			std::cerr<<"Variable type \""<<variableType<<"\" of variable \"Variable\" not supported in \"ListCut\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}
		int* variableAddr = data.getAddr<int>(variableStr);

		return (new antok::cuts::ListCut(shortName, longName, abbreviation, result, list, variableAddr, mode));

	};

	antok::Cut* __generateIsNotNANCut(const YAML::Node& cut,
	                                  const std::string& shortName,
	                                  const std::string& longName,
	                                  const std::string& abbreviation,
	                                  bool* const result)
	{

		using antok::YAMLUtils::hasNodeKey;

		if(not (hasNodeKey(cut, "Variable"))) {
			std::cerr<<"One of the required arguments (\"Variable\") for \"IsNotNAN\" cut \""<<shortName<<"\" is missing."<<std::endl;
			return 0;
		}

		std::string variableName = antok::YAMLUtils::getString(cut["Variable"]);
		if(variableName == "") {
			std::cerr<<"Could not convert \"IsNotNAN\" cut \""<<shortName<<"\"'s \"Variable\" to std::string."<<std::endl;
			return 0;
		}

		antok::Data& data = antok::ObjectManager::instance()->getData();
		std::string typeName = data.getType(variableName);

		antok::Cut* antokCut = nullptr;
		if(typeName == "double") {
			antokCut = __getIsNotNANCut<double>(cut, shortName, longName, abbreviation, result);
		} else if (typeName == "Long64_t") {
			antokCut = __getIsNotNANCut<Long64_t>(cut, shortName, longName, abbreviation, result);
		} else {
			std::cerr<<"Variable type \""<<typeName<<"\" of variable \"" << variableName << "\" not supported in \"IsNotNAN\" cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		return antokCut;

	};


	antok::Cut* __generateCut(const YAML::Node& cut,
	                          const std::string& shortName,
	                          const std::string& longName,
	                          const std::string& abbreviation,
	                          bool*& result)
	{
		std::string cutName = antok::YAMLUtils::getString(cut["Name"]);
		if(cutName == "") {
			std::cerr<<"Could not get the cut's \"Cut\"->\"Name\" for cut \""<<shortName<<"\"."<<std::endl;
			return 0;
		}

		result = new bool();
		antok::Cut* antokCut = nullptr;
		if(cutName == "Range") {
			antokCut = __generateRangeCut(cut, shortName, longName, abbreviation, result);
		} else if (cutName == "Equality") {
			antokCut = __generateEqualityCut(cut, shortName, longName, abbreviation, result);
		} else if (cutName == "Elliptic") {
			antokCut = __generateEllipticCut(cut, shortName, longName, abbreviation, result);
		} else if (cutName == "TriggerMask") {
			antokCut = __generateTriggerMaskCut(cut, shortName, longName, abbreviation, result);
		} else if (cutName == "Group") {
			antokCut = __generateGroupCut(cut, shortName, longName, abbreviation, result);
		} else if (cutName == "NoCut") {
			antokCut = new antok::cuts::NoCut(shortName, longName, abbreviation, result);
		} else if (cutName == "ListCut") {
			antokCut = __generateListCut(cut, shortName, longName, abbreviation, result);
		} else if (cutName == "IsNotNAN") {
			antokCut = __generateIsNotNANCut(cut, shortName, longName, abbreviation, result);
		} else {
			std::cerr<<"Cut \""<<cutName<<"\" not supported."<<std::endl;
			delete result;
			return 0;
		}
		if(antokCut == nullptr) {
			delete result;
			return 0;
		}
		return antokCut;
	}
}

bool antok::generators::generateCut(const YAML::Node& cutEntry, antok::Cut*& antokCut, bool*& result) {

	using antok::YAMLUtils::hasNodeKey;

	std::string shortName = antok::YAMLUtils::getString(cutEntry["ShortName"]);
	std::string longName = antok::YAMLUtils::getString(cutEntry["LongName"]);
	std::string abbreviation = antok::YAMLUtils::getString(cutEntry["Abbreviation"]);

	if(shortName == "" or longName == "" or abbreviation == "") {
		std::cerr<<"Did not find one of the cut's names (needed are \"ShortName\", \"LongName\" and \"Abbreviation\")."<<std::endl;
		return false;
	}

	if(not hasNodeKey(cutEntry, "Cut")) {
		std::cerr<<"Cut \""<<shortName<<"\" does not have required entry \"Cut\"."<<std::endl;
		return false;
	}

	const YAML::Node& cut = cutEntry["Cut"];

	antokCut = __generateCut(cut, shortName, longName, abbreviation, result);

	if(antokCut == nullptr) {
		return false;
	}

	return true;

};

