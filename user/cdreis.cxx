#include <fstream>

#include "cdreis.h"
#include "cdreis_functions.hpp"
#include "constants.h"
#include "data.h"
#include "functions.hpp"
#include "generators_functions.h"
#include "yaml_utils.hpp"

using antok::generators::nmbArgsIsExactly;
using antok::generators::functionArgumentHandler;
using antok::generators::functionArgumentHandlerConst;
using antok::generators::getFunctionArgumentHandlerErrorMsg;
using antok::YAMLUtils::hasNodeKey;


antok::Function*
antok::user::cdreis::getUserFunction(const YAML::Node&         function,
                                     std::vector<std::string>& quantityNames,
                                     int                       index)
{
	const std::string& functionName = antok::YAMLUtils::getString(function["Name"]);
	if        (functionName == "GetRecoilLorentzVec") {
		return antok::user::cdreis::generateGetRecoilLorentzVec             (function, quantityNames, index);
	} else if (functionName == "getPhotonLorentzVecs") {
		return antok::user::cdreis::generateGetPhotonLorentzVecs            (function, quantityNames, index);
	} else if (functionName == "getVectorLorentzVectorAttributes") {
		return antok::user::cdreis::generateGetVectorLorentzVectorAttributes(function, quantityNames, index);
	} else if (functionName == "getCleanedEcalClusters") {
		return antok::user::cdreis::generateGetCleanedEcalClusters          (function, quantityNames, index);
	} else if (functionName == "getPi0Pair") {
		return antok::user::cdreis::generateGetPi0Pair                      (function, quantityNames, index);
	} else if (functionName == "getOmega") {
		return antok::user::cdreis::generateGetOmega                        (function, quantityNames, index);
	} else if (functionName == "getECALCorrectedEnergy") {
		return antok::user::cdreis::generateGetECALCorrectedEnergy          (function, quantityNames, index);
	} else if (functionName == "getECALCorrectedTiming") {
		return antok::user::cdreis::generateGetECALCorrectedTiming          (function, quantityNames, index);
	} else if (functionName == "getPhotonPairParticles") {
		return antok::user::cdreis::generateGetPhotonPairParticles          (function, quantityNames, index);
	} else if (functionName == "getKinematicFittingMass") {
		return antok::user::cdreis::generateGetKinematicFittingMass         (function, quantityNames, index);
	} else if (functionName == "getThreePionCombinationMass") {
		return antok::user::cdreis::generateGetThreePionCombinationMass     (function, quantityNames, index);
	} else if (functionName == "getVector3VectorAttributes") {
		return antok::user::cdreis::generateGetVector3VectorAttributes      (function, quantityNames, index);
	} else if (functionName == "getECALVariables") {
		return antok::user::cdreis::generateGetECALVariables                (function, quantityNames, index);
	}
	return nullptr;
}


antok::Function*
antok::user::cdreis::generateGetRecoilLorentzVec(const YAML::Node&               function,
                                                 const std::vector<std::string>& quantityNames,
                                                 const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args
		= {{"BeamLorentzVec", "TLorentzVector"},
		   {"XLorentzVec"   , "TLorentzVector"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs = {{"RecoilMass", 0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<TLorentzVector>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetRecoilLorentzVec(*data.getAddr<TLorentzVector>(args[0].first),  // BeamLV
	                                                               *data.getAddr<TLorentzVector>(args[1].first),  // XLV
	                                                               constArgs["RecoilMass"],
	                                                               *data.getAddr<TLorentzVector>(quantityName));  // ResultRecoilLV
}


//TODO sync improved argument names of function with YAML strings
antok::Function*
antok::user::cdreis::generateGetPhotonLorentzVecs(const YAML::Node&               function,
                                                  const std::vector<std::string>& quantityNames,
                                                  const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 2)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args
		= {{"VectorPos",  "std::vector<TVector3>"},
		   {"VectorE",    "std::vector<double>"},
		   {"VectorTime", "std::vector<double>"},
		   {"PV",         "TVector3"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs = {{"RangeECAL1", 0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& ResultLVs         = quantityNames[0];
	const std::string& ResultECALIndices = quantityNames[1];
	antok::Data&       data              = antok::ObjectManager::instance()->getData();
	if (not data.insert<std::vector<TLorentzVector>>(ResultLVs)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames, ResultLVs);
		return nullptr;
	}
	if (not data.insert<std::vector<int>>(ResultECALIndices)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames, ResultECALIndices);
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetPhotonLorentzVecs(*data.getAddr<std::vector<TVector3>>(args[0].first),  // Positions
	                                                                *data.getAddr<std::vector<double>>(args[1].first),  // Energies
	                                                                *data.getAddr<std::vector<double>>(args[2].first),  // Times
	                                                                *data.getAddr<TVector3>             (args[3].first),  // PV
	                                                                constArgs["RangeECAL1"],                            // zECAL1
	                                                                *data.getAddr<std::vector<TLorentzVector>>(ResultLVs),
	                                                                *data.getAddr<std::vector<int>>           (ResultECALIndices));
}


antok::Function*
antok::user::cdreis::generateGetCleanedEcalClusters(const YAML::Node&               function,
                                                    const std::vector<std::string>& quantityNames,
                                                    const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 6)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args
		= {{"VectorX", "std::vector<double>"},
		   {"VectorY", "std::vector<double>"},
		   {"VectorZ", "std::vector<double>"},
		   {"VectorE", "std::vector<double>"},
		   {"VectorT", "std::vector<double>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs
		= {{"RangeECAL1",           0},
		   {"ThresholdEnergyECAL1", 0},
		   {"ThresholdTimingECAL1", 0},
		   {"ThresholdEnergyECAL2", 0},
		   {"ThresholdTimingECAL2", 0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	for (size_t i = 0; i < quantityNames.size() - 1; ++i) {
		if (not data.insert<std::vector<double>>(quantityNames[i])) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return nullptr;
		}
	}
	if (not data.insert<std::vector<int>>(quantityNames[5])) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[5]);
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetCleanedEcalClusters(*data.getAddr<std::vector<double>>(args[0].first),      // xPositions
	                                                                  *data.getAddr<std::vector<double>>(args[1].first),      // yPositions
	                                                                  *data.getAddr<std::vector<double>>(args[2].first),      // zPositions
	                                                                  *data.getAddr<std::vector<double>>(args[3].first),      // Energies
	                                                                  *data.getAddr<std::vector<double>>(args[4].first),      // Times
	                                                                  constArgs["RangeECAL1"],                                // zECAL1
	                                                                  constArgs["ThresholdEnergyECAL1"],                      // ThresholdEnergyECAL1
	                                                                  constArgs["ThresholdTimingECAL1"],                      // WindowTimingECAL1
	                                                                  constArgs["ThresholdEnergyECAL2"],                      // ThresholdEnergyECAL2
	                                                                  constArgs["ThresholdTimingECAL2"],                      // WindowTimingECAL2
	                                                                  *data.getAddr<std::vector<double>>(quantityNames[0]),   // ResultXPositions
	                                                                  *data.getAddr<std::vector<double>>(quantityNames[1]),   // ResultYPositions
	                                                                  *data.getAddr<std::vector<double>>(quantityNames[2]),   // ResultZPositions
	                                                                  *data.getAddr<std::vector<double>>(quantityNames[3]),   // ResultEnergies
	                                                                  *data.getAddr<std::vector<double>>(quantityNames[4]),   // ResultTimes
	                                                                  *data.getAddr<std::vector<int>>   (quantityNames[5]));  // ResultECALIndices
}


antok::Function*
antok::user::cdreis::generateGetVectorLorentzVectorAttributes(const YAML::Node&               function,
                                                              const std::vector<std::string>& quantityNames,
                                                              const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 7)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args = {{"VectorLV", "std::vector<TLorentzVector>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	for (size_t i = 0; i < quantityNames.size(); ++i) {
		if (not data.insert<std::vector<double>>(quantityNames[i])) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return nullptr;
		}
	}

	return new antok::user::cdreis::functions::GetVectorLorentzVectorAttributes(*data.getAddr<std::vector<TLorentzVector>>(args[0].first),  // LVs
	                                                                            *data.getAddr<std::vector<double>>(quantityNames[0]),       // ResultXComponents
	                                                                            *data.getAddr<std::vector<double>>(quantityNames[1]),       // ResultYComponents
	                                                                            *data.getAddr<std::vector<double>>(quantityNames[2]),       // ResultZComponents
	                                                                            *data.getAddr<std::vector<double>>(quantityNames[3]),       // ResultEnergies
	                                                                            *data.getAddr<std::vector<double>>(quantityNames[4]),       // ResultThetas
	                                                                            *data.getAddr<std::vector<double>>(quantityNames[5]),       // ResultPhis
	                                                                            *data.getAddr<std::vector<double>>(quantityNames[6]));      // ResultMags
}


antok::Function*
antok::user::cdreis::generateGetPi0Pair(const YAML::Node&               function,
                                        const std::vector<std::string>& quantityNames,
                                        const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 5)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args
		= {{"VectorLV" , "std::vector<TLorentzVector>"},
		   {"ECALIndex", "std::vector<int>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs
		= {{"Mass",            0},
		   {"ECALResolution",  0},
		   {"ECAL1Resolution", 0},
		   {"ECAL2Resolution", 0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& ResultPi0PairLVs  = quantityNames[0];
	const std::string& ResultPi0LV_0     = quantityNames[1];
	const std::string& ResultPi0LV_1     = quantityNames[2];
	const std::string& ResultGoodPi0Pair = quantityNames[3];
	const std::string& ResultECALIndices = quantityNames[4];
	antok::Data&       data              = antok::ObjectManager::instance()->getData();
	if (not data.insert<std::vector<TLorentzVector>>(ResultPi0PairLVs)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultPi0PairLVs);
		return nullptr;
	}
	if (not data.insert<TLorentzVector>(ResultPi0LV_0)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultPi0LV_0);
		return nullptr;
	}
	if (not data.insert<TLorentzVector>(ResultPi0LV_1)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultPi0LV_1);
		return nullptr;
	}
	if (not data.insert<int>(ResultGoodPi0Pair)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultGoodPi0Pair);
		return nullptr;
	}
	if (not data.insert<std::vector<int>>(ResultECALIndices)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultECALIndices);
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetPi0Pair(*data.getAddr<std::vector<TLorentzVector>>(args[0].first),  // PhotonLVs
	                                                      *data.getAddr<std::vector<int>>           (args[1].first),  // ECALIndices
	                                                      constArgs["Mass"],                                          // Pi0Mass
	                                                      constArgs["ECALResolution"],                                // MassWindowECALMixed
	                                                      constArgs["ECAL1Resolution"],                               // MassWindowECAL1
	                                                      constArgs["ECAL2Resolution"],                               // MassWindowECAL2
	                                                      *data.getAddr<std::vector<TLorentzVector>>(ResultPi0PairLVs),
	                                                      *data.getAddr<TLorentzVector>             (ResultPi0LV_0),
	                                                      *data.getAddr<TLorentzVector>             (ResultPi0LV_1),
	                                                      *data.getAddr<int>                        (ResultGoodPi0Pair),
	                                                      *data.getAddr<std::vector<int>>           (ResultECALIndices));
}


antok::Function*
antok::user::cdreis::generateGetOmega(const YAML::Node&               function,
                                      const std::vector<std::string>& quantityNames,
                                      const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 4)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args
		= {{"Pi0_0"     , "TLorentzVector"},
		   {"Pi0_1"     , "TLorentzVector"},
		   {"Scattered0", "TLorentzVector"},
		   {"Scattered1", "TLorentzVector"},
		   {"Scattered2", "TLorentzVector"},
		   {"Charge0"   , "int"},
		   {"Charge1"   , "int"},
		   {"Charge2"   , "int"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs
		= {{"Mass",            0},
		   {"ResolutionOmega", 0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	const std::string& ResultOmegaLV          = quantityNames[0];
	const std::string& ResultAccepted         = quantityNames[1];
	const std::string& ResultNotUsedPi0LV     = quantityNames[2];
	const std::string& ResultNotUsedPiMinusLV = quantityNames[3];
	antok::Data&       data                   = antok::ObjectManager::instance()->getData();
	if (not data.insert<TLorentzVector>(ResultOmegaLV)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultOmegaLV);
		return nullptr;
	}
	if (not data.insert<int>(ResultAccepted)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultAccepted);
		return nullptr;
	}
	if (not data.insert<TLorentzVector>(ResultNotUsedPi0LV)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultNotUsedPi0LV);
		return nullptr;
	}
	if (not data.insert<TLorentzVector>(ResultNotUsedPiMinusLV)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultNotUsedPiMinusLV);
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetOmega(*data.getAddr<TLorentzVector>(args[0].first),  // Pi0LV_0
	                                                    *data.getAddr<TLorentzVector>(args[1].first),  // Pi0LV_1
	                                                    *data.getAddr<TLorentzVector>(args[2].first),  // ChargedPartLV_0
	                                                    *data.getAddr<TLorentzVector>(args[3].first),  // ChargedPartLV_1
	                                                    *data.getAddr<TLorentzVector>(args[4].first),  // ChargedPartLV_2
	                                                    *data.getAddr<int>           (args[5].first),  // Charge_0
	                                                    *data.getAddr<int>           (args[6].first),  // Charge_1
	                                                    *data.getAddr<int>           (args[7].first),  // Charge_2
	                                                    constArgs["Mass"],                             // OmegaMass,
	                                                    constArgs["ResolutionOmega"],                  // MassWindowOmega,
	                                                    *data.getAddr<TLorentzVector>(ResultOmegaLV),
	                                                    *data.getAddr<int>           (ResultAccepted),
	                                                    *data.getAddr<TLorentzVector>(ResultNotUsedPi0LV),
	                                                    *data.getAddr<TLorentzVector>(ResultNotUsedPiMinusLV));
}


antok::Function*
antok::user::cdreis::generateGetECALCorrectedEnergy(const YAML::Node&               function,
                                                    const std::vector<std::string>& quantityNames,
                                                    const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args
		= {{"Energy"   , "std::vector<double>"},
		   {"ClusterZ" , "std::vector<double>"},
		   {"RunNumber", "int"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgsDouble = {{"RangeECAL1", 0}};
	if (not functionArgumentHandlerConst<double>(constArgsDouble, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	std::map<std::string, std::string> constArgsString = {{"Calibration", ""}};
	if (not functionArgumentHandlerConst<std::string>(constArgsString, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	// Read energy corrections from file
	const std::string& quantityName        = quantityNames[0];
	const std::string& CorrectionsFileName = constArgsString["Calibration"];
	std::map<int, std::pair<double, double>> Corrections;
	{
		std::ifstream CorrectionsFile;
		CorrectionsFile.open(CorrectionsFileName);
		if (not CorrectionsFile) {
			std::cerr << "Could not open file at '" << CorrectionsFileName << "' for 'Calibration' in function '" << function["Name"] << "' "
			          << "which is required for calculation of variable '" << quantityName << "'" << std::endl;
			return nullptr;
		}
		int    runNumber;
		double a, b;
		while (CorrectionsFile >> runNumber >> a >> b) {
			Corrections[runNumber] = std::pair<double, double>(a, b);
		}
		if (not CorrectionsFile.eof()) {
			std::cerr << "ERROR: Invalid ECAL energy-correction entries at end of file '" << CorrectionsFileName << "'; "
			          << "last good entry for run " << runNumber << "." << std::endl;
			return nullptr;
		}
		CorrectionsFile.close();
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	if (not data.insert<std::vector<double>>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetECALCorrectedEnergy(*data.getAddr<std::vector<double>>(args[0].first),  // Energies
	                                                                  *data.getAddr<std::vector<double>>(args[1].first),  // zPositions
	                                                                  constArgsDouble["RangeECAL1"],                       // zECAL1
	                                                                  *data.getAddr<int>                (args[2].first),  // RunNumber
	                                                                  Corrections,
	                                                                  *data.getAddr<std::vector<double>>(quantityName));  // ResultCorrectedEnergies
}


antok::Function*
antok::user::cdreis::generateGetECALCorrectedTiming(const YAML::Node&               function,
                                                    const std::vector<std::string>& quantityNames,
                                                    const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args
		= {{"Timing"  , "std::vector<double>"},
		   {"Energy"  , "std::vector<double>"},
		   {"ClusterZ", "std::vector<double>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgsDouble = {{"RangeECAL1", 0}};
	if (not functionArgumentHandlerConst<double>(constArgsDouble, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	std::map<std::string, std::string> constArgsString = {{"Calibration", ""}};
	if (not functionArgumentHandlerConst<std::string>(constArgsString, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	// Read time-correction coefficients from file
	const std::string& quantityName        = quantityNames[0];
	const std::string& CalibrationFileName = constArgsString["Calibration"];
	std::map<std::string, std::vector<double>> CalibrationCoeffs;
	{
		std::ifstream CalibrationFile;
		CalibrationFile.open(CalibrationFileName);
		if (not CalibrationFile) {
			std::cerr << "Could not open file at '" << CalibrationFileName << "' for 'Calibration' in function '" << function["Name"] << "' "
			          << "which is required for calculation of variable '" << quantityName << "'" << std::endl;
			return nullptr;
		}
		std::string ECALName;
		double      a, b, c, d, e, f, g;
		while (CalibrationFile >> ECALName >> a >> b >> c >> d >> e >> f >> g) {
			CalibrationCoeffs[ECALName] = {a, b, c, d, e, f, g};
		}
		if (not CalibrationFile.eof()) {
			std::cerr << "ERROR: Invalid ECAL time-calibration entries at end of file '" << CalibrationFileName << "'; "
			          << "last good entry for ECAL '" << ECALName << "'." << std::endl;
			return nullptr;
		}
		CalibrationFile.close();
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	if (not data.insert<std::vector<double>>(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetECALCorrectedTiming(*data.getAddr<std::vector<double>>(args[0].first),  // Times
	                                                                  *data.getAddr<std::vector<double>>(args[1].first),  // Energies
	                                                                  *data.getAddr<std::vector<double>>(args[2].first),  // zPositions
	                                                                  constArgsDouble["RangeECAL1"],                       // zECAL1
	                                                                  CalibrationCoeffs,
	                                                                  *data.getAddr<std::vector<double>>(quantityName));  // ResultCorrectedTimes
}


antok::Function*
antok::user::cdreis::generateGetPhotonPairParticles(const YAML::Node&               function,
                                                    const std::vector<std::string>& quantityNames,
                                                    const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 2)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args
		= {{"Photons"  , "std::vector<TLorentzVector>"},
		   {"ECALIndex", "std::vector<int>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs
		= {{"ECALResolution",  0},
		   {"ECAL1Resolution", 0},
		   {"ECAL2Resolution", 0},
		   {"Mass",            0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& ResultParticleLVs  = quantityNames[0];
	const std::string& ResultHasParticles = quantityNames[1];
	antok::Data&       data               = antok::ObjectManager::instance()->getData();
	if (not data.insert<std::vector<TLorentzVector>>(ResultParticleLVs)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultParticleLVs);
		return nullptr;
	}
	if (not data.insert<int>(ResultHasParticles)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultHasParticles);
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetPhotonPairParticles(*data.getAddr<std::vector<TLorentzVector>>(args[0].first),  // PhotonLVs
	                                                                  constArgs["Mass"],                                          // NominalMass,
	                                                                  constArgs["ECALResolution"],                                // MassWindowECALMixed,
	                                                                  constArgs["ECAL1Resolution"],                               // MassWindowECAL1,
	                                                                  constArgs["ECAL2Resolution"],                               // MassWindowECAL2,
	                                                                  *data.getAddr<std::vector<int>>           (args[1].first),  // ECALindices
	                                                                  *data.getAddr<std::vector<TLorentzVector>>(ResultParticleLVs),
	                                                                  *data.getAddr<int>                        (ResultHasParticles));
}


antok::Function*
antok::user::cdreis::generateGetKinematicFittingMass(const YAML::Node&               function,
                                                     const std::vector<std::string>& quantityNames,
                                                     const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 10)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args
		= {{"ClusterPositions",      "std::vector<TVector3>"},
		   {"ClusterPositionsError", "std::vector<TVector3>"},
		   {"VertexPosition",        "TVector3"},
		   {"ClusterEnergies",       "std::vector<double>"},
		   {"ClusterEnergiesError",  "std::vector<double>"},
		   {"ClusterIndices",        "std::vector<int>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, int> constArgsInt = {{"ErrorEstimateType", 0}};
	if (not functionArgumentHandlerConst<int>(constArgsInt, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	std::map<std::string, double> constArgsDouble
		= {{"PrecisionGoal", 0},
		   {"Mass",          0}};
	if (not functionArgumentHandlerConst<double>(constArgsDouble, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& ResultLorentzVectors = quantityNames[0];
	const std::string& ResultChi2s          = quantityNames[1];
	const std::string& ResultCLs            = quantityNames[2];
	const std::string& ResultSuccess        = quantityNames[3];
	antok::Data&       data                 = antok::ObjectManager::instance()->getData();
	if (not data.insert<std::vector<TLorentzVector>>(ResultLorentzVectors)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultLorentzVectors);
		return nullptr;
	}
	if (not data.insert<std::vector<double>>(ResultChi2s)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultChi2s);
		return nullptr;
	}
	for (size_t i = 4; i < 10; ++i) {
		if (not data.insert<std::vector<double>>(quantityNames[i])) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return nullptr;
		}
	}
	if (not data.insert<std::vector<double>>(ResultCLs)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultCLs);
		return nullptr;
	}
	if (not data.insert<int>(ResultSuccess)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(ResultSuccess);
		return nullptr;
	}

	//TODO check address mapping
	return new antok::user::cdreis::functions::GetKinematicFittingMass(*data.getAddr<TVector3>             (args[2].first),  // VertexPosition
	                                                                   *data.getAddr<std::vector<TVector3>>(args[0].first),  // ClusterPositions
	                                                                   *data.getAddr<std::vector<TVector3>>(args[1].first),  // ClusterPositionErrors
	                                                                   *data.getAddr<std::vector<double>>  (args[3].first),  // ClusterEnergies
	                                                                   *data.getAddr<std::vector<double>>  (args[4].first),  // ClusterEnergieErrors
	                                                                   *data.getAddr<std::vector<int>>     (args[5].first),  // ClusterIndices
	                                                                   constArgsDouble["Mass"],                              // Mass
	                                                                   constArgsDouble["PrecisionGoal"],                     // Window,
	                                                                   constArgsInt   ["ErrorEstimateType"],                 // EnergyErrorType,
	                                                                   *data.getAddr<std::vector<TLorentzVector>>(ResultLorentzVectors),
	                                                                   *data.getAddr<std::vector<double>>        (ResultChi2s),
	                                                                   *data.getAddr<std::vector<double>>        (quantityNames[4]),  // _ResultPullsX0
	                                                                   *data.getAddr<std::vector<double>>        (quantityNames[5]),  // _ResultPullsY0
	                                                                   *data.getAddr<std::vector<double>>        (quantityNames[6]),  // _ResultPullsE0
	                                                                   *data.getAddr<std::vector<double>>        (quantityNames[7]),  // _ResultPullsX1
	                                                                   *data.getAddr<std::vector<double>>        (quantityNames[8]),  // _ResultPullsY1
	                                                                   *data.getAddr<std::vector<double>>        (quantityNames[9]),  // _ResultPullsE1
	                                                                   *data.getAddr<std::vector<double>>        (ResultCLs),
	                                                                   *data.getAddr<int>                        (ResultSuccess));
}


antok::Function*
antok::user::cdreis::generateGetThreePionCombinationMass(const YAML::Node&               function,
                                                         const std::vector<std::string>& quantityNames,
                                                         const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string> > args
		= {{"Pi0_0"     , "TLorentzVector"},
		   {"Pi0_1"     , "TLorentzVector"},
		   {"Scattered0", "TLorentzVector"},
		   {"Scattered1", "TLorentzVector"},
		   {"Scattered2", "TLorentzVector"},
		   {"Charge0"   , "int"},
		   {"Charge1"   , "int"},
		   {"Charge2"   , "int"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, int> constArgs = {{"UseSquared", 0}};
	if (not functionArgumentHandlerConst<int>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	const std::string& quantityName = quantityNames[0];
	antok::Data&       data         = antok::ObjectManager::instance()->getData();
	if (not data.insert<std::vector<double> >(quantityName)) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetThreePionCombinationMass(*data.getAddr<TLorentzVector>(args[0].first),       // Pi0LV_0
	                                                                       *data.getAddr<TLorentzVector>(args[1].first),       // Pi0LV_1
	                                                                       *data.getAddr<TLorentzVector>(args[2].first),       // ChargedPartLV_0
	                                                                       *data.getAddr<TLorentzVector>(args[3].first),       // ChargedPartLV_1
	                                                                       *data.getAddr<TLorentzVector>(args[4].first),       // ChargedPartLV_2
	                                                                       *data.getAddr<int>           (args[5].first),       // Charge_0
	                                                                       *data.getAddr<int>           (args[6].first),       // Charge_1
	                                                                       *data.getAddr<int>           (args[7].first),       // Charge_2
	                                                                       constArgs["UseSquared"],                            // UseMassSquared,
	                                                                       *data.getAddr<std::vector<double>>(quantityName));  // Result
}

antok::Function*
antok::user::cdreis::generateGetVector3VectorAttributes(const YAML::Node&               function,
                                                              const std::vector<std::string>& quantityNames,
                                                              const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 3)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string>> args = {{"Vector3s", "std::vector<TVector3>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	for (size_t i = 0; i < quantityNames.size(); ++i) {
		if (not data.insert<std::vector<double>>(quantityNames[i])) {
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return nullptr;
		}
	}

	return new antok::user::cdreis::functions::GetVector3VectorAttributes(*data.getAddr<std::vector<TVector3>>(args[0].first),        // TVector3s
                                                                          *data.getAddr<std::vector<double>>(quantityNames[0]),       // ResultXComponents
                                                                          *data.getAddr<std::vector<double>>(quantityNames[1]),       // ResultYComponents
                                                                          *data.getAddr<std::vector<double>>(quantityNames[2]));      // ResultZComponents
}

antok::Function*
antok::user::cdreis::generateGetECALVariables(const YAML::Node&               function,
                                              const std::vector<std::string>& quantityNames,
                                              const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 6)) {
		return nullptr;
	}

	// Get input variables
	std::vector<std::pair<std::string, std::string> > args
		= {{"ECAL_clusterIndex"      , "std::vector<double>"},
		   {"PhotonVecsECAL_Vec"            , "std::vector<TLorentzVector>"},
		   {"ECAL_clusterPos"        , "std::vector<TVector3>"},
		   {"ECAL_clusterPosError"   , "std::vector<TVector3>"},
		   {"ECAL_clusterEnergy"     , "std::vector<double>"},
		   {"ECAL_clusterEnergyError", "std::vector<double>"},
		   {"ECAL_clusterT"     		, "std::vector<double>"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments

	std::map<std::string, int> constArgsInt = {{"ECALIndex", 0}};
		if (not functionArgumentHandlerConst<int>(constArgsInt, function)) {
			std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
			return nullptr;
		}

	// Register output variables
	antok::Data&       data         = antok::ObjectManager::instance()->getData();

	if (not data.insert<std::vector<TLorentzVector>>(quantityNames[0])) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[0]);
		return nullptr;
	}
	if (not data.insert<std::vector<TVector3>>(quantityNames[1])) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[1]);
		return nullptr;
	}
	if (not data.insert<std::vector<TVector3>>(quantityNames[2])) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[2]);
		return nullptr;
	}
	if (not data.insert<std::vector<double>>(quantityNames[3])) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[3]);
		return nullptr;
	}
	if (not data.insert<std::vector<double>>(quantityNames[4])) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[4]);
		return nullptr;
	}
	if (not data.insert<std::vector<double>>(quantityNames[5])) {
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[5]);
		return nullptr;
	}


	return new antok::user::cdreis::functions::getECALVariables( *data.getAddr<std::vector<double>>(args[0].first),            // ClusterIndex
																 *data.getAddr<std::vector<TLorentzVector>>(args[1].first),    // PhotonVec
															 	 *data.getAddr<std::vector<TVector3>>(args[2].first),          // ClusterPos
																 *data.getAddr<std::vector<TVector3>>(args[3].first),          // ClusterPosError
																 *data.getAddr<std::vector<double>>(args[4].first),            // ClusterE
																 *data.getAddr<std::vector<double>>(args[5].first),            // ClusterEError
																 *data.getAddr<std::vector<double>>(args[6].first),            // ClusterT
																 constArgsInt["ECALIndex"],                                    // SelectedECALIndex
																 *data.getAddr<std::vector<TLorentzVector>>(quantityNames[0]), // ResultPhotonVec
																 *data.getAddr<std::vector<TVector3>>(quantityNames[1]),       // ResultClusterPos
																 *data.getAddr<std::vector<TVector3>>(quantityNames[2]),       // ResultClusterPosError
																 *data.getAddr<std::vector<double>>(quantityNames[3]),         // ResultClusterE
																 *data.getAddr<std::vector<double>>(quantityNames[4]),         // ResultClusterEError
																 *data.getAddr<std::vector<double>>(quantityNames[5]));        // ResultClusterT
}
