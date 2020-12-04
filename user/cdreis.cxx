#include <fstream>

#include "cdreis.h"
#include "cdreis_functions.hpp"
#include "constants.h"
#include "functions.hpp"
#include "generators_functions.h"
#include "yaml_utils.hpp"

using antok::generators::functionArgumentHandler;
using antok::generators::functionArgumentHandlerConst;
using antok::generators::getFunctionArgumentHandlerErrorMsg;
using antok::generators::nmbArgsIsExactly;
using antok::YAMLUtils::hasNodeKey;

// type alias to save some typing
template <class T>
using vecPairString = std::vector<std::pair<std::string, T>>;


antok::Function*
antok::user::cdreis::getUserFunction(const YAML::Node&               function,
                                     const std::vector<std::string>& quantityNames,
                                     int                             index)
{
	const std::string& functionName = antok::YAMLUtils::getString(function["Name"]);
	if        (functionName == "getVector3VectorAttributes") {
		return antok::user::cdreis::generateGetVector3VectorAttributes       (function, quantityNames, index);
	} else if (functionName == "getVectorLorentzVectorAttributes") {
		return antok::user::cdreis::generateGetVectorLorentzVectorAttributes (function, quantityNames, index);
	} else if (functionName == "getNominalMassDifferences") {
		return antok::user::cdreis::generateGetNominalMassDifferences        (function, quantityNames, index);
	} else if (functionName == "getRecoilLorentzVec") {
		return antok::user::cdreis::generateGetRecoilLorentzVec              (function, quantityNames, index);
	} else if (functionName == "getECALCorrectedEnergy") {
		return antok::user::cdreis::generateGetECALCorrectedEnergy           (function, quantityNames, index);
	} else if (functionName == "getECALCorrectedTiming") {
		return antok::user::cdreis::generateGetECALCorrectedTiming           (function, quantityNames, index);
	} else if (functionName == "getCleanedEcalClusters") {
		return antok::user::cdreis::generateGetCleanedEcalClusters           (function, quantityNames, index);
	} else if (functionName == "getECALVariables") {
		return antok::user::cdreis::generateGetECALVariables                 (function, quantityNames, index);
	} else if (functionName == "getPhotonLorentzVecs") {
		return antok::user::cdreis::generateGetPhotonLorentzVecs             (function, quantityNames, index);
	} else if (functionName == "getPhotonPairParticles") {
		return antok::user::cdreis::generateGetPhotonPairParticles           (function, quantityNames, index);
	} else if (functionName == "getPi0Pair") {
		return antok::user::cdreis::generateGetPi0Pair                       (function, quantityNames, index);
	} else if (functionName == "getKinematicFittingMass") {
		return antok::user::cdreis::generateGetKinematicFittingMass          (function, quantityNames, index);
	} else if (functionName == "getOmega") {
		return antok::user::cdreis::generateGetOmega                         (function, quantityNames, index);
	} else if (functionName == "getFittedOmegaMassVsPrecisionGoal") {
		return antok::user::cdreis::generateGetFittedOmegaMassVsPrecisionGoal(function, quantityNames, index);
	} else if (functionName == "getThreePionCombinationMass") {
		return antok::user::cdreis::generateGetThreePionCombinationMass      (function, quantityNames, index);
	}
	return nullptr;
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
	vecPairString<std::string> args = {{"Vector3s", "std::vector<TVector3>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes
		= {"std::vector<double>",   // ResultXComponents
		   "std::vector<double>",   // ResultYComponents
		   "std::vector<double>"};  // ResultYComponents
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetVector3VectorAttributes(
		*data.getAddr<std::vector<TVector3>>(args[0].first),   // TVector3s
		*data.getAddr<std::vector<double>>(quantityNames[0]),  // ResultXComponents
		*data.getAddr<std::vector<double>>(quantityNames[1]),  // ResultYComponents
		*data.getAddr<std::vector<double>>(quantityNames[2])   // ResultZComponents
	);
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
	vecPairString<std::string> args = {{"VectorLV", "std::vector<TLorentzVector>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes
		= {"std::vector<double>",   // ResultXComponents
		   "std::vector<double>",   // ResultYComponents
		   "std::vector<double>",   // ResultZComponents
		   "std::vector<double>",   // ResultEnergies
		   "std::vector<double>",   // ResultThetas
		   "std::vector<double>",   // ResultPhis
		   "std::vector<double>"};  // ResultMags
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetVectorLorentzVectorAttributes(
			*data.getAddr<std::vector<TLorentzVector>>(args[0].first),  // LVs
			*data.getAddr<std::vector<double>>(quantityNames[0]),       // ResultXComponents
			*data.getAddr<std::vector<double>>(quantityNames[1]),       // ResultYComponents
			*data.getAddr<std::vector<double>>(quantityNames[2]),       // ResultZComponents
			*data.getAddr<std::vector<double>>(quantityNames[3]),       // ResultEnergies
			*data.getAddr<std::vector<double>>(quantityNames[4]),       // ResultThetas
			*data.getAddr<std::vector<double>>(quantityNames[5]),       // ResultPhis
			*data.getAddr<std::vector<double>>(quantityNames[6])        // ResultMags
	);
}

antok::Function*
antok::user::cdreis::generateGetNominalMassDifferences(const YAML::Node&               function,
                                                       const std::vector<std::string>& quantityNames,
                                                       const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args = {{"VectorLV", "std::vector<TLorentzVector>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs = {{"NominalMass", 0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

    // Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes = {"std::vector<double>"};  // ResultMassDifferences
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::getNominalMassDifferences(
		*data.getAddr<std::vector<TLorentzVector>>(args[0].first),  // VectorLV
		constArgs["NominalMass"],                                   // NominalMass
		*data.getAddr<std::vector<double>>(quantityNames[0])        // ResultMassDifferences
	);
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
	vecPairString<std::string> args
		= {{"BeamLorentzVec", "TLorentzVector"},
		   {"XLorentzVec",    "TLorentzVector"}};
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
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes = {"TLorentzVector"};  // ResultRecoilLV
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetRecoilLorentzVec(
		*data.getAddr<TLorentzVector>(args[0].first),    // BeamLV
		*data.getAddr<TLorentzVector>(args[1].first),    // XLV
		constArgs["RecoilMass"],
		*data.getAddr<TLorentzVector>(quantityNames[0])  // ResultRecoilLV
	);
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
	vecPairString<std::string> args
		= {{"Energies",           "std::vector<double>"},
		   {"ECALClusterIndices", "std::vector<int>"},
		   {"RunNumber",          "int"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
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
	const std::vector<std::string> outputVarTypes = {"std::vector<double>"};  // ResultCorrectedEnergies
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetECALCorrectedEnergy(
		*data.getAddr<std::vector<double>>(args[0].first),    // Energies
		*data.getAddr<std::vector<int>>   (args[1].first),    // ClusterIndices
		*data.getAddr<int>                (args[2].first),    // RunNumber
		Corrections,
		*data.getAddr<std::vector<double>>(quantityNames[0])  // ResultCorrectedEnergies
	);
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
	vecPairString<std::string> args
		= {{"Times",              "std::vector<double>"},
		   {"Energies",           "std::vector<double>"},
		   {"ECALClusterIndices", "std::vector<int>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
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
	const std::vector<std::string> outputVarTypes = {"std::vector<double>"};  // ResultCorrectedTimes
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetECALCorrectedTiming(
		*data.getAddr<std::vector<double>>(args[0].first),    // Times
		*data.getAddr<std::vector<double>>(args[1].first),    // Energies
		*data.getAddr<std::vector<int>>   (args[2].first),    // ClusterIndices
		CalibrationCoeffs,
		*data.getAddr<std::vector<double>>(quantityNames[0])  // ResultCorrectedTimes
	);
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
	vecPairString<std::string> args
		= {{"Positions",          "std::vector<TVector3>"},
		   {"PositionVariances",  "std::vector<TVector3>"},
		   {"Energies",           "std::vector<double>"},
		   {"EnergyVariances",    "std::vector<double>"},
		   {"Times",              "std::vector<double>"},
		   {"ECALClusterIndices", "std::vector<int>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs
		= {{"ECAL1ThresholdEnergy", 0},
		   {"ECAL1ThresholdTiming", 0},
		   {"ECAL2ThresholdEnergy", 0},
		   {"ECAL2ThresholdTiming", 0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes
		 = {"std::vector<TVector3>",  // ResultPositions
		    "std::vector<TVector3>",  // ResultPositionVariances
		    "std::vector<double>",    // ResultEnergies
		    "std::vector<double>",    // ResultEnergyVariances
		    "std::vector<double>",    // ResultTimes
		    "std::vector<int>"};      // ResultECALIndices
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetCleanedEcalClusters(
		*data.getAddr<std::vector<TVector3>>(args[0].first),     // Positions
		*data.getAddr<std::vector<TVector3>>(args[1].first),     // PositionVariances
		*data.getAddr<std::vector<double>>  (args[2].first),     // Energies
		*data.getAddr<std::vector<double>>  (args[3].first),     // EnergyVariances
		*data.getAddr<std::vector<double>>  (args[4].first),     // Times
		*data.getAddr<std::vector<int>>     (args[5].first),     // ECAL Cluster Indices
		constArgs["ECAL1ThresholdEnergy"],                       // ECAL1ThresholdEnergy
		constArgs["ECAL1ThresholdTiming"],                       // ECAL1ThresholdTiming
		constArgs["ECAL2ThresholdEnergy"],                       // ECAL2ThresholdEnergy
		constArgs["ECAL2ThresholdTiming"],                       // ECAL2ThresholdTiming
		*data.getAddr<std::vector<TVector3>>(quantityNames[0]),  // ResultPositions
		*data.getAddr<std::vector<TVector3>>(quantityNames[1]),  // ResultPositionVariances
		*data.getAddr<std::vector<double>>  (quantityNames[2]),  // ResultEnergies
		*data.getAddr<std::vector<double>>  (quantityNames[3]),  // ResultEnergyVariances
		*data.getAddr<std::vector<double>>  (quantityNames[4]),  // ResultTimes
		*data.getAddr<std::vector<int>>     (quantityNames[5])   // ResultECALClusterIndices
	);
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
	vecPairString<std::string> args
		= {{"ECALClusterIndices", "std::vector<int>"},
		   {"Positions",          "std::vector<TVector3>"},
		   {"PositionVariances",  "std::vector<TVector3>"},
		   {"Energies",           "std::vector<double>"},
		   {"EnergyVariances",    "std::vector<double>"},
		   {"Times",              "std::vector<double>"}};
	if (not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, int> constArgsInt = {{"SelectedECALIndex", 0}};
		if (not functionArgumentHandlerConst<int>(constArgsInt, function)) {
			std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
			return nullptr;
		}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes
		= {"std::vector<int>",       // ResultECALClusterIndex
		   "std::vector<TVector3>",  // ResultClusterPosition
		   "std::vector<TVector3>",  // ResultClusterPositionVariance
		   "std::vector<double>",    // ResultClusterEnergy
		   "std::vector<double>",    // ResultClusterEnergyVariances
		   "std::vector<double>"};   // ResultClusterTime
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::getECALVariables(
		*data.getAddr<std::vector<int>>     (args[0].first),     // ECALClusterIndices
		*data.getAddr<std::vector<TVector3>>(args[1].first),     // Positions
		*data.getAddr<std::vector<TVector3>>(args[2].first),     // PositionVariances
		*data.getAddr<std::vector<double>>  (args[3].first),     // Energies
		*data.getAddr<std::vector<double>>  (args[4].first),     // EnergyVariances
		*data.getAddr<std::vector<double>>  (args[5].first),     // Times
		constArgsInt["SelectedECALIndex"],                       // SelectedECALIndex
		*data.getAddr<std::vector<int>>     (quantityNames[0]),  // ResultECALClusterIndices
		*data.getAddr<std::vector<TVector3>>(quantityNames[1]),  // ResultPositions
		*data.getAddr<std::vector<TVector3>>(quantityNames[2]),  // ResultPositionVariances
		*data.getAddr<std::vector<double>>  (quantityNames[3]),  // ResultEnergies
		*data.getAddr<std::vector<double>>  (quantityNames[4]),  // ResultEnergyVariances
		*data.getAddr<std::vector<double>>  (quantityNames[5])   // ResultTimes
	);
}


antok::Function*
antok::user::cdreis::generateGetPhotonLorentzVecs(const YAML::Node&               function,
                                                  const std::vector<std::string>& quantityNames,
                                                  const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 1)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args
		= {{"Positions",          "std::vector<TVector3>"},
		   {"Energies",           "std::vector<double>"},
		   {"Times",              "std::vector<double>"},
		   {"ECALClusterIndices", "std::vector<int>"},
		   {"PV",                 "TVector3"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes = {"std::vector<TLorentzVector>"};  // ResultPhotonLVs
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetPhotonLorentzVecs(
		*data.getAddr<std::vector<TVector3>>(args[0].first),          // Positions
		*data.getAddr<std::vector<double>>  (args[1].first),          // Energies
		*data.getAddr<std::vector<double>>  (args[2].first),          // Times
		*data.getAddr<std::vector<int>>     (args[3].first),          // ECALClusterIndices
		*data.getAddr<TVector3>             (args[4].first),          // VertexPosition
		*data.getAddr<std::vector<TLorentzVector>>(quantityNames[0])  // ResultPhotonLVs
	);
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
	vecPairString<std::string> args
		= {{"PhotonLVs",          "std::vector<TLorentzVector>"},
		   {"ECALClusterIndices", "std::vector<int>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs
		= {{"NominalMass",         0},
		   {"ECALMixedMassWindow", 0},
		   {"ECAL1MassWindow",     0},
		   {"ECAL2MassWindow",     0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes
		= {"std::vector<TLorentzVector>",  // ResultParticleLVs
		   "int"};                         // ResultHasParticles
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetPhotonPairParticles(
		*data.getAddr<std::vector<TLorentzVector>>(args[0].first),     // PhotonLVs
		*data.getAddr<std::vector<int>>           (args[1].first),     // ECALClusterIndices
		constArgs["NominalMass"],                                      // NominalMass
		constArgs["ECALMixedMassWindow"],                              // ECALMixedMassWindow
		constArgs["ECAL1MassWindow"],                                  // ECAL1MassWindow
		constArgs["ECAL2MassWindow"],                                  // ECAL2MassWindow
		*data.getAddr<std::vector<TLorentzVector>>(quantityNames[0]),  // ResultParticleLVs
		*data.getAddr<int>                        (quantityNames[1])   // ResultHasParticles
	);
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
	vecPairString<std::string> args
		= {{"PhotonLVs",          "std::vector<TLorentzVector>"},
		   {"ECALClusterIndices", "std::vector<int>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs
		= {{"Pi0Mass",             0},
		   {"ECALMixedMassWindow", 0},
		   {"ECAL1MassWindow",     0},
		   {"ECAL2MassWindow",     0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes
		= {"std::vector<TLorentzVector>",  // ResultPi0PairLVs
		   "TLorentzVector",               // ResultPi0LV_0
		   "TLorentzVector",               // ResultPi0LV_1
		   "int",                          // ResultGoodPi0Pair
		   "std::vector<int>"};            // ResultECALClusterIndices
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetPi0Pair(
		*data.getAddr<std::vector<TLorentzVector>>(args[0].first),     // PhotonLVs
		*data.getAddr<std::vector<int>>           (args[1].first),     // ECALClusterIndices
		constArgs["Pi0Mass"],                                          // Pi0Mass
		constArgs["ECALMixedMassWindow"],                              // ECALMixedMassWindow
		constArgs["ECAL1MassWindow"],                                  // ECAL1MassWindow
		constArgs["ECAL2MassWindow"],                                  // ECAL2MassWindow
		*data.getAddr<std::vector<TLorentzVector>>(quantityNames[0]),  // ResultPi0PairLVs
		*data.getAddr<TLorentzVector>             (quantityNames[1]),  // ResultPi0LV_0
		*data.getAddr<TLorentzVector>             (quantityNames[2]),  // ResultPi0LV_1
		*data.getAddr<int>                        (quantityNames[3]),  // ResultGoodPi0Pair
		*data.getAddr<std::vector<int>>           (quantityNames[4])   // ResultECALClusterIndices
	);
}


antok::Function*
antok::user::cdreis::generateGetKinematicFittingMass(const YAML::Node&               function,
                                                     const std::vector<std::string>& quantityNames,
                                                     const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 11)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args
		= {{"VertexPosition",           "TVector3"},
		   {"ClusterPositions",         "std::vector<TVector3>"},
		   {"ClusterPositionsVariance", "std::vector<TVector3>"},
		   {"ClusterEnergies",          "std::vector<double>"},
		   {"ClusterEnergiesVariance",  "std::vector<double>"},
		   {"ClusterIndices",           "std::vector<int>"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, int> constArgsInt = {{"WhichEnergyVariance", 0}};
	if (not functionArgumentHandlerConst<int>(constArgsInt, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}
	std::map<std::string, double> constArgsDouble
		= {{"Mass",           0},
		   {"MassLowerLimit", 0},
		   {"MassUpperLimit", 0},
		   {"PrecisionGoal",  0}};
	if (not functionArgumentHandlerConst<double>(constArgsDouble, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes
		= {"std::vector<TLorentzVector>",  // ResultLorentzVectors
		   "std::vector<double>",          // ResultChi2s
		   "std::vector<double>",          // ResultPValues
		   "std::vector<int>",             // ResultNmbIterations
		   "int",                          // ResultSuccess
		   "std::vector<double>",          // ResultPullsX0
		   "std::vector<double>",          // ResultPullsY0
		   "std::vector<double>",          // ResultPullsE0
		   "std::vector<double>",          // ResultPullsX1
		   "std::vector<double>",          // ResultPullsY1
		   "std::vector<double>"};         // ResultPullsE1
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}


	return new antok::user::cdreis::functions::GetKinematicFittingMass(
		*data.getAddr<TVector3>             (args[0].first),           // VertexPosition
		*data.getAddr<std::vector<TVector3>>(args[1].first),           // ClusterPositions
		*data.getAddr<std::vector<TVector3>>(args[2].first),           // ClusterPositionVariances
		*data.getAddr<std::vector<double>>  (args[3].first),           // ClusterEnergies
		*data.getAddr<std::vector<double>>  (args[4].first),           // ClusterEnergieVariances
		*data.getAddr<std::vector<int>>     (args[5].first),           // ClusterIndices
		constArgsDouble["Mass"],                                       // Mass
		constArgsDouble["MassLowerLimit"],                             // MassLowerLimit,
		constArgsDouble["MassUpperLimit"],                             // MassUpperLimit,
		constArgsDouble["PrecisionGoal"],                              // PrecisionGoal
		constArgsInt   ["WhichEnergyVariance"],                        // WhichEnergyVariance,
		*data.getAddr<std::vector<TLorentzVector>>(quantityNames[0]),  // ResultLorentzVectors
		*data.getAddr<std::vector<double>>        (quantityNames[1]),  // ResultChi2s
		*data.getAddr<std::vector<double>>        (quantityNames[2]),  // ResultPValues
		*data.getAddr<std::vector<int>>           (quantityNames[3]),  // ResultNmbIterations
		*data.getAddr<int>                        (quantityNames[4]),  // ResultSuccess
		*data.getAddr<std::vector<double>>        (quantityNames[5]),  // ResultPullsX0
		*data.getAddr<std::vector<double>>        (quantityNames[6]),  // ResultPullsY0
		*data.getAddr<std::vector<double>>        (quantityNames[7]),  // ResultPullsE0
		*data.getAddr<std::vector<double>>        (quantityNames[8]),  // ResultPullsX1
		*data.getAddr<std::vector<double>>        (quantityNames[9]),  // ResultPullsY1
		*data.getAddr<std::vector<double>>        (quantityNames[10])  // ResultPullsE1
	);
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
	vecPairString<std::string> args
		= {{"Pi0LV_0",         "TLorentzVector"},
		   {"Pi0LV_1",         "TLorentzVector"},
		   {"ChargedPartLV_0", "TLorentzVector"},
		   {"ChargedPartLV_1", "TLorentzVector"},
		   {"ChargedPartLV_2", "TLorentzVector"},
		   {"Charge_0",        "int"},
		   {"Charge_1",        "int"},
		   {"Charge_2",        "int"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, double> constArgs
		= {{"OmegaMass",       0},
		   {"OmegaMassWindow", 0}};
	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes
		= {"TLorentzVector",   // ResultOmegaLV
		   "int",              // ResultAccepted
		   "TLorentzVector",   // ResultNotUsedPi0LV
		   "TLorentzVector"};  // ResultNotUsedPiMinusLV
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetOmega(
		*data.getAddr<TLorentzVector>(args[0].first),     // Pi0LV_0
		*data.getAddr<TLorentzVector>(args[1].first),     // Pi0LV_1
		*data.getAddr<TLorentzVector>(args[2].first),     // ChargedPartLV_0
		*data.getAddr<TLorentzVector>(args[3].first),     // ChargedPartLV_1
		*data.getAddr<TLorentzVector>(args[4].first),     // ChargedPartLV_2
		*data.getAddr<int>           (args[5].first),     // Charge_0
		*data.getAddr<int>           (args[6].first),     // Charge_1
		*data.getAddr<int>           (args[7].first),     // Charge_2
		constArgs["OmegaMass"],                           // OmegaMass,
		constArgs["OmegaMassWindow"],                     // MassWindowOmega,
		*data.getAddr<TLorentzVector>(quantityNames[0]),  // ResultOmegaLV
		*data.getAddr<int>           (quantityNames[1]),  // ResultAccepted
		*data.getAddr<TLorentzVector>(quantityNames[2]),  // ResultNotUsedPi0LV
		*data.getAddr<TLorentzVector>(quantityNames[3])   // ResultNotUsedPiMinusLV
	);
}


antok::Function*
antok::user::cdreis::generateGetFittedOmegaMassVsPrecisionGoal(const YAML::Node&               function,
                                                               const std::vector<std::string>& quantityNames,
                                                               const int                       index)
{
	if (not nmbArgsIsExactly(function, quantityNames.size(), 3)) {
		return nullptr;
	}

	// Get input variables
	vecPairString<std::string> args
		= {{"VertexPosition",           "TVector3"},
		   {"ChargedPartLV_0",          "TLorentzVector"},
		   {"ChargedPartLV_1",          "TLorentzVector"},
		   {"ChargedPartLV_2",          "TLorentzVector"},
		   {"Charge_0",                 "int"},
		   {"Charge_1",                 "int"},
		   {"Charge_2",                 "int"},
		   {"ClusterPositions",         "std::vector<TVector3>"},
		   {"ClusterPositionsVariance", "std::vector<TVector3>"},
		   {"ClusterEnergies",          "std::vector<double>"},
		   {"ClusterEnergiesVariance",  "std::vector<double>"},
		   {"ClusterIndices",           "std::vector<int>"}};
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
	std::map<std::string, double> constArgs
		= {{"PiMass",                  0},
		   {"PiMassLowerLimit",        0},
		   {"PiMassUpperLimit",        0},
		   {"PrecisionGoalLowerLimit", 0},
		   {"PrecisionGoalUpperLimit", 0},
		   {"OmegaMass",               0},
		   {"OmegaMasswindow",         0}};

	if (not functionArgumentHandlerConst<double>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

    // Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes
		= {"std::vector<double>",   // ResultPrecisionGoals
		   "std::vector<int>",      // ResultAcceptedOmegas
		   "std::vector<double>"};  // ResultOmegaMasses
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetFittedOmegaMassVsPrecisionGoal(
		*data.getAddr<TVector3>             (args[0].first),   // VertexPosition
		*data.getAddr<TLorentzVector>       (args[1].first),   // ChargedPartLV_0
		*data.getAddr<TLorentzVector>       (args[2].first),   // ChargedPartLV_1
		*data.getAddr<TLorentzVector>       (args[3].first),   // ChargedPartLV_2
		*data.getAddr<int>                  (args[4].first),   // Charge_0
		*data.getAddr<int>                  (args[5].first),   // Charge_1
		*data.getAddr<int>                  (args[6].first),   // Charge_2
		*data.getAddr<std::vector<TVector3>>(args[7].first),   // ClusterPositions
		*data.getAddr<std::vector<TVector3>>(args[8].first),   // ClusterPositionVariances
		*data.getAddr<std::vector<double>>  (args[9].first),   // ClusterEnergies
		*data.getAddr<std::vector<double>>  (args[10].first),  // ClusterEnergieVariances
		*data.getAddr<std::vector<int>>     (args[11].first),  // ClusterIndices
		constArgs["PiMass"],                                   // PiMass
		constArgs["PiMassLowerLimit"],                         // PiMassLowerLimit
		constArgs["PiMassUpperLimit"],                         // PiMassUpperLimit
		constArgs["PrecisionGoalLowerLimit"],                  // PrecisionGoalLowerLimit
		constArgs["PrecisionGoalUpperLimit"],                  // PrecisionGoalUpperLimit
		constArgs["ErrorEstimateType"],                        // ErrorEstimateType
		constArgs["OmegaMass"],                                // OmegaMass
		constArgs["OmegaMasswindow"],                          // OmegaMasswindow
		*data.getAddr<std::vector<double>>(quantityNames[0]),  // ResultPrecisionGoals
		*data.getAddr<std::vector<int>>   (quantityNames[1]),  // ResultAccpetedOmegas
		*data.getAddr<std::vector<double>>(quantityNames[2])   // ResultOmegaMasses
	);
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
	vecPairString<std::string> args
		= {{"Pi0_0",           "TLorentzVector"},
		   {"Pi0_1",           "TLorentzVector"},
		   {"ChargedPartLV_0", "TLorentzVector"},
		   {"ChargedPartLV_1", "TLorentzVector"},
		   {"ChargedPartLV_2", "TLorentzVector"},
		   {"Charge_0",        "int"},
		   {"Charge_1",        "int"},
		   {"Charge_2",        "int"}};
	if (not functionArgumentHandler(args, function, index)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Get constant arguments
	std::map<std::string, int> constArgs = {{"UseMassSquared", 0}};
	if (not functionArgumentHandlerConst<int>(constArgs, function)) {
		std::cerr << getFunctionArgumentHandlerErrorMsg(quantityNames);
		return nullptr;
	}

	// Register output variables
	antok::Data& data = antok::ObjectManager::instance()->getData();
	const std::vector<std::string> outputVarTypes = {"std::vector<double>"};  // Result Masses or squared Masses
	if (not registerOutputVarTypes(data, quantityNames, outputVarTypes)) {
		return nullptr;
	}

	return new antok::user::cdreis::functions::GetThreePionCombinationMass(
		*data.getAddr<TLorentzVector>(args[0].first),         // Pi0LV_0
		*data.getAddr<TLorentzVector>(args[1].first),         // Pi0LV_1
		*data.getAddr<TLorentzVector>(args[2].first),         // ChargedPartLV_0
		*data.getAddr<TLorentzVector>(args[3].first),         // ChargedPartLV_1
		*data.getAddr<TLorentzVector>(args[4].first),         // ChargedPartLV_2
		*data.getAddr<int>           (args[5].first),         // Charge_0
		*data.getAddr<int>           (args[6].first),         // Charge_1
		*data.getAddr<int>           (args[7].first),         // Charge_2
		constArgs["UseMassSquared"],                          // UseMassSquared,
		*data.getAddr<std::vector<double>>(quantityNames[0])  // Result
	);
}


bool
antok::user::cdreis::registerOutputVarTypes(antok::Data&                    data,
                                            const std::vector<std::string>& quantityNames,
                                            const std::vector<std::string>& outputVarTypes)
{
	for (size_t i = 0; i < outputVarTypes.size(); ++i) {
		const std::string& typeName     = outputVarTypes[i];
		const std::string& quantityName = quantityNames[i];
		if (typeName == "int") {
			if (not data.insert<int>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
		if (typeName == "Long64_t") {
			if (not data.insert<Long64_t>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
		if (typeName == "double") {
			if (not data.insert<double>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
		if (typeName == "TVector3") {
			if (not data.insert<TVector3>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
		if (typeName == "TLorentzVector") {
			if (not data.insert<TLorentzVector>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
		if (typeName == "std::vector<int>") {
			if (not data.insert<std::vector<int>>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
		if (typeName == "std::vector<Long64_t>") {
			if (not data.insert<std::vector<Long64_t>>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
		if (typeName == "std::vector<double>") {
			if (not data.insert<std::vector<double>>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
		if (typeName == "std::vector<TVector3>") {
			if (not data.insert<std::vector<TVector3>>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
		if (typeName == "std::vector<TLorentzVector>") {
			if (not data.insert<std::vector<TLorentzVector>>(quantityName)) {
				std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityName);
				return false;
			}
		}
	}
	return true;
}
