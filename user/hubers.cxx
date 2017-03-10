
#include<hubers.h>

#include<constants.h>
#include<data.h>
#include<functions.hpp>
#include<generators_functions.h>
#include<hubers_functions.hpp>
#include<yaml_utils.hpp>
#include<iostream>
#include<fstream>

antok::Function* antok::user::hubers::getUserFunction(const YAML::Node& function,
                                                      std::vector<std::string>& quantityNames,
                                                      int index)
{
	std::string functionName = antok::YAMLUtils::getString(function["Name"]);
	antok::Function* antokFunctionPtr = 0;
	if(functionName == "sqrt")
		antokFunctionPtr = antok::user::hubers::generateSqrt(function, quantityNames, index);
	else if(functionName == "frac")
		antokFunctionPtr = antok::user::hubers::generateFrac(function, quantityNames, index);
	else if(functionName == "thetaRICH")
		antokFunctionPtr = antok::user::hubers::generateThetaRICH(function, quantityNames, index);
	else if(functionName == "getPt")
		antokFunctionPtr = antok::user::hubers::generateGetPt(function, quantityNames, index);
	else if(functionName == "EnforceEConservation")
		antokFunctionPtr = antok::user::hubers::generateEnforceEConservation(function, quantityNames, index);
	else if(functionName == "BeamNN")
		antokFunctionPtr = antok::user::hubers::generateGetNeuronalBeam(function, quantityNames, index);
	else if(functionName == "Theta")
		antokFunctionPtr = antok::user::hubers::generateGetTheta(function, quantityNames, index);
	else if(functionName == "ThetaZ")
		antokFunctionPtr = antok::user::hubers::generateGetThetaZCut(function, quantityNames, index);
	else if(functionName == "BadSpill")
		antokFunctionPtr = antok::user::hubers::generateGetBadSpill(function, quantityNames, index);
	else if(functionName == "Shift")
		antokFunctionPtr = antok::user::hubers::generateGetShifted(function, quantityNames, index);
	else if(functionName == "ScaleCluster")
		antokFunctionPtr = antok::user::hubers::generateGetScaledCluster(function, quantityNames, index);
	else if(functionName == "ClusterPosCor")
		antokFunctionPtr = antok::user::hubers::generateGetClusterPosCor(function, quantityNames, index);
	else if(functionName == "ExtrapNeutral")
		antokFunctionPtr = antok::user::hubers::generateExtrapNeutral(function, quantityNames, index);
	else if(functionName == "CleanClusters")
		antokFunctionPtr = antok::user::hubers::generateGetCleanedClusters(function, quantityNames, index);
	else if(functionName == "MaximumCluster")
		antokFunctionPtr = antok::user::hubers::generateGetMaximumCluster(function, quantityNames, index);
	else if(functionName == "getNeutralLorentzVec")
		antokFunctionPtr = antok::user::hubers::generateGetNeutralLorentzVec(function, quantityNames, index);
	else if(functionName == "FormFactor")
		antokFunctionPtr = antok::user::hubers::generateGetFormFactor(function, quantityNames, index);
	else if(functionName == "CutBgTracks")
		antokFunctionPtr = antok::user::hubers::generateGetBgTrackCut(function, quantityNames, index);
	else if(functionName == "GetClosestPi0")
		antokFunctionPtr = antok::user::hubers::generateGetClosestPi0(function, quantityNames, index);
	else if(functionName == "GetClosestPi0Pi0")
		antokFunctionPtr = antok::user::hubers::generateGetClosestPi0Pi0(function, quantityNames, index);
	else if(functionName == "GetRunSpill")
		antokFunctionPtr = antok::user::hubers::generateRunSpill(function, quantityNames, index);
	return antokFunctionPtr;
}

antok::Function* antok::user::hubers::generateSqrt(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* argAddr = data.getAddr<double>(args[0].first);

	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::Sqrt(argAddr, data.getAddr<double>(quantityName)));
};

antok::Function* antok::user::hubers::generateFrac(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Numerator", "double"));
	args.push_back(std::pair<std::string, std::string>("Denominator", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* arg1Addr = data.getAddr<double>(args[0].first);
	double* arg2Addr = data.getAddr<double>(args[1].first);

	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::Frac(arg1Addr, arg2Addr, data.getAddr<double>(quantityName)));
};

antok::Function* antok::user::hubers::generateThetaRICH(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Mom", "double"));
	args.push_back(std::pair<std::string, std::string>("Theta", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* arg1Addr = data.getAddr<double>(args[0].first);
	double* arg2Addr = data.getAddr<double>(args[1].first);

	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::thetaRICHcut(arg1Addr, arg2Addr, data.getAddr<double>(quantityName)));
};

antok::Function* antok::user::hubers::generateGetPt(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 1) {
		std::cerr<<"Need 1 names for function \"GetPt\""<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;

	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* beamLVAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* xLVAddr = data.getAddr<TLorentzVector>(args[1].first);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return (new antok::user::hubers::functions::GetPt(xLVAddr, beamLVAddr, quantityAddrs[0]));
};

antok::Function* antok::user::hubers::generateEnforceEConservation(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("LVBeam", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("LVPion", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("LVGamma", "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	TLorentzVector* beamAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* pionAddr = data.getAddr<TLorentzVector>(args[1].first);
	TLorentzVector* neutralAddr = data.getAddr<TLorentzVector>(args[2].first);
	if(not data.insert<TLorentzVector>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	double* massAddr;
	if(antok::YAMLUtils::hasNodeKey(function, "mass"))
		massAddr = antok::YAMLUtils::getAddress<double>(function["mass"]);
	else
		massAddr = new double(0);

	int* modeAddr;
	if(antok::YAMLUtils::hasNodeKey(function, "mode"))
		modeAddr = antok::YAMLUtils::getAddress<int>(function["mode"]);
	else
		modeAddr = new int(0);
	return (new antok::user::hubers::functions::EnforceEConservation(beamAddr, pionAddr, neutralAddr, massAddr, modeAddr,
	                                                                 data.getAddr<TLorentzVector>(quantityName)
	                                                                )
	       );
}

antok::Function* antok::user::hubers::generateGetNeuronalBeam(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 2) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<" needed two\"."<<std::endl;
		return 0;
	}
	std::string quantityNameD = quantityNames[0];
	std::string quantityNameLV = quantityNames[1];


	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("X" , "double"));
	args.push_back(std::pair<std::string, std::string>("Y" , "double"));
	args.push_back(std::pair<std::string, std::string>("dX", "double"));
	args.push_back(std::pair<std::string, std::string>("dY", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* xAddr  = data.getAddr<double>(args[0].first);
	double* yAddr  = data.getAddr<double>(args[1].first);
	double* dxAddr = data.getAddr<double>(args[2].first);
	double* dyAddr = data.getAddr<double>(args[3].first);

	if(not data.insert<double>(quantityNameD)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	if(not data.insert<TLorentzVector>(quantityNameLV)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	int* yearAddr;
	if(antok::YAMLUtils::hasNodeKey(function, "year")) {
		yearAddr = antok::YAMLUtils::getAddress<int>(function["year"]);
	}
	else {
		return 0;
	}

	try
	{
		function["beam2009"].as<std::string>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"beam2009\" in function \"GetNeuronalBeam\" not found" << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"beam2009\" in function \"GetNeuronalBeam\" should be of type std::string (variable \""
		          << quantityNames[0] << " and " << quantityNames[1] << "\")." << std::endl;
		return 0;
	}
	std::string *Calibration2009 = new std::string();
	(*Calibration2009)           = function["beam2009"].as<std::string>();

	std::ifstream      configFile2009;
	configFile2009.unsetf(std::ios_base::skipws);

	unsigned long int linecount2009 = (unsigned long int)std::count( std::istream_iterator<char>(configFile2009),std::istream_iterator<char>(),'\n');
	std::vector<double> calibration2009;
	calibration2009.reserve(linecount2009);

	configFile2009.open(*Calibration2009);
	if( !configFile2009.is_open() )
	{
		return 0;
	}
	double a09;
	while( configFile2009 >> a09 )
	{
		calibration2009.push_back(a09);
	}
	if( not configFile2009.eof() )
	{
		std::cout << "ERROR: Invalid input correction for value " << a09 << std::endl;
	}
	configFile2009.close();

	try
	{
		function["beam2012"].as<std::string>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"beam2012\" in function \"GetNeuronalBeam\" not found" << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"beam2012\" in function \"GetNeuronalBeam\" should be of type std::string (variable \""
		          << quantityNames[0] << " and " << quantityNames[1] << "\")." << std::endl;
		return 0;
	}
	std::string *Calibration2012 = new std::string();
	(*Calibration2012)           = function["beam2012"].as<std::string>();

	std::ifstream      configFile2012;
	configFile2012.unsetf(std::ios_base::skipws);

	unsigned long int linecount2012 = (unsigned long int)std::count( std::istream_iterator<char>(configFile2012),std::istream_iterator<char>(),'\n');
	std::vector<double> calibration2012;
	calibration2012.reserve(linecount2012);

	configFile2012.open(*Calibration2012);
	if( !configFile2012.is_open() )
	{
		return 0;
	}
	double a12;
	while( configFile2012 >> a12 )
	{
		calibration2012.push_back(a12);
	}
	if( not configFile2012.eof() )
	{
		std::cout << "ERROR: Invalid input correction for value " << a12 << std::endl;
	}
	configFile2012.close();

	return (new antok::user::hubers::functions::GetNeuronalBeam(xAddr, yAddr, dxAddr, calibration2009, calibration2012, dyAddr,
	                                                            data.getAddr<double>(quantityNameD),
	                                                            data.getAddr<TLorentzVector>(quantityNameLV),
	                                                            yearAddr
	                                                           )
	       );
}


//***********************************
//Calculates the angle theta between
//two TLorentzVectors
//***********************************
antok::Function* antok::user::hubers::generateGetTheta(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];


	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("OutLorentzVec",  "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	TLorentzVector* beamLVAddr  = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* outLVAddr   = data.getAddr<TLorentzVector>(args[1].first);
	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetTheta(beamLVAddr, outLVAddr, data.getAddr<double>(quantityName)));
}

//***********************************
//Calculates the condition for a
//theta dependend Z cut
//***********************************
antok::Function* antok::user::hubers::generateGetThetaZCut(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];


	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Theta" , "double"));
	args.push_back(std::pair<std::string, std::string>("Z" , "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* thetaAddr  = data.getAddr<double>(args[0].first);
	double* zAddr  = data.getAddr<double>(args[1].first);

	double *zMeanAddr;
	try {
		function["ZMean"].as<double>();
	} catch(const YAML::TypedBadConversion<double>& e) {
		std::cerr<<"Argument \"ZMean\" in function \""<<function["Name"]<<"\" should be of type double (variable \""<<"ZMean\")."<<std::endl;
		return 0;
	}
	zMeanAddr = new double();
	(*zMeanAddr) = function["ZMean"].as<double>();


	if(not data.insert<int>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetThetaZCut(zAddr, thetaAddr, zMeanAddr, data.getAddr<int>(quantityName)));
}

//***********************************
//result is 1 if a run is in a bad
//spill list
//***********************************
antok::Function* antok::user::hubers::generateGetBadSpill(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector< std::pair<int,int> > *badSpillList = new std::vector< std::pair<int,int> >;

	std::string fileName = antok::YAMLUtils::getString(function["fileName"]);
	std::ifstream file;
	file.open(fileName.c_str());
	char line[100];
	int run,spill;
	if(file.is_open()){
		while(file.getline(line,100)){
			int ret = sscanf(line,"%d %d",&run,&spill);
			assert( ret == 2 );
			badSpillList->push_back(std::make_pair(run,spill));
		}
	}
	else{
		std::cerr<<function["Name"]<<" could not open BadSpillList file: "<<fileName<<std::endl;
		return 0;

	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Run"  , "int"));
	args.push_back(std::pair<std::string, std::string>("Spill", "int"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	int* runAddr   = data.getAddr<int>(args[0].first);
	int* spillAddr = data.getAddr<int>(args[1].first);
	if(not data.insert<int>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetBadSpill(runAddr, spillAddr, badSpillList, data.getAddr<int>(quantityName)));
}

//***********************************
//Shifts std::vectors
//***********************************
antok::Function* antok::user::hubers::generateGetShifted(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 1) {
		std::cerr<<function["Name"]<<" needs one quantity\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];


	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Vector" , "std::vector<double>"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double>* VectorAddr  = data.getAddr<std::vector<double> >(args[0].first);


	std::vector<std::string> possiblyConstArgs;
	possiblyConstArgs.push_back("Offset");

	for(unsigned int i = 0; i < possiblyConstArgs.size(); ++i) {
		if(not antok::YAMLUtils::hasNodeKey(function, possiblyConstArgs[i])) {
			std::cerr<<"Argument \""<<possiblyConstArgs[i]<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
	}

	double* offsetAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[0]]);


	if(not data.insert<std::vector<double> >(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetShifted(VectorAddr, offsetAddr,
	                                                       data.getAddr<std::vector<double> >(quantityName)
	                                                      )
	       );
}

//***********************************
// Scale Energy of a cluster depending
// linear method in case of MC
// posdep method in case of RD
//***********************************
antok::Function* antok::user::hubers::generateGetScaledCluster(const YAML::Node& function, std::vector<std::string>& quantityNames, int index){

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("X", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("Y", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("E", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("Method", "int"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double>* XAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double>* YAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double>* EAddr = data.getAddr<std::vector<double> >(args[2].first);
	int* methodAddr = data.getAddr<int>(args[3].first);

	std::vector<std::string> possiblyConstArgs;
	possiblyConstArgs.push_back("Threshold");

	for(unsigned int i = 0; i < possiblyConstArgs.size(); ++i) {
		if(not antok::YAMLUtils::hasNodeKey(function, possiblyConstArgs[i])) {
			std::cerr<<"Argument \""<<possiblyConstArgs[i]<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
	}

	double* thresholdAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[0]]);

	if(not data.insert<std::vector<double> >(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetScaledCluster(XAddr, YAddr, EAddr, methodAddr, thresholdAddr,
	                                                             data.getAddr<std::vector<double> >(quantityName)
	                                                            )
	       );
};

//***********************************
//cleans calorimeter clusters
//and merges them dependend onn their distance
//***********************************
antok::Function* antok::user::hubers::generateGetCleanedClusters(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 5) {
		std::cerr<<function["Name"]<<" needs five quantities\"."<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorX" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorY" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorZ" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorT" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorE" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("trackX" , "double"));
	args.push_back(std::pair<std::string, std::string>("trackY" , "double"));
	args.push_back(std::pair<std::string, std::string>("trackT" , "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double>* VectorXAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double>* VectorYAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double>* VectorZAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double>* VectorTAddr = data.getAddr<std::vector<double> >(args[3].first);
	std::vector<double>* VectorEAddr = data.getAddr<std::vector<double> >(args[4].first);
	double* trackXAddr = data.getAddr<double>(args[5].first);
	double* trackYAddr = data.getAddr<double>(args[6].first);
	double* trackTAddr = data.getAddr<double>(args[7].first);

	std::vector<std::string> possiblyConstArgs;
	possiblyConstArgs.push_back("mergeDist");
	possiblyConstArgs.push_back("timeThreshold");

	for(unsigned int i = 0; i < possiblyConstArgs.size(); ++i) {
		if(not antok::YAMLUtils::hasNodeKey(function, possiblyConstArgs[i])) {
			std::cerr<<"Argument \""<<possiblyConstArgs[i]<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
	}
	double* mergeDistAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[0]]);
	double* timeThresholdAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[1]]);

	std::string resultVecX = quantityNames[0];
	std::string resultVecY = quantityNames[1];
	std::string resultVecZ = quantityNames[2];
	std::string resultVecT = quantityNames[3];
	std::string resultVecE = quantityNames[4];

	for(unsigned int i = 0; i < 5; ++i){
		if(not data.insert<std::vector<double> >(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return 0;
		}
	}
	return (new antok::user::hubers::functions::GetCleanedClusters(VectorXAddr, VectorYAddr, VectorZAddr, VectorTAddr, VectorEAddr,
	                                                              trackXAddr, trackYAddr, trackTAddr,
	                                                              mergeDistAddr, timeThresholdAddr,
	                                                              data.getAddr<std::vector<double> >(resultVecX),
	                                                              data.getAddr<std::vector<double> >(resultVecY),
	                                                              data.getAddr<std::vector<double> >(resultVecZ),
	                                                              data.getAddr<std::vector<double> >(resultVecT),
	                                                              data.getAddr<std::vector<double> >(resultVecE)
	                                                             )
	       );
}

//***********************************
//gets highest energetic calorimeter cluster
//***********************************
antok::Function* antok::user::hubers::generateGetMaximumCluster(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 6) {
		std::cerr<<function["Name"]<<" needs six quantities\"."<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorX" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorY" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorZ" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorT" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorE" , "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("trackX" , "double"));
	args.push_back(std::pair<std::string, std::string>("trackY" , "double"));
	args.push_back(std::pair<std::string, std::string>("trackT" , "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double>* VectorXAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double>* VectorYAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double>* VectorZAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double>* VectorTAddr = data.getAddr<std::vector<double> >(args[3].first);
	std::vector<double>* VectorEAddr = data.getAddr<std::vector<double> >(args[4].first);
	double* trackXAddr = data.getAddr<double>(args[5].first);
	double* trackYAddr = data.getAddr<double>(args[6].first);
	double* trackTAddr = data.getAddr<double>(args[7].first);

	std::string maximumX = quantityNames[0];
	std::string maximumY = quantityNames[1];
	std::string maximumZ = quantityNames[2];
	std::string maximumT = quantityNames[3];
	std::string maximumE = quantityNames[4];
	std::string NClus    = quantityNames[5];

	for(unsigned int i=0; i<5; ++i){
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return 0;
		}
	}
	for(unsigned int i=5; i<6; ++i){
		if(not data.insert<int>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return 0;
		}
	}
	return (new antok::user::hubers::functions::GetMaximumCluster(VectorXAddr, VectorYAddr, VectorZAddr,
	                                                              VectorTAddr, VectorEAddr,
	                                                              trackXAddr, trackYAddr, trackTAddr,
	                                                              data.getAddr<double>(maximumX),
	                                                              data.getAddr<double>(maximumY),
	                                                              data.getAddr<double>(maximumZ),
	                                                              data.getAddr<double>(maximumT),
	                                                              data.getAddr<double>(maximumE),
	                                                              data.getAddr<int>(NClus)
	                                                             )
	       );
}

//***********************************
//gets LorentzVector for cluster
//produced in a  vertex with coordinates X/Y/Z
//***********************************
antok::Function* antok::user::hubers::generateGetNeutralLorentzVec(const YAML::Node& function, std::vector<std::string>& quantityNames, int index){
	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("X", "double"));
	args.push_back(std::pair<std::string, std::string>("Y", "double"));
	args.push_back(std::pair<std::string, std::string>("Z", "double"));
	args.push_back(std::pair<std::string, std::string>("E", "double"));
	args.push_back(std::pair<std::string, std::string>("xPV", "double"));
	args.push_back(std::pair<std::string, std::string>("yPV", "double"));
	args.push_back(std::pair<std::string, std::string>("zPV", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* xAddr= data.getAddr<double>(args[0].first);
	double* yAddr= data.getAddr<double>(args[1].first);
	double* zAddr= data.getAddr<double>(args[2].first);
	double* eAddr= data.getAddr<double>(args[3].first);
	double* xPVAddr= data.getAddr<double>(args[4].first);
	double* yPVAddr= data.getAddr<double>(args[5].first);
	double* zPVAddr= data.getAddr<double>(args[6].first);


	if(not data.insert<TLorentzVector>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetNeutralLorentzVec(xAddr, yAddr, zAddr, eAddr, xPVAddr, yPVAddr, zPVAddr,
	                                                                 data.getAddr<TLorentzVector>(quantityName)
	                                                                )
	       );

};

//***********************************
//Calculates the Form Factor correction
//for Nickel from -Q2-MCTRUTH
//returns true/false
//a lot of things are hard coded for Nickel
//***********************************
antok::Function* antok::user::hubers::generateGetFormFactor(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Arg", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* argAddr = data.getAddr<double>(args[0].first);

	if(not data.insert<int>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::FormFactor(argAddr, data.getAddr<int>(quantityName)));
};


//***********************************
//Calculates bgTrack Cut
//returns true/false
//***********************************
antok::Function* antok::user::hubers::generateGetBgTrackCut(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("evTime", "double"));
	args.push_back(std::pair<std::string, std::string>("tracksP", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("tracksT", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("tracksTSigma", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("tracksZfirst", "std::vector<double>"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	double* evTimeAddr = data.getAddr<double>(args[0].first);
	std::vector<double>* tracksPAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double>* tracksTAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double>* tracksTSigmaAddr = data.getAddr<std::vector<double> >(args[3].first);
	std::vector<double>* tracksZfirstAddr = data.getAddr<std::vector<double> >(args[4].first);

	if(not data.insert<int>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::BgTracks(evTimeAddr, tracksPAddr, tracksTAddr,
	                                                     tracksTSigmaAddr, tracksZfirstAddr,
	                                                     data.getAddr<int>(quantityName)
	                                                    )
	       );
};

//***********************************
//Gets best pi0 pair
//gives an LV and the mass
//***********************************
antok::Function* antok::user::hubers::generateGetClosestPi0(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 2) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<" needed two\"."<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorE", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorX", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorY", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorZ", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("xPV", "double"));
	args.push_back(std::pair<std::string, std::string>("yPV", "double"));
	args.push_back(std::pair<std::string, std::string>("zPV", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double>* VectorEAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double>* VectorXAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double>* VectorYAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double>* VectorZAddr = data.getAddr<std::vector<double> >(args[3].first);
	double* xPVAddr = data.getAddr<double>(args[4].first);
	double* yPVAddr = data.getAddr<double>(args[5].first);
	double* zPVAddr = data.getAddr<double>(args[6].first);

	if(not data.insert<TLorentzVector>(quantityNames[0])) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[0]);
		return 0;
	}
	if(not data.insert<double>(quantityNames[1])) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[1]);
		return 0;
	}


		if(not antok::YAMLUtils::hasNodeKey(function, "selectedMass")) {
			std::cerr<<"Argument \""<<"selectedMass"<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
		double* selectedMassAddr = antok::YAMLUtils::getAddress<double>(function["selectedMass"]);

	return (new antok::user::hubers::functions::GetClosestPi0(VectorEAddr, VectorXAddr, VectorYAddr, VectorZAddr,
	                                                              xPVAddr, yPVAddr, zPVAddr, selectedMassAddr,
	                                                              data.getAddr<TLorentzVector>(quantityNames[0]),
	                                                              data.getAddr<double>(quantityNames[1])
	                                                             )
	       );
}

//***********************************
//Gets best pi0pi0 pair
//gives an LV and the mass
//***********************************
antok::Function* antok::user::hubers::generateGetClosestPi0Pi0(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 4) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<" needed four\"."<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorE", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorX", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorY", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("VectorZ", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("xPV", "double"));
	args.push_back(std::pair<std::string, std::string>("yPV", "double"));
	args.push_back(std::pair<std::string, std::string>("zPV", "double"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double>* VectorEAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double>* VectorXAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double>* VectorYAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double>* VectorZAddr = data.getAddr<std::vector<double> >(args[3].first);
	double* xPVAddr = data.getAddr<double>(args[4].first);
	double* yPVAddr = data.getAddr<double>(args[5].first);
	double* zPVAddr = data.getAddr<double>(args[6].first);

	if(not data.insert<TLorentzVector>(quantityNames[0])) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[0]);
		return 0;
	}
	if(not data.insert<double>(quantityNames[1])) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[1]);
		return 0;
	}
	if(not data.insert<TLorentzVector>(quantityNames[2])) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[2]);
		return 0;
	}
	if(not data.insert<double>(quantityNames[3])) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[3]);
		return 0;
	}


		if(not antok::YAMLUtils::hasNodeKey(function, "selectedMass")) {
			std::cerr<<"Argument \""<<"selectedMass"<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
		double* selectedMassAddr = antok::YAMLUtils::getAddress<double>(function["selectedMass"]);

	return (new antok::user::hubers::functions::GetClosestPi0Pi0(VectorEAddr, VectorXAddr, VectorYAddr, VectorZAddr,
	                                                              xPVAddr, yPVAddr, zPVAddr, selectedMassAddr,
	                                                              data.getAddr<TLorentzVector>(quantityNames[0]),
	                                                              data.getAddr<double>(quantityNames[1]),
	                                                              data.getAddr<TLorentzVector>(quantityNames[2]),
	                                                              data.getAddr<double>(quantityNames[3])
	                                                             )
	       );
}

//***********************************
//Calculates the condition for a
//theta dependend Z cut
//***********************************
antok::Function* antok::user::hubers::generateRunSpill(const YAML::Node& function, std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() > 1) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];


	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Run" , "int"));
	args.push_back(std::pair<std::string, std::string>("Spill" , "int"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	int* runAddr  = data.getAddr<int>(args[0].first);
	int* spillAddr  = data.getAddr<int>(args[1].first);

// 	double *FactorAddr;
// // 	try {
// // 		function["Factor"].as<double>();
// // 	} catch(const YAML::TypedBadConversion<double>& e) {
// // 		std::cerr<<"Argument \"Factor\" in function \""<<function["Name"]<<"\" should be of type double (variable \""<<"Factor\")."<<std::endl;
// // 		return 0;
// // 	}
// 	FactorAddr = new double();
// // 	(*FactorAddr) = function["ZMean"].as<double>();
//   *FactorAddr =200;
	std::vector<std::string> possiblyConstArgs;
	possiblyConstArgs.push_back("Factor");

	for(unsigned int i = 0; i < possiblyConstArgs.size(); ++i) {
		if(not antok::YAMLUtils::hasNodeKey(function, possiblyConstArgs[i])) {
			std::cerr<<"Argument \""<<possiblyConstArgs[i]<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
	}

	double* FactorAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[0]]);


	if(not data.insert<double>(quantityName)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::hubers::functions::GetRunSpill(runAddr, spillAddr, FactorAddr, data.getAddr<double>(quantityName)));
}

//***********************************
// Scale Energy of a cluster depending
// linear method in case of MC
// posdep method in case of RD
//***********************************
antok::Function* antok::user::hubers::generateGetClusterPosCor(const YAML::Node& function, std::vector<std::string>& quantityNames, int index){

	if(quantityNames.size() > 2) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityNameX = quantityNames[0];
	std::string quantityNameY = quantityNames[1];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("X", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("Y", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("Xic", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("Yic", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("E", "std::vector<double>"));
// 	args.push_back(std::pair<std::string, std::string>("Method", "int"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double>* XAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double>* YAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double>* XicAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double>* YicAddr = data.getAddr<std::vector<double> >(args[3].first);
	std::vector<double>* EAddr = data.getAddr<std::vector<double> >(args[4].first);
	int* methodAddr = 0;//data.getAddr<int>(args[3].first);

	std::vector<std::string> possiblyConstArgs;
	possiblyConstArgs.push_back("Threshold");

	for(unsigned int i = 0; i < possiblyConstArgs.size(); ++i) {
		if(not antok::YAMLUtils::hasNodeKey(function, possiblyConstArgs[i])) {
			std::cerr<<"Argument \""<<possiblyConstArgs[i]<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
	}

	double* thresholdAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[0]]);

	if(not data.insert<std::vector<double> >(quantityNameX)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	if(not data.insert<std::vector<double> >(quantityNameY)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
  }
	return (new antok::user::hubers::functions::GetClusterPosCor(XAddr, YAddr, XicAddr, YicAddr, EAddr, methodAddr, thresholdAddr,
	                                                             data.getAddr<std::vector<double> >(quantityNameX),
	                                                             data.getAddr<std::vector<double> >(quantityNameY)
	                                                            )
	       );
};

antok::Function* antok::user::hubers::generateExtrapNeutral(const YAML::Node& function, std::vector<std::string>& quantityNames, int index){

	if(quantityNames.size() > 2) {
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}
	std::string quantityNameX = quantityNames[0];
	std::string quantityNameY = quantityNames[1];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("Pi", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("Beam", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("X", "double"));
	args.push_back(std::pair<std::string, std::string>("Y", "double"));
// 	args.push_back(std::pair<std::string, std::string>("E", "std::vector<double>"));
// 	args.push_back(std::pair<std::string, std::string>("Method", "int"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

  TLorentzVector* piLV = data.getAddr<TLorentzVector>(args[0].first);
  TLorentzVector* beamLV = data.getAddr<TLorentzVector>(args[1].first);
	double* XAddr = data.getAddr<double>(args[2].first);
	double* YAddr = data.getAddr<double>(args[3].first);

	std::vector<std::string> possiblyConstArgs;
	possiblyConstArgs.push_back("Z");

	for(unsigned int i = 0; i < possiblyConstArgs.size(); ++i) {
		if(not antok::YAMLUtils::hasNodeKey(function, possiblyConstArgs[i])) {
			std::cerr<<"Argument \""<<possiblyConstArgs[i]<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
	}

	double* ZAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[0]]);

	if(not data.insert<double>(quantityNameX)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}
	if(not data.insert<double >(quantityNameY)) {
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
  }
	return (new antok::user::hubers::functions::ExtrapNeutral(XAddr, YAddr, ZAddr, piLV, beamLV,
	                                                             data.getAddr<double>(quantityNameX),
	                                                             data.getAddr<double>(quantityNameY)
	                                                            )
	       );
};





double antok::user::hubers::IntraCellX( double cx ){
	static const double fMeanECALX = 0.0;
	const double xbound = 3.060 + 9.575 - fMeanECALX - 3.5 * 3.83;  // x-position of cell boundary closest to 0
	return fabs( remainder( cx - xbound, 3.83 ) );
}

double antok::user::hubers::IntraCellY( double cy ){
	static const double fMeanECALY = 0.00;
	const double ybound = 0.150 + 5.745 - fMeanECALY - 1.5 * 3.83;  // y-position of cell boundary closest to 0
	return fabs( remainder( cy - ybound, 3.83 ) );
}

double antok::user::hubers::PEDepGammaCorrection( double Egamma, double cx, double cy ){
	// obtained from SelectorHelper::FitPosEnergyCorrection()
	const double x190 = Egamma/190.;

	const double p0 = -9.57965 - 6.42201 * atan( 7.38429*(x190-0.749578) )  +0.74;
	const double p1 =  1.89692 + 1.2888  * atan( 8.79757*(x190-0.738527) );
	const double p2 =  1.61223 + 1.13902 * atan( 9.43193*(x190-0.759991) );

	const double p4 = 2.57235  - x190 * 15.9715;
	const double p5 = 0.214072 - x190 *  0.202193;
	const double p6 = 1.;
	const double p7 = 0.97;

	double x = IntraCellX( cx );
	double y = IntraCellY( cy );

	double retval = Egamma -
	       ( p0 + p1 * antok::sqr(x - 3.83/2.) +
	         p2 * antok::sqr(y - 3.83/2.) +
	         p4 * exp( -0.5 * ( (antok::sqr(x - p6) + antok::sqr(y - p7)) / antok::sqr(p5) ) ) );
			 return retval;

}

double antok::user::hubers::LinearGammaCorrection( double Egamma ) {
	static const double fLinEcorr0 = -1.65+2;//-1.21919-0.8;
	static const double fLinEcorr1 = 0.033045;
// 	static const double fLinEcorr0 = -0.859936;
// 	static const double fLinEcorr1 = 5.76984;
	const double corr = fLinEcorr0 + fLinEcorr1 * Egamma;
	Egamma -= corr;
	return Egamma ;
}
