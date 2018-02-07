
#include<cdreis.h>

#include<constants.h>
#include<data.h>
#include<functions.hpp>
#include<generators_functions.h>
#include<cdreis_functions.hpp>
#include<yaml_utils.hpp>

#include "fstream"

antok::Function *antok::user::cdreis::getUserFunction(const YAML::Node &function,
                                                      std::vector<std::string> &quantityNames,
                                                      int index) {
	std::string functionName = antok::YAMLUtils::getString(function["Name"]);
	antok::Function *antokFunctionPtr = 0;
	if( functionName == "GetRecoilLorentzVec" )
		antokFunctionPtr = antok::user::cdreis::generateGetRecoilLorentzVec(function, quantityNames, index);
	else if( functionName == "getPhotonLorentzVecs")
		antokFunctionPtr = antok::user::cdreis::generateGetPhotonLorentzVecs(function, quantityNames, index);
	else if( functionName == "getVectorLorentzVectorAttributes" )
		antokFunctionPtr = antok::user::cdreis::generateGetVectorLorentzVectorAttributes(function, quantityNames, index);
	else if( functionName == "getCleanedEcalClusters" )
		antokFunctionPtr = antok::user::cdreis::generateGetCleanedEcalClusters(function, quantityNames, index);
	else if( functionName == "getPi0Pair" )
		antokFunctionPtr = antok::user::cdreis::generateGetPi0Pair(function, quantityNames, index);
	else if( functionName == "getOmega" )
		antokFunctionPtr = antok::user::cdreis::generateGetOmega(function, quantityNames, index);
	else if( functionName == "getECALCorrectedEnergy" )
		antokFunctionPtr = antok::user::cdreis::generateGetECALCorrectedEnergy(function, quantityNames, index);
	else if( functionName == "getECALCorrectedTiming" )
		antokFunctionPtr = antok::user::cdreis::generateGetECALCorrectedTiming(function, quantityNames, index);
	else if( functionName == "getPhotonPairParticles" )
		antokFunctionPtr = antok::user::cdreis::generateGetPhotonPairParticles(function, quantityNames, index);

	return antokFunctionPtr;
}

antok::Function *antok::user::cdreis::generateGetRecoilLorentzVec( const YAML::Node         &function,
                                                                   std::vector<std::string> &quantityNames,
                                                                   int index )
{
	if (quantityNames.size() > 1)
	{
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>( "BeamLorentzVec", "TLorentzVector") );
	args.push_back(std::pair<std::string, std::string>( "XLorentzVec"   , "TLorentzVector") );

	if (not antok::generators::functionArgumentHandler(args, function, index))
	{
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	try
	{
		function["RecoilMass"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"RecoilMass\" in function \"GetRecoilLorentzVec\" not found for calculation of variable \"" << quantityName << "\"" << std::endl;
		return 0;
	}
	catch (const YAML::TypedBadConversion<double> &e)
	{
		std::cerr << "Argument \"RecoilMass\" in function \"GetRecoilLorentzVec\" should be of type double variable for calculation of variable \"" << quantityName << "\"" << std::endl;
		return 0;
	}
	double *RecoilMass = new double();
	(*RecoilMass) = function["RecoilMass"].as<double>();

	antok::Data &data = antok::ObjectManager::instance()->getData();

	TLorentzVector *BeamLorentzVec = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector *XLorentzVec    = data.getAddr<TLorentzVector>(args[1].first);

	if (not data.insert<TLorentzVector>(quantityName))
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	return (new antok::user::cdreis::functions::GetRecoilLorentzVec( BeamLorentzVec,
	                                                                 XLorentzVec,
	                                                                 RecoilMass,
	                                                                 data.getAddr<TLorentzVector>(quantityName)));
};

antok::Function *antok::user::cdreis::generateGetPhotonLorentzVecs( const YAML::Node         &function,
                                                                    std::vector<std::string> &quantityNames,
                                                                    int index )
{
	if (quantityNames.size() > 2)
	{
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>( "VectorX"   , "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "VectorY"   , "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "VectorZ"   , "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "VectorE"   , "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "VectorTime", "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "xPV"       , "double"              ));
	args.push_back(std::pair<std::string, std::string>( "yPV"       , "double"              ));
	args.push_back(std::pair<std::string, std::string>( "zPV"       , "double"              ));

	if (not antok::generators::functionArgumentHandler(args, function, index))
	{
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<double> *VectorXAddr    = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double> *VectorYAddr    = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double> *VectorZAddr    = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double> *VectorEAddr    = data.getAddr<std::vector<double> >(args[3].first);
	std::vector<double> *VectorTimeAddr = data.getAddr<std::vector<double> >(args[4].first);
	double *xPVAddr                     = data.getAddr<double>              (args[5].first);
	double *yPVAddr                     = data.getAddr<double>              (args[6].first);
	double *zPVAddr                     = data.getAddr<double>              (args[7].first);

	std::string resultVec       = quantityNames[0];
	std::string resultECALIndex = quantityNames[1];

	if (not data.insert<std::vector<TLorentzVector> >(resultVec))
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[0]);
		return 0;
	}

	if (not data.insert<std::vector<int> >(resultECALIndex))
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[1]);
		return 0;
	}

	try
	{
		function["RangeECAL1"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"RangeECAL1\" in function \"GetPhotonLorentzVecs\" not found for calculation of variables \"" << quantityNames[0] << "\" and \""
		                                                                                                                      << quantityNames[1] << "\"" << std::endl;
		return 0;
	}
	catch (const YAML::TypedBadConversion<double> &e)
	{
		std::cerr << "Argument \"RangeECAL1\" in function \"GetPhotonLorentzVecs\" should be of type double for calculation of variables \"" << quantityNames[0] << "\" and \""
		                                                                                                                                     << quantityNames[1] << "\"" << std::endl;
		return 0;
	}
	double *RangeECAL1 = new double();
	(*RangeECAL1) = function["RangeECAL1"].as<double>();

	return (new antok::user::cdreis::functions::GetPhotonLorentzVecs( VectorXAddr,
	                                                                  VectorYAddr,
	                                                                  VectorZAddr,
	                                                                  VectorEAddr,
	                                                                  VectorTimeAddr,
	                                                                  xPVAddr,
	                                                                  yPVAddr,
	                                                                  zPVAddr,
	                                                                  RangeECAL1,
	                                                                  data.getAddr<std::vector<TLorentzVector> >(resultVec),
	                                                                  data.getAddr<std::vector<int> >           (resultECALIndex))
	);
};


antok::Function *antok::user::cdreis::generateGetCleanedEcalClusters( const YAML::Node         &function,
                                                                      std::vector<std::string> &quantityNames,
                                                                      int index )
{
	if (quantityNames.size() > 6)
	{
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>( "VectorX", "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "VectorY", "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "VectorZ", "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "VectorE", "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "VectorT", "std::vector<double>" ));

	if (not antok::generators::functionArgumentHandler(args, function, index))
	{
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<double> *VectorXAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double> *VectorYAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double> *VectorZAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double> *VectorEAddr = data.getAddr<std::vector<double> >(args[3].first);
	std::vector<double> *VectorTAddr = data.getAddr<std::vector<double> >(args[4].first);

	std::string resultVectorX     = quantityNames[0];
	std::string resultVectorY     = quantityNames[1];
	std::string resultVectorZ     = quantityNames[2];
	std::string resultVectorE     = quantityNames[3];
	std::string resultVectorT     = quantityNames[4];
	std::string resultVectorIndex = quantityNames[5];

	for( unsigned int i = 0; i < quantityNames.size() - 1; i++ )
	{
		if( not data.insert<std::vector<double> >(quantityNames[i]) )
		{
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return 0;
		}
	}
	if( not data.insert<std::vector<int> >(quantityNames[5]) )
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[5]);
		return 0;
	}

	try
	{
		function["RangeECAL1"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"RangeECAL1\" in function \"GetCleanedEcalClusters\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                        << quantityNames[1] << "\", "    << "\""
		                                                                                                                        << quantityNames[2] << "\", "    << "\""
		                                                                                                                        << quantityNames[3] << "\", "    << "\""
		                                                                                                                        << quantityNames[4] << "\" and " << "\""
		                                                                                                                        << quantityNames[5] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"RangeECAL1\" in function \"GetCleanedEcalClusters\" should be type of double for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                                       << quantityNames[1] << "\", "    << "\""
		                                                                                                                                       << quantityNames[2] << "\", "    << "\""
		                                                                                                                                       << quantityNames[3] << "\", "    << "\""
		                                                                                                                                       << quantityNames[4] << "\" and " << "\""
		                                                                                                                                       << quantityNames[5] << "\""      << std::endl;
		return 0;
	}
	double *RangeECAL1 = new double();
	(*RangeECAL1)      = function["RangeECAL1"].as<double>();

	try
	{
		function["ThresholdEnergyECAL1"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ThresholdEnergyECAL1\" in function \"GetCleanedEcalClusters\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                                  << quantityNames[1] << "\", "    << "\""
		                                                                                                                                  << quantityNames[2] << "\", "    << "\""
		                                                                                                                                  << quantityNames[3] << "\", "    << "\""
		                                                                                                                                  << quantityNames[4] << "\" and " << "\""
		                                                                                                                                  << quantityNames[5] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"ThresholdEnergyECAL1\" in function \"GetCleanedEcalClusters\" should be type of double for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[1] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[2] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[3] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[4] << "\" and " << "\""
		                                                                                                                                                 << quantityNames[5] << "\""      << std::endl;;
		return 0;
	}
	double *ThresholdEnergyECAL1 = new double();
	(*ThresholdEnergyECAL1)      = function["ThresholdEnergyECAL1"].as<double>();

	try
	{
		function["ThresholdTimingECAL1"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ThresholdTimingECAL1\" in function \"GetCleanedEcalClusters\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                                  << quantityNames[1] << "\", "    << "\""
		                                                                                                                                  << quantityNames[2] << "\", "    << "\""
		                                                                                                                                  << quantityNames[3] << "\", "    << "\""
		                                                                                                                                  << quantityNames[4] << "\" and " << "\""
		                                                                                                                                  << quantityNames[5] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"ThresholdTimingECAL1\" in function \"GetCleanedEcalClusters\" should be type of double for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[1] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[2] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[3] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[4] << "\" and " << "\""
		                                                                                                                                                 << quantityNames[5] << "\""      << std::endl;;
		return 0;
	}
	double *ThresholdTimingECAL1 = new double();
	(*ThresholdTimingECAL1)      = function["ThresholdTimingECAL1"].as<double>();

	try
	{
		function["ThresholdEnergyECAL2"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ThresholdEnergyECAL2\" in function \"GetCleanedEcalClusters\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                                  << quantityNames[1] << "\", "    << "\""
		                                                                                                                                  << quantityNames[2] << "\", "    << "\""
		                                                                                                                                  << quantityNames[3] << "\", "    << "\""
		                                                                                                                                  << quantityNames[4] << "\" and " << "\""
		                                                                                                                                  << quantityNames[5] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"ThresholdEnergyECAL2\" in function \"GetCleanedEcalClusters\" should be type of double for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[1] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[2] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[3] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[4] << "\" and " << "\""
		                                                                                                                                                 << quantityNames[5] << "\""      << std::endl;;
		return 0;
	}
	double *ThresholdEnergyECAL2 = new double();
	(*ThresholdEnergyECAL2)      = function["ThresholdEnergyECAL2"].as<double>();

	try
	{
		function["ThresholdTimingECAL2"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ThresholdTimingECAL2\" in function \"GetCleanedEcalClusters\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                                  << quantityNames[1] << "\", "    << "\""
		                                                                                                                                  << quantityNames[2] << "\", "    << "\""
		                                                                                                                                  << quantityNames[3] << "\", "    << "\""
		                                                                                                                                  << quantityNames[4] << "\" and " << "\""
		                                                                                                                                  << quantityNames[5] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"ThresholdTimingECAL2\" in function \"GetCleanedEcalClusters\" should be type of double for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[1] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[2] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[3] << "\", "    << "\""
		                                                                                                                                                 << quantityNames[4] << "\" and " << "\""
		                                                                                                                                                 << quantityNames[5] << "\""      << std::endl;;
		return 0;
		return 0;
	}
	double *ThresholdTimingECAL2 = new double();
	(*ThresholdTimingECAL2)      = function["ThresholdTimingECAL2"].as<double>();


	return( new antok::user::cdreis::functions::GetCleanedEcalClusters( VectorXAddr,
	                                                                    VectorYAddr,
	                                                                    VectorZAddr,
	                                                                    VectorEAddr,
	                                                                    VectorTAddr,
	                                                                    RangeECAL1,
	                                                                    ThresholdEnergyECAL1,
	                                                                    ThresholdTimingECAL1,
	                                                                    ThresholdEnergyECAL2,
	                                                                    ThresholdTimingECAL2,
	                                                                    data.getAddr<std::vector<double> >(resultVectorX),
	                                                                    data.getAddr<std::vector<double> >(resultVectorY),
	                                                                    data.getAddr<std::vector<double> >(resultVectorZ),
	                                                                    data.getAddr<std::vector<double> >(resultVectorE),
	                                                                    data.getAddr<std::vector<double> >(resultVectorT),
	                                                                    data.getAddr<std::vector<int> >   (resultVectorIndex))
	);
};


antok::Function *antok::user::cdreis::generateGetVectorLorentzVectorAttributes( const YAML::Node         &function,
                                                                                std::vector<std::string> &quantityNames,
                                                                                int                      index )
{
	if (quantityNames.size() > 7)
	{
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("VectorLV", "std::vector<TLorentzVector>"));

	if (not antok::generators::functionArgumentHandler(args, function, index))
	{
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<TLorentzVector> *VectorLVAddr = data.getAddr<std::vector<TLorentzVector> >(args[0].first);

	std::string resultVecX     = quantityNames[0];
	std::string resultVecY     = quantityNames[1];
	std::string resultVecZ     = quantityNames[2];
	std::string resultVecE     = quantityNames[3];
	std::string resultVecTheta = quantityNames[4];
	std::string resultVecPhi   = quantityNames[5];
	std::string resultVecMag   = quantityNames[6];

	for (unsigned int i = 0; i < quantityNames.size(); i++)
	{
		if (not data.insert<std::vector<double> >(quantityNames[i]))
		{
			std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames[i]);
			return 0;
		}
	}

	return (new antok::user::cdreis::functions::GetVectorLorentzVectorAttributes( VectorLVAddr,
	                                                                              data.getAddr<std::vector<double> >(resultVecX),
	                                                                              data.getAddr<std::vector<double> >(resultVecY),
	                                                                              data.getAddr<std::vector<double> >(resultVecZ),
	                                                                              data.getAddr<std::vector<double> >(resultVecE),
	                                                                              data.getAddr<std::vector<double> >(resultVecTheta),
	                                                                              data.getAddr<std::vector<double> >(resultVecPhi),
	                                                                              data.getAddr<std::vector<double> >(resultVecMag))
	);
};

antok::Function *antok::user::cdreis::generateGetPi0Pair( const YAML::Node         &function,
                                                          std::vector<std::string> &quantityNames,
                                                          int                       index )
{
	if (quantityNames.size() > 4)
	{
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>( "VectorLV" , "std::vector<TLorentzVector>" ));
	args.push_back(std::pair<std::string, std::string>( "ECALIndex", "std::vector<int>"            ));

	if( not antok::generators::functionArgumentHandler(args, function, index) )
	{
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<TLorentzVector> *VectorLVAddr = data.getAddr<std::vector<TLorentzVector> >(args[0].first);
	std::vector<int> *ECALIndexAddr           = data.getAddr<std::vector<int> >           (args[1].first);

	std::string resultVecLV    = quantityNames[0];
	std::string resultVecLV0   = quantityNames[1];
	std::string resultVecLV1   = quantityNames[2];
	std::string resultGoodPair = quantityNames[3];

	if (not data.insert<std::vector<TLorentzVector> >(resultVecLV))
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultVecLV);
		return 0;
	}

	if (not data.insert<TLorentzVector>(resultVecLV0))
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultVecLV0);
		return 0;
	}

	if (not data.insert<TLorentzVector>(resultVecLV1))
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultVecLV1);
		return 0;
	}

	if (not data.insert<int>(resultGoodPair))
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultGoodPair);
		return 0;
	}

	try
	{
		function["Mass"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"Mass\" in function \"GetPi0Pair\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                      << quantityNames[1] << "\", "    << "\""
		                                                                                                      << quantityNames[2] << "\" and " << "\""
		                                                                                                      << quantityNames[3] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"ECALResolution\" in function \"GetRecoilLorentzVec\" should be of type double for calculation of variables \"" << quantityNames[0] << "\", "    <<  "\""
		                                                                                                                                        << quantityNames[1] << "\", "    <<  "\""
		                                                                                                                                        << quantityNames[2] << "\" and " <<  "\""
		                                                                                                                                        << quantityNames[3] << "\""      << std::endl;
		return 0;
	}
	double *Mass = new double();
	(*Mass)      = function["Mass"].as<double>();

	try
	{
		function["ECALResolution"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ECALResolution\" in function \"GetPi0Pair\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                << quantityNames[1] << "\", "    << "\""
		                                                                                                                << quantityNames[2] << "\" and " << "\""
		                                                                                                                << quantityNames[3] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<double> &exception )
	{
		std::cerr << "Argument \"ECALResolution\" in function \"GetRecoilLorentzVec\" should be of type double for calculatiom of variables \"" << quantityNames[0] << "\", "    <<  "\""
		                                                                                                                                        << quantityNames[1] << "\", "    <<  "\""
		                                                                                                                                        << quantityNames[2] << "\" and " <<  "\""
		                                                                                                                                        << quantityNames[3] << "\""      << std::endl;
		return 0;
	}
	double* ECALResolution = new double();
	(*ECALResolution)      = function["ECALResolution"].as<double>();

	try
	{
		function["ECAL1Resolution"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ECAL1Resolution\" in function \"GetPi0Pair\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                 << quantityNames[1] << "\", "    << "\""
		                                                                                                                 << quantityNames[2] << "\" and " << "\""
		                                                                                                                 << quantityNames[3] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<double> &exception )
	{
		std::cerr << "Argument \"ECAL1Resolution\" in function \"GetRecoilLorentzVec\" should be of type double for calculatiom of variables \"" << quantityNames[0] << "\", "    <<  "\""
		                                                                                                                                         << quantityNames[1] << "\", "    <<  "\""
		                                                                                                                                         << quantityNames[2] << "\" and " <<  "\""
		                                                                                                                                         << quantityNames[3] << "\""      << std::endl;
		return 0;
	}
	double* ECAL1Resolution = new double();
	(*ECAL1Resolution)      = function["ECAL1Resolution"].as<double>();

	try
	{
		function["ECAL2Resolution"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ECAL2Resolution\" in function \"GetPi0Pair\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                 << quantityNames[1] << "\", "    << "\""
		                                                                                                                 << quantityNames[2] << "\" and " << "\""
		                                                                                                                 << quantityNames[3] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<double> &exception )
	{
		std::cerr << "Argument \"ECAL2Resolution\" in function \"GetPi0Pair\" not found for calculation of variables \"" << quantityNames[0] << "\", "    << "\""
		                                                                                                                 << quantityNames[1] << "\", "    << "\""
		                                                                                                                 << quantityNames[2] << "\" and " << "\""
		                                                                                                                 << quantityNames[3] << "\""      << std::endl;
		return 0;
	}
	double* ECAL2Resolution = new double();
	(*ECAL2Resolution)      = function["ECAL2Resolution"].as<double>();

	return (new antok::user::cdreis::functions::GetPi0Pair(VectorLVAddr,
	                                                       ECALIndexAddr,
	                                                       Mass,
	                                                       ECALResolution,
	                                                       ECAL1Resolution,
	                                                       ECAL2Resolution,
	                                                       data.getAddr<std::vector<TLorentzVector> >(resultVecLV),
	                                                       data.getAddr<TLorentzVector>(resultVecLV0),
	                                                       data.getAddr<TLorentzVector>(resultVecLV1),
	                                                       data.getAddr<int>(resultGoodPair))
	);
};

antok::Function* antok::user::cdreis::generateGetOmega( const YAML::Node         &function,
                                                        std::vector<std::string> &quantityNames,
                                                        int index )
{
	if(quantityNames.size() > 4)
	{
		std::cerr<<"Too many names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>( "Pi0_0"     , "TLorentzVector" ));
	args.push_back(std::pair<std::string, std::string>( "Pi0_1"     , "TLorentzVector" ));
	args.push_back(std::pair<std::string, std::string>( "Scattered0", "TLorentzVector" ));
	args.push_back(std::pair<std::string, std::string>( "Scattered1", "TLorentzVector" ));
	args.push_back(std::pair<std::string, std::string>( "Scattered2", "TLorentzVector" ));
	args.push_back(std::pair<std::string, std::string>( "Charge0"   , "int"            ));
	args.push_back(std::pair<std::string, std::string>( "Charge1"   , "int"            ));
	args.push_back(std::pair<std::string, std::string>( "Charge2"   , "int"            ));

	if(not antok::generators::functionArgumentHandler(args, function, index))
	{
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	TLorentzVector* Pi0_0Addr      = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* Pi0_1Addr      = data.getAddr<TLorentzVector>(args[1].first);
	TLorentzVector* Scattered0Addr = data.getAddr<TLorentzVector>(args[2].first);
	TLorentzVector* Scattered1Addr = data.getAddr<TLorentzVector>(args[3].first);
	TLorentzVector* Scattered2Addr = data.getAddr<TLorentzVector>(args[4].first);

	int* Charged0Addr = data.getAddr<int>(args[5].first);
	int* Charged1Addr = data.getAddr<int>(args[6].first);
	int* Charged2Addr = data.getAddr<int>(args[7].first);

	std::string resultOmegaLV        = quantityNames[0];
	std::string resultAccepted       = quantityNames[1];
	std::string resultNotUsedPi0     = quantityNames[2];
	std::string resultNotUsedPiMinus = quantityNames[3];

	if(not data.insert<TLorentzVector>(quantityNames[0]))
	{
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[0]);
		return 0;
	}

	if(not data.insert<int>(quantityNames[1]))
	{
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[1]);
		return 0;
	}

	if(not data.insert<TLorentzVector>(quantityNames[2]))
	{
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[2]);
		return 0;
	}

	if(not data.insert<TLorentzVector>(quantityNames[3]))
	{
		std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames[3]);
		return 0;
	}

	try
	{
		function["Mass"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"Mass\" in function \"GetOmega\" not found for calculation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                    << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"Mass\" in function \"GetOmega\" should be type of double for calculation of variables" << quantityNames[0] << "\" and " << "\""
		                                                                                                                << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	double *Mass = new double();
	(*Mass)      = function["Mass"].as<double>();

	try
	{
		function["ResolutionOmega"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ResolutionOmega\" in function \"GetOmega\" not found for calculation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                               << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<double> &exception )
	{
		std::cerr << "Argument \"ResolutionOmega\" in function \"GetOmega\" should be type of double for calculation of variables" << quantityNames[0] << "\" and " << "\""
		                                                                                                                           << quantityNames[1] << "\""      << std::endl;

		return 0;
	}
	double* ResolutionOmega = new double();
	(*ResolutionOmega)      = function["ResolutionOmega"].as<double>();


	return (new antok::user::cdreis::functions::GetOmega( Pi0_0Addr,
	                                                      Pi0_1Addr,
	                                                      Scattered0Addr,
	                                                      Scattered1Addr,
	                                                      Scattered2Addr,
	                                                      Charged0Addr,
	                                                      Charged1Addr,
	                                                      Charged2Addr,
	                                                      Mass,
	                                                      ResolutionOmega,
	                                                      data.getAddr<TLorentzVector>(quantityNames[0]),
	                                                      data.getAddr<int>           (quantityNames[1]),
	                                                      data.getAddr<TLorentzVector>(quantityNames[2]),
	                                                      data.getAddr<TLorentzVector>(quantityNames[3]) )
	);
};

antok::Function * antok::user::cdreis::generateGetECALCorrectedEnergy( const YAML::Node         &function,
                                                                       std::vector<std::string> &quantityNames,
                                                                       int                       index )
{
	if( quantityNames.size() > 1 )
	{
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>( "Energy"   , "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "ClusterZ" , "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "RunNumber", "int"                 ));

	if( not antok::generators::functionArgumentHandler(args, function, index) )
	{
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<double>* Energy    = data.getAddr<std::vector<double>>(args[0].first);
	std::vector<double>* ClusterZ  = data.getAddr<std::vector<double>>(args[1].first);
	int*                 RunNumber = data.getAddr<int>                (args[2].first);

	if( not data.insert<std::vector<double>>(quantityName) )
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	try
	{
		function["Calibration"].as<std::string>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"Calibration\" in function \"GetECALCorrectedEnergy\" not found for calculation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                         << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"Calibration\" in function \"GetECALCorrectedEnergy\" should be type of double for calculation of variables " << quantityNames[0] << "\" and " << "\""
		                                                                                                                                      << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	std::string *Calibration = new std::string();
	(*Calibration)           = function["Calibration"].as<std::string>();

	std::map<int, std::pair<double,double>> correction;
	std::ifstream                           configFile;
	configFile.open(*Calibration);
	if( !configFile.is_open() )
	{
		return 0;
	}
	int runNumber;
	double a, b;
	while( configFile >> runNumber >> a >> b )
	{
		correction[runNumber] = std::pair<double,double>(a,b);
	}
	if( not configFile.eof() )
	{
		std::cerr << "ERROR: Invalid input correction for run " << runNumber << std::endl;
		return 0;
	}
	configFile.close();

	try
	{
		function["RangeECAL1"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"RangeECAL1\" in function \"GetECALCorrectedEnergy\" not found for caluclation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                        << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"RangeECAL1\" in function \"GetECALCorrectedEnergy\" should be type of double for calculation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                                       << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	double *RangeECAL1 = new double();
	(*RangeECAL1)      = function["RangeECAL1"].as<double>();

	return (new antok::user::cdreis::functions::GetECALCorrectedEnergy( Energy,
	                                                                    ClusterZ,
	                                                                    RangeECAL1,
	                                                                    RunNumber,
	                                                                    correction,
	                                                                    data.getAddr<std::vector<double>>(quantityName)));
};

antok::Function * antok::user::cdreis::generateGetECALCorrectedTiming( const YAML::Node&         function,
                                                                       std::vector<std::string>& quantityNames,
                                                                       int                       index )
{
	if( quantityNames.size() > 1 )
	{
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}
	std::string quantityName = quantityNames[0];

	std::vector<std::pair<std::string, std::string> > args;

	args.push_back(std::pair<std::string, std::string>( "Timing"  , "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "Energy"  , "std::vector<double>" ));
	args.push_back(std::pair<std::string, std::string>( "ClusterZ", "std::vector<double>" ));

	if( not antok::generators::functionArgumentHandler(args, function, index) )
	{
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<double>* Timing   = data.getAddr<std::vector<double>>(args[0].first);
	std::vector<double>* Energy   = data.getAddr<std::vector<double>>(args[1].first);
	std::vector<double>* ClusterZ = data.getAddr<std::vector<double>>(args[2].first);

	if( not data.insert<std::vector<double>>(quantityName) )
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(quantityNames);
		return 0;
	}

	try
	{
		function["Calibration"].as<std::string>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"Calibration\" in function \"GetECALCorrectedTiming\" not found for calculation of variable \"" << quantityName << "\"" << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"Calibration\" in function \"GetECALCorrectedTiming\" should be type of string for calculation of variable \"" << quantityName << "\"" << std::endl;
		return 0;
	}
	std::string* Calibration = new std::string();
	(*Calibration)           = function["Calibration"].as<std::string>();

	std::map<std::string, std::vector<double>> correctionValues;
	std::ifstream configFile;
	configFile.open(*Calibration);
	if( !configFile.is_open() )
	{
		return 0;
	}
	std::string name;
	double a, b, c, d, e, f, g;
	while( configFile >> name >> a >> b >> c >> d >> e >> f >> g )
	{
		std::vector<double> values;
		values.resize(7, double() );
		values[0] = a;
		values[1] = b;
		values[2] = c;
		values[3] = d;
		values[4] = e;
		values[5] = f;
		values[6] = g;

		correctionValues[name] = values;
	}
	if( not configFile.eof() )
	{
		std::cerr << "ERROR: Invalid input correction for " << name << std::endl;
		return 0;
	}
	configFile.close();

	try
	{
		function["RangeECAL1"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"RangeECAL1\" in function \"GetECALCorrectedTiming\" not found for calculation of variable \"" << quantityName << "\"" << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<std::string> &exception )
	{
		std::cerr << "Argument \"Calibration\" in function \"GetECALCorrectedTiming\" should be type of double for calculation of variable \"" << quantityName << "\"" << std::endl;
		return 0;
	}
	double *RangeECAL1 = new double();
	(*RangeECAL1)      = function["RangeECAL1"].as<double>();

	return (new antok::user::cdreis::functions::GetECALCorrectedTiming( Timing,
	                                                                    Energy,
	                                                                    ClusterZ,
	                                                                    RangeECAL1,
	                                                                    correctionValues,
	                                                                    data.getAddr<std::vector<double>>(quantityName)));
};

antok::Function *antok::user::cdreis::generateGetPhotonPairParticles( const YAML::Node&         function,
                                                                      std::vector<std::string>& quantityNames,
                                                                      int                       index )
{
	if (quantityNames.size() > 2)
	{
		std::cerr << "Too many names for function \"" << function["Name"] << "\"." << std::endl;
		return 0;
	}
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>( "Photons"        , "std::vector<TLorentzVector>" ));
	args.push_back(std::pair<std::string, std::string>( "ECALIndex"      , "std::vector<int>"            ));


	if( not antok::generators::functionArgumentHandler(args, function, index) )
	{
		std::cerr << antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data &data = antok::ObjectManager::instance()->getData();

	std::vector<TLorentzVector> *Photons         = data.getAddr<std::vector<TLorentzVector> >(args[0].first);
	std::vector<int>            *ECALIndex       = data.getAddr<std::vector<int> >           (args[1].first);

	std::string resultParticles    = quantityNames[0];
	std::string resultHasParticles = quantityNames[1];

	if( not data.insert<std::vector<TLorentzVector> >(resultParticles) )
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultParticles);
		return 0;
	}
	if( not data.insert<int>(resultHasParticles) )
	{
		std::cerr << antok::Data::getVariableInsertionErrorMsg(resultHasParticles);
		return 0;
	}

	try
	{
		function["ECALResolution"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ECALResolution\" in function \"GetPhotonPairParticles\" not found for calculation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                            << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<double> &exception )
	{
		std::cerr << "Argument \"ECALResolution\" in function \"GetPhotonPairParticles\" should be type of double  for calculation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                                            << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	double* ECALResolution = new double();
	(*ECALResolution)      = function["ECALResolution"].as<double>();

	try
	{
		function["ECAL1Resolution"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ECAL1Resolution\" in function \"GetPhotonPairParticles\" not found for calculation of variables " << quantityNames[0] << "\" and " << "\""
		                                                                                                                           << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<double> &exception )
	{
		std::cerr << "Argument \"ECAL1Resolution\" in function \"GetPhotonPairParticles\" should be type of double  for calculation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                                             << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	double* ECAL1Resolution = new double();
	(*ECAL1Resolution)      = function["ECAL1Resolution"].as<double>();

	try
	{
		function["ECAL2Resolution"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"ECAL2Resolution\" in function \"GetPhotonPairParticles\" not found for calculation of variables  \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                              << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<double> &exception )
	{
		std::cerr << "Argument \"ECAL2Resolution\" in function \"GetPhotonPairParticles\" should be type of double  for calculation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                                             << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	double* ECAL2Resolution = new double();
	(*ECAL2Resolution)      = function["ECAL2Resolution"].as<double>();

	try
	{
		function["Mass"].as<double>();
	}
	catch( const YAML::InvalidNode &exception )
	{
		std::cerr << "Argument \"Mass\" in function \"GetPhotonPairParticles\" not found for calculation of variables  \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                   << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	catch( const YAML::TypedBadConversion<double> &exception )
	{
		std::cerr << "Argument \"Mass\" in function \"GetPhotonPairParticles\" should be type of double  for calculation of variables \"" << quantityNames[0] << "\" and " << "\""
		                                                                                                                                  << quantityNames[1] << "\""      << std::endl;
		return 0;
	}
	double* Mass = new double();
	(*Mass)      = function["Mass"].as<double>();

	return (new antok::user::cdreis::functions::GetPhotonPairParticles( Photons,
	                                                                    Mass,
	                                                                    ECALResolution,
	                                                                    ECAL1Resolution,
	                                                                    ECAL2Resolution,
	                                                                    ECALIndex,
	                                                                    data.getAddr<std::vector<TLorentzVector> >(resultParticles),
	                                                                    data.getAddr<int>(resultHasParticles)));
};
