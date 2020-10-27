
#include<kbicker.h>

#include<TLorentzRotation.h>

#include<constants.h>
#include<data.h>
#include<functions.hpp>
#include<generators_functions.h>
#include<kbicker_functions.hpp>
#include<rpd_helper_helper.h>
#include<yaml_utils.hpp>

antok::Function* antok::user::kbicker::getUserFunction(const YAML::Node&               function,
                                                       const std::vector<std::string>& quantityNames,
                                                       int                             index)
{
	std::string functionName = antok::YAMLUtils::getString(function["Name"]);
	Function* antokFunctionPtr = 0;
	if(functionName == "getRpdExpectedHitsParameters") {
		antokFunctionPtr = antok::user::kbicker::generateGetRpdExpectedHitsParameters(function, quantityNames, index);
	} else if(functionName == "getRpdPhi") {
		antokFunctionPtr = antok::user::kbicker::generateGetRpdPhi(function, quantityNames, index);
	} else if(functionName == "handleExtraTracks") {
		antokFunctionPtr = antok::user::kbicker::generateGetCutOnExtraTracks(function, quantityNames, index);
	}
	return antokFunctionPtr;
}

antok::Function* antok::user::kbicker::generateGetRpdExpectedHitsParameters(const YAML::Node& function, const std::vector<std::string>& quantityNames, int index) {

	using antok::YAMLUtils::hasNodeKey;

	if(quantityNames.size() != 4) {
		std::cerr<<"Need 4 names for function \""<<function["Name"]<<"\"."<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("Vertex", "TVector3"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* beamLorentzVecAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* XLorentzVecAddr = data.getAddr<TLorentzVector>(args[1].first);
	TVector3* vertexAddr = data.getAddr<TVector3>(args[2].first);

	std::vector<std::string> possiblyConstArgs;
	possiblyConstArgs.push_back("XOffset");
	possiblyConstArgs.push_back("YOffset");
	possiblyConstArgs.push_back("XAngle");
	possiblyConstArgs.push_back("YAngle");

	for(unsigned int i = 0; i < possiblyConstArgs.size(); ++i) {
		if(not hasNodeKey(function, possiblyConstArgs[i])) {
			std::cerr<<"Argument \""<<possiblyConstArgs[i]<<"\" not found (required for function \""<<function["Name"]<<"\")."<<std::endl;
			return 0;
		}
	}

	double* xOffsetAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[0]]);
	double* yOffsetAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[1]]);
	double* xAngleAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[2]]);
	double* yAngleAddr = antok::YAMLUtils::getAddress<double>(function[possiblyConstArgs[3]]);

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	return (new antok::user::kbicker::functions::GetRPDExpectedHitsParameters(beamLorentzVecAddr,
	                                                                          XLorentzVecAddr,
	                                                                          vertexAddr,
	                                                                          xOffsetAddr,
	                                                                          yOffsetAddr,
	                                                                          xAngleAddr,
	                                                                          yAngleAddr,
	                                                                          quantityAddrs[0],
	                                                                          quantityAddrs[1],
	                                                                          quantityAddrs[2],
	                                                                          quantityAddrs[3]));

}

antok::Function* antok::user::kbicker::generateGetRpdPhi(const YAML::Node& function, const std::vector<std::string>& quantityNames, int index)
{

	if(not (quantityNames.size() == 2 or quantityNames.size() == 4)) {
		std::cerr<<"Need 2 names for function \"getRpdPhi\""<<std::endl;
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::string method = antok::YAMLUtils::getString(function["Method"]);
	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("BeamLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("RPDProtonLorentzVec", "TLorentzVector"));
	args.push_back(std::pair<std::string, std::string>("XLorentzVec", "TLorentzVector"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	TLorentzVector* beamLVAddr = data.getAddr<TLorentzVector>(args[0].first);
	TLorentzVector* RPDprotLVAddr = data.getAddr<TLorentzVector>(args[1].first);
	TLorentzVector* xLVAddr = data.getAddr<TLorentzVector>(args[2].first);

	TVector3* vertexAddr = 0;

	if(method == "Prediction") {
		args.clear();
		args.push_back(std::pair<std::string, std::string>("Vertex", "TVector3"));
		if(not antok::generators::functionArgumentHandler(args, function, index)) {
			std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
			std::cerr<<"\"Prediction\" method requires \"Vertex\"."<<std::endl;
			return 0;
		}
		vertexAddr = data.getAddr<TVector3>(args[0].first);
	}


	int methodSwitch = -1;
	if(method == "Projection") {
		if(quantityNames.size() != 2) {
			std::cerr<<"\"getRpdPhi\" with method \"Projection\" only works with 2 calculated quantities, found "
			         <<quantityNames.size()<<" when calculating variables \"["<<std::endl;
			for(unsigned int i = 0; i < quantityNames.size()-1; ++i) {
				std::cerr<<quantityNames[i]<<", ";
			}
			std::cerr<<quantityNames[quantityNames.size()-1]<<"]\"."<<std::endl;
			return 0;
		}
		methodSwitch = 0;
	} else if (method == "Rotation") {
		methodSwitch = 1;
	} else if (method == "Prediction") {
		methodSwitch = 2;
	} else {
		if(method == "") {
			std::cerr<<"Could not convert \"Method\" to std::string in function \"getRpdPhi\" when calculating variables \"["<<std::endl;
		} else {
			std::cerr<<"Method \""<<method<<"\" is not supported for function \"getRpdPhi\" when calculating variables \"["<<std::endl;
		}
		for(unsigned int i = 0; i < quantityNames.size()-1; ++i) {
			std::cerr<<quantityNames[i]<<", ";
		}
		std::cerr<<quantityNames[quantityNames.size()-1]<<"]\"."<<std::endl;
		return 0;
	}

	std::vector<double*> quantityAddrs;
	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(not data.insert<double>(quantityNames[i])) {
			std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
			return 0;
		}
		quantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
	}

	double* quantityAddrs2 = 0;
	double* quantityAddrs3 = 0;

	if(quantityNames.size() == 4) {
		quantityAddrs2 = quantityAddrs[2];
		quantityAddrs3 = quantityAddrs[3];
	}
	return (new antok::user::kbicker::functions::GetRpdPhi(beamLVAddr,
	                                                       RPDprotLVAddr,
	                                                       xLVAddr,
	                                                       quantityAddrs[0],
	                                                       quantityAddrs[1],
	                                                       methodSwitch,
	                                                       quantityAddrs2,
	                                                       quantityAddrs3,
	                                                       vertexAddr));

};

antok::Function* antok::user::kbicker::generateGetCutOnExtraTracks(const YAML::Node& function, const std::vector<std::string>& quantityNames, int index)
{

	if(quantityNames.size() != 3) {
		std::cerr<<"Need 2 names for function \"handleExtraTracks\""<<std::endl;
		return 0;
	}

	std::vector<std::pair<std::string, std::string> > args;
	args.push_back(std::pair<std::string, std::string>("TrackTimes", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("TrackTimeSigmas", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("TrackNHits", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("TrackZFirst", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("TrackZLast", "std::vector<double>"));
	args.push_back(std::pair<std::string, std::string>("TrackQP", "std::vector<double>"));

	if(not antok::generators::functionArgumentHandler(args, function, index)) {
		std::cerr<<antok::generators::getFunctionArgumentHandlerErrorMsg(quantityNames);
		return 0;
	}

	antok::Data& data = antok::ObjectManager::instance()->getData();

	std::vector<double>* trackTimesAddr = data.getAddr<std::vector<double> >(args[0].first);
	std::vector<double>* trackTimeSigmasAddr = data.getAddr<std::vector<double> >(args[1].first);
	std::vector<double>* trackNHitsAddr = data.getAddr<std::vector<double> >(args[2].first);
	std::vector<double>* trackZFirstAddr = data.getAddr<std::vector<double> >(args[3].first);
	std::vector<double>* trackZLastAddr = data.getAddr<std::vector<double> >(args[4].first);
	std::vector<double>* trackQAddr = data.getAddr<std::vector<double> >(args[5].first);

	std::vector<double*> doubleQuantityAddrs;
	std::vector<int*> intQuantityAddrs;

	const unsigned int NUMBER_OF_INT_QUANTITIES = 1;

	for(unsigned int i = 0; i < quantityNames.size(); ++i) {
		if(i >= quantityNames.size() - NUMBER_OF_INT_QUANTITIES) {
			if(not data.insert<int>(quantityNames[i])) {
				std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
				return 0;
			}
			intQuantityAddrs.push_back(data.getAddr<int>(quantityNames[i]));
		} else {
			if(not data.insert<double>(quantityNames[i])) {
				std::cerr<<antok::Data::getVariableInsertionErrorMsg(quantityNames, quantityNames[i]);
				return 0;
			}
			doubleQuantityAddrs.push_back(data.getAddr<double>(quantityNames[i]));
		}
	}

	return (new antok::user::kbicker::functions::GetCutOnExtraTracks(trackTimesAddr,
	                                                                 trackTimeSigmasAddr,
	                                                                 trackNHitsAddr,
	                                                                 trackZFirstAddr,
	                                                                 trackZLastAddr,
	                                                                 trackQAddr,
	                                                                 doubleQuantityAddrs[0],
	                                                                 intQuantityAddrs[0],
	                                                                 doubleQuantityAddrs[1]));

}


void antok::user::kbicker::getRPDDeltaPhiResProjection(const TLorentzVector& pBeam,
                                                       const TLorentzVector& pProton,
                                                       const TLorentzVector& pX,
                                                       double& delta_phi, double& res)
{

	double rpd_Phi = pProton.Phi() / TMath::TwoPi();
	if(rpd_Phi < 0.) rpd_Phi += 1.;
	int hit_slat = (int)((rpd_Phi*24.) + 0.5);
	if(hit_slat >= 24) hit_slat -= 24;
	if(hit_slat % 2) res = 3.75;
	else res = 7.5;
	res = (res/180.)*TMath::Pi();
	res = std::sqrt(res*res + 0.067260*0.067260);

	TVector3 proton3 = pProton.Vect();
	TVector3 pTot3 = pX.Vect();
	TVector3 pBeam3 = pBeam.Vect();
	TVector3 a = pTot3 - (pTot3.Dot(pBeam3.Unit()))*pBeam3.Unit();
	TVector3 b = proton3 - (proton3.Dot(pBeam3.Unit()))*pBeam3.Unit();
	TVector3 n = a.Cross(b);
	if(n.Dot(TVector3(0.,0.,1.)) > 0) {
		n *= -1.;
	}
	double y = ((n.Cross(a)).Dot(b)) / (((a.Cross(b)).Cross(a)).Mag()*b.Mag());
	double x = (a.Unit()).Dot(b.Unit());
	delta_phi = TMath::ATan2(y, x)-TMath::Pi();
	if(delta_phi < (-1.*TMath::Pi())) {
		delta_phi += TMath::TwoPi();
	}

};

void antok::user::kbicker::getRPDDeltaPhiResRotation(const TLorentzVector& pBeam,
                                                     const TLorentzVector& pProton,
                                                     const TLorentzVector& pX,
                                                     double& delta_phi, double& res,
                                                     double& phiProton, double& phiX)
{

	double rpd_Phi = pProton.Phi();
	rpd_Phi /= TMath::TwoPi();

	if(rpd_Phi < 0.) rpd_Phi += 1.;
	int hit_slat = (int)((rpd_Phi*24.) + 0.5);
	if(hit_slat >= 24) hit_slat -= 24;
	if(hit_slat % 2) res = 3.75;
	else res = 7.5;
	res = (res/180.)*TMath::Pi();

	TVector3 z(0.,0.,1.);
	double angle = pBeam.Vect().Angle(z);
	TVector3 ra = pBeam.Vect().Cross(z);
	TRotation r;
	r.Rotate(angle, ra);
	TLorentzRotation Tr(r);

	TLorentzVector proton_lab(pProton);
	TLorentzVector pX_lab(pX);
	proton_lab.Transform(Tr);
	pX_lab.Transform(Tr);

	phiProton = proton_lab.Phi();
	phiX = pX_lab.Phi();

	delta_phi = std::fabs(phiX - phiProton) - TMath::Pi();
	res = std::sqrt(res*res + 0.067260*0.067260);

};

void antok::user::kbicker::getRPDDeltaPhiResRotation(const TLorentzVector& pBeam,
                                                     const TLorentzVector& pProton,
                                                     const TLorentzVector& pX,
                                                     double& delta_phi, double& res)
{
	double phiProton, phiX;
	getRPDDeltaPhiResRotation(pBeam, pProton, pX, delta_phi, res, phiProton, phiX);
}

void antok::user::kbicker::getRPDDeltaPhiResPrediction(const TLorentzVector& pBeam,
                                                       const TLorentzVector& pProton,
                                                       const TLorentzVector& pX,
                                                       const TVector3& vertex,
                                                       double& delta_phi, double& likelihood,
                                                       double& phiProton, double& phiX)
{

	const TVector3& beamVector = pBeam.Vect();
	const TVector3& XVector = pX.Vect();

	TVector3 proton = beamVector - XVector;

	phiX = proton.Phi();
	phiProton = pProton.Phi();

	delta_phi = phiX - phiProton;

	if(delta_phi < -TMath::Pi()) {
		delta_phi += TMath::TwoPi();
	}
	if(delta_phi > TMath::Pi()) {
		delta_phi -= TMath::TwoPi();
	}

	double vertexX = vertex.X();
	double vertexY = vertex.Y();
	likelihood = RpdHelperHelper::getInstance()->getLikelihood(delta_phi,
	                                                           phiProton,
	                                                           vertexX,
	                                                           vertexY,
	                                                           pX.M(),
	                                                           proton.Mag());

};

void antok::user::kbicker::getRPDDeltaPhiResPrediction(const TLorentzVector& pBeam,
                                                     const TLorentzVector& pProton,
                                                     const TLorentzVector& pX,
                                                     const TVector3& vertex,
                                                     double& delta_phi, double& likelihood)
{
	double phiProton, phiX;
	getRPDDeltaPhiResPrediction(pBeam, pProton, pX, vertex, delta_phi, likelihood, phiProton, phiX);
}

void antok::user::kbicker::getRPDExpectedHitsParameters(const TLorentzVector& pBeam,
                                                        const TLorentzVector& pX,
                                                        const TVector3& vertex,
                                                        const double& xOffset,
                                                        const double& yOffset,
                                                        const double& xAngle,
                                                        const double& yAngle,
                                                        double& rpdPhiRingA,
                                                        double& rpdPhiRingB,
                                                        double& rpdZRingA,
                                                        double& rpdZRingB)
{

	static const double INNER_RING_RADIUS = antok::Constants::RPDRingARadius();
	static const double OUTER_RING_RADIUS = antok::Constants::RPDRingBRadius();
	static const double RPD_CENTER_Z_POSITION = -36.2; // cm, from COMGEANT's geom_hadron_2008.ffr

	const TVector3& beamVector = pBeam.Vect();
	const TVector3& XVector = pX.Vect();

	TVector3 proton = beamVector - XVector;

	TRotation rotation;
	rotation.SetToIdentity();
	rotation.RotateX(yAngle);
	rotation.RotateY(xAngle);
	rotation.Invert();

	TVector3 transformedVertex(vertex);
	const TVector3 translation(xOffset, yOffset, RPD_CENTER_Z_POSITION);
	transformedVertex -= translation;
	transformedVertex.Transform(rotation);

	proton.Transform(rotation);

	bool error = false;
	const double a = proton.X()*proton.X() + proton.Y()*proton.Y();
	const double b = 2 * (proton.X()*transformedVertex.X() + proton.Y()*transformedVertex.Y());
	double c = transformedVertex.X()*transformedVertex.X() + transformedVertex.Y()*transformedVertex.Y() - INNER_RING_RADIUS*INNER_RING_RADIUS;
	double n;
	try {
		n = antok::utils::getPositiveSolutionOfQuadraticEquation(a, b, c);
	} catch (const int& number) {
		if(number != 42) {
			std::cerr<<"In antok::getRPDExpectedHitsParameters: Something went very "
			         <<"wrong when calculating the RPD hit coordinates for ring A. Aborting..."<<std::endl;
			throw;
		}
		error = true;
	}
	TVector3 coordinatesRingA = transformedVertex + n*proton;
	rpdZRingA = coordinatesRingA.Z();
	rpdPhiRingA = coordinatesRingA.Phi();
	if(error) {
		rpdZRingA = 0.;
		rpdPhiRingA = -10.;
		error = false;
	}
	c = transformedVertex.X()*transformedVertex.X() + transformedVertex.Y()*transformedVertex.Y() - OUTER_RING_RADIUS*OUTER_RING_RADIUS;
	try {
		n = antok::utils::getPositiveSolutionOfQuadraticEquation(a, b, c);
	} catch (const int& number) {
		if(number != 42) {
			std::cerr<<"In antok::getRPDExpectedHitsParameters: Something went very "
			         <<"wrong when calculating the RPD hit coordinates for ring B. Aborting..."<<std::endl;
			throw;
		}
		error = true;
	}
	TVector3 coordinatesRingB = transformedVertex + n*proton;
	rpdZRingB = coordinatesRingB.Z();
	rpdPhiRingB = coordinatesRingB.Phi();
	if(error) {
		rpdZRingB = 0.;
		rpdPhiRingB = -10.;
	}

}

bool antok::user::kbicker::extraTracksCut(std::vector<double> trackTimes,
                                          std::vector<double> trackTimeSigmas,
                                          std::vector<double> trackNHits,
                                          std::vector<double> trackZFirst,
                                          std::vector<double> trackZLast,
                                          std::vector<double> trackQP)
{

	return false;

}
