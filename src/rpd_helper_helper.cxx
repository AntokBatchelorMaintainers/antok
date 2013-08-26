#include<rpd_helper_helper.h>

#include<iostream>

#include<TVector2.h>

#include<constants.h>

antok::RpdHelperHelper* antok::RpdHelperHelper::_instance = 0;

antok::RpdHelperHelper* antok::RpdHelperHelper::getInstance() {

	if(not _instance) {
		_instance = new antok::RpdHelperHelper();
	}
	return _instance;

}

antok::RpdHelperHelper::RpdHelperHelper() {

	RPD_SLAB_PHI_ANGLES.push_back(-0.196);
	RPD_SLAB_PHI_ANGLES.push_back(0.002);
	RPD_SLAB_PHI_ANGLES.push_back(0.194);

	std::vector<double> leftSideLower;
	leftSideLower.push_back(-0.0600894);
	leftSideLower.push_back(-0.0211335);
	leftSideLower.push_back(-0.0789747);
	antok::FitParameters leftSideLowerParams = antok::FitParameters(leftSideLower);
	std::vector<double> leftSideUpper;
	leftSideUpper.push_back(0.0640313);
	leftSideUpper.push_back(-0.00169554);
	leftSideUpper.push_back(-0.0122071);
	antok::FitParameters leftSideUpperParams = antok::FitParameters(leftSideUpper);

	_leftSlabLimitsParameters = std::pair<antok::FitParameters,
	                                      antok::FitParameters>(leftSideLowerParams, leftSideUpperParams);

	std::vector<double> middleLower;
	middleLower.push_back(-0.131045);
	middleLower.push_back(-0.0016428);
	middleLower.push_back(-0.0122572);
	antok::FitParameters middleLowerParams = antok::FitParameters(middleLower);
	std::vector<double> middleUpper;
	middleUpper.push_back(0.128913);
	middleUpper.push_back(0.00157375);
	middleUpper.push_back(-0.0125813);
	antok::FitParameters middleUpperParams = antok::FitParameters(middleUpper);

	_middleSlabLimitsParameters = std::pair<antok::FitParameters,
	                                        antok::FitParameters>(middleLowerParams, middleUpperParams);

	std::vector<double> rightSideLower;
	rightSideLower.push_back(-0.0646319);
	rightSideLower.push_back(0.00171409);
	rightSideLower.push_back(-0.0124187);
	antok::FitParameters rightSideLowerParams = antok::FitParameters(rightSideLower);
	std::vector<double> rightSideUpper;
	rightSideUpper.push_back(0.0587703);
	rightSideUpper.push_back(0.0213982);
	rightSideUpper.push_back(-0.0799849);
	antok::FitParameters rightSideUpperParams = antok::FitParameters(rightSideUpper);

	_rightSlabLimitsParameters = std::pair<antok::FitParameters,
	                                       antok::FitParameters>(rightSideLowerParams, rightSideUpperParams);


}

std::pair<double, double> antok::RpdHelperHelper::getLimits(const double& rpdPhi,
                                                            const TVector2& vertexXY)
{

	std::cout<<"before transform: phi="<<rpdPhi<<" x="<<vertexXY.X()<<" y="<<vertexXY.Y()<<std::endl;

	std::pair<double, TVector2> transformedVariabes = rotateRpdAndVertex(rpdPhi, vertexXY);
	const double& correctedProtonPhi = transformedVariabes.first;
	const TVector2& correctedVertexXY = transformedVariabes.second;

	std::cout<<"after transform: phi="<<correctedProtonPhi<<" x="<<correctedVertexXY.X()<<" y="<<correctedVertexXY.Y()<<std::endl;

	std::pair<double, double> retval = std::pair<double, double>(0., 0.);

	const double RPD_SLAB_STEP = 0.025;

	unsigned int slabType = 0;
	for(; correctedProtonPhi > RPD_SLAB_PHI_ANGLES[slabType] + RPD_SLAB_STEP; ++slabType);
	switch(slabType) {
		case LEFTSLAB:
			std::cout<<"returning stuff for left slab"<<std::endl;
			retval.first = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _leftSlabLimitsParameters.first);
			retval.second = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _leftSlabLimitsParameters.second);
			break;
		case MIDDLESLAB:
			std::cout<<"returning stuff for middle slab"<<std::endl;
			retval.first = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _middleSlabLimitsParameters.first);
			retval.second = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _middleSlabLimitsParameters.second);
			break;
		case RIGHTSLAB:
			std::cout<<"returning stuff for right slab"<<std::endl;
			retval.first = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _rightSlabLimitsParameters.first);
			retval.second = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _rightSlabLimitsParameters.second);
			break;
		default:
			throw;
	}

	return retval;

}

double antok::RpdHelperHelper::evaluatePlaneFunction(const double& x,
                                                     const double& y,
                                                     const antok::FitParameters& fitParameters)
{
	return fitParameters.getParameter(0) +
	       fitParameters.getParameter(1) * x +
	       fitParameters.getParameter(2) * y;
}

std::pair<double, TVector2> antok::RpdHelperHelper::rotateRpdAndVertex(const double& rpdPhi,
                                                                       const TVector2 vertexXY)
{
	const double rpdPeriod = antok::Constants::RPDPeriod();

	TVector2 correctedVertexXY = vertexXY;
	double correctedProtonPhi = rpdPhi;
	double correctedVertexPhi = correctedVertexXY.Phi();

	while(correctedProtonPhi > 0.25) {
		correctedProtonPhi -= rpdPeriod;
		correctedVertexPhi -= rpdPeriod;

	}
	while(correctedProtonPhi < -0.25) {
		correctedProtonPhi += rpdPeriod;
		correctedVertexPhi += rpdPeriod;
	}

	correctedVertexXY.SetMagPhi(vertexXY.Mod(), correctedVertexPhi);

	return std::pair<double, TVector2>(correctedProtonPhi, correctedVertexXY);

}
