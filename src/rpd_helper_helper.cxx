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
	leftSideLower.push_back(-0.05928);
	leftSideLower.push_back(-0.0204924);
	leftSideLower.push_back(-0.0791312);
	antok::FitParameters leftSideLowerParams = antok::FitParameters(leftSideLower);
	std::vector<double> leftSideUpper;
	leftSideUpper.push_back(0.0642057);
	leftSideUpper.push_back(-0.00171496);
	leftSideUpper.push_back(-0.0122408);
	antok::FitParameters leftSideUpperParams = antok::FitParameters(leftSideUpper);

	_leftSlabLimitsParameters = std::pair<antok::FitParameters,
	                                      antok::FitParameters>(leftSideLowerParams, leftSideUpperParams);

	std::vector<double> middleLower;
	middleLower.push_back(-0.130212);
	middleLower.push_back(-0.00156496);
	middleLower.push_back(-0.0124321);
	antok::FitParameters middleLowerParams = antok::FitParameters(middleLower);
	std::vector<double> middleUpper;
	middleUpper.push_back(0.129598);
	middleUpper.push_back(0.00164405);
	middleUpper.push_back(-0.0124098);
	antok::FitParameters middleUpperParams = antok::FitParameters(middleUpper);

	_middleSlabLimitsParameters = std::pair<antok::FitParameters,
	                                        antok::FitParameters>(middleLowerParams, middleUpperParams);

	std::vector<double> rightSideLower;
	rightSideLower.push_back(-0.0642652);
	rightSideLower.push_back(0.00158913);
	rightSideLower.push_back(-0.011811);
	antok::FitParameters rightSideLowerParams = antok::FitParameters(rightSideLower);
	std::vector<double> rightSideUpper;
	rightSideUpper.push_back(0.0587655);
	rightSideUpper.push_back(0.0207241);
	rightSideUpper.push_back(-0.0791871);
	antok::FitParameters rightSideUpperParams = antok::FitParameters(rightSideUpper);

	_rightSlabLimitsParameters = std::pair<antok::FitParameters,
	                                       antok::FitParameters>(rightSideLowerParams, rightSideUpperParams);


}

std::pair<double, double> antok::RpdHelperHelper::getLimits(const double& rpdPhi,
                                                            const TVector2& vertexXY)
{

	std::pair<double, TVector2> transformedVariabes = rotateRpdAndVertex(rpdPhi, vertexXY);
	const double& correctedProtonPhi = transformedVariabes.first;
	const TVector2& correctedVertexXY = transformedVariabes.second;

	std::pair<double, double> retval = std::pair<double, double>(0., 0.);

	const double RPD_SLAB_STEP = 0.025;

	unsigned int slabType = 0;
	for(; correctedProtonPhi > RPD_SLAB_PHI_ANGLES[slabType] + RPD_SLAB_STEP; ++slabType);
	switch(slabType) {
		case LEFTSLAB:
			retval.first = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _leftSlabLimitsParameters.first);
			retval.second = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _leftSlabLimitsParameters.second);
			break;
		case MIDDLESLAB:
			retval.first = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _middleSlabLimitsParameters.first);
			retval.second = evaluatePlaneFunction(correctedVertexXY.X(), correctedVertexXY.Y(), _middleSlabLimitsParameters.second);
			break;
		case RIGHTSLAB:
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
