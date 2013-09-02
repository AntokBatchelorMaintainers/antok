#include<rpd_helper_helper.h>

#include<TF1.h>

#include<constants.h>

#include<TApplication.h>
#include<cmath>

antok::RpdHelperHelper* antok::RpdHelperHelper::_instance = 0;

antok::RpdHelperHelper* antok::RpdHelperHelper::getInstance() {

	if(not _instance) {
		_instance = new antok::RpdHelperHelper();
	}
	return _instance;

}

double antok::RpdHelperHelper::getLikelihood(const double& deltaPhi,
                                             const double& rpdProtonPhi,
                                             const double& vertexX,
                                             const double& vertexY,
                                             const double& xMass,
                                             const double& predictedProtonMom) const
{
	const double inBetweenThreshold = 0.1;
	bool inBetweenSlabs = (std::fabs(std::fabs(rpdProtonPhi)-0.261799388) < inBetweenThreshold) or   //  15 deg
	                      (std::fabs(std::fabs(rpdProtonPhi)-0.785398163) < inBetweenThreshold) or   //  45 deg
	                      (std::fabs(std::fabs(rpdProtonPhi)-1.30899694) < inBetweenThreshold) or    //  75 deg
	                      (std::fabs(std::fabs(rpdProtonPhi)-1.83259571) < inBetweenThreshold) or    // 105 deg
	                      (std::fabs(std::fabs(rpdProtonPhi)-2.35619449019) < inBetweenThreshold) or // 135 deg
	                      (std::fabs(std::fabs(rpdProtonPhi)-2.87979327) < inBetweenThreshold);      // 165 deg
	if(inBetweenSlabs) {
		return -10.;
	}
	TVector2 vertexXY(vertexX, vertexY);
	std::pair<double, TVector2> transformedVariabes = rotateRpdAndVertex(rpdProtonPhi, vertexXY);
	__correctedProtonPhi = transformedVariabes.first;
	__correctedVertexXY = transformedVariabes.second;

	std::vector<double> limits = getAllLimits();
	unsigned int slabType = getHitSlab(__correctedProtonPhi);

	bool print = false;
	if(slabType > 2) {
		std::cerr<<"getHitSlab failed, returned slabType="<<slabType<<". Setting likelihood to 0."<<std::endl;
		return 0;
	}

	double limitFirst, limitSecond;
	switch(slabType) {
		case LEFTSLAB:
			limitFirst = limits[0];
			limitSecond = limits[1];
			break;
		case MIDDLESLAB:
			limitFirst = limits[2];
			limitSecond = limits[3];
			break;
		case RIGHTSLAB:
			limitFirst = limits[4];
			limitSecond = limits[5];
			break;
	}
	if((limitSecond - limitFirst) < 0.) {
		return 0.;
	}
	std::pair<double, double> sharpSigmas = getSharpSigmas(xMass, predictedProtonMom);
	if(sharpSigmas.first < 1e-9) {
		sharpSigmas.first = 1e-9;
	}
	if(sharpSigmas.second < 1e-9) {
		sharpSigmas.second = 1e-9;
	}
	double broadSigmaScalingFactor = getBroadSigmasScalingFactor(xMass, predictedProtonMom);
	if(broadSigmaScalingFactor < 0.) {
		broadSigmaScalingFactor = 0.00000001;
	}
	std::pair<double, double> broadSigmas = getBroadSigmas(xMass, predictedProtonMom);

	_likelihoodFunction->SetParameter(0, limitSecond); // limit right
	_likelihoodFunction->SetParameter(1, -limitFirst); // limit left
	_likelihoodFunction->SetParameter(2, sharpSigmas.second); // sharp sigma right
	_likelihoodFunction->SetParameter(3, sharpSigmas.first); // sharp sigma left
	_likelihoodFunction->SetParameter(4, broadSigmaScalingFactor); // relative contribution broad sigmas
	_likelihoodFunction->SetParameter(5, broadSigmas.second); // broad sigma right
	_likelihoodFunction->SetParameter(6, broadSigmas.first); // broad sigma left
	double likelihood = getScalingFactor(broadSigmaScalingFactor) * _likelihoodFunction->Eval(deltaPhi);

	if(likelihood < 0. or likelihood > 1.) {
		std::cout<<"-------------------------------------------------------------"<<std::endl;
		std::cout<<"rpdProtonPhi-165 = "<<std::fabs(std::fabs(rpdProtonPhi)-2.87979327)<<std::endl;
		std::cout<<"rpdProtonPhi-135 = "<<std::fabs(std::fabs(rpdProtonPhi)-2.35619449019)<<std::endl;
		std::cout<<"rpdProtonPhi-105 = "<<std::fabs(std::fabs(rpdProtonPhi)-1.83259571)<<std::endl;
		std::cout<<"rpdProtonPhi-75  = "<<std::fabs(std::fabs(rpdProtonPhi)-1.30899694)<<std::endl;
		std::cout<<"rpdProtonPhi-45  = "<<std::fabs(std::fabs(rpdProtonPhi)-0.785398163397)<<std::endl;
		std::cout<<"rpdProtonPhi-15  = "<<std::fabs(std::fabs(rpdProtonPhi)-0.261799388)<<std::endl;
		std::cout<<"rpdProtonPhi     = "<<rpdProtonPhi<<std::endl;
		std::cout<<"corRpdProtonPhi  = "<<__correctedProtonPhi<<std::endl;
		std::cout<<"slabType         = "<<slabType<<std::endl;
		std::cout<<"deltaPhi         = "<<deltaPhi<<std::endl;
		std::cout<<std::endl;
		std::cout<<"allLimits        = ["<<limits[0];
		for(unsigned int i = 1; i <  limits.size(); ++i) std::cout<<", "<<limits[i];
		std::cout<<"]"<<std::endl;
		std::cout<<std::endl;
		std::cout<<"vertex x         = "<<vertexX<<std::endl;
		std::cout<<"vertex y         = "<<vertexY<<std::endl;
		std::cout<<"m                = "<<xMass<<std::endl;
		std::cout<<"q                = "<<predictedProtonMom<<std::endl;
		std::cout<<std::endl;
		std::cout<<"limit left       = "<<limitFirst<<std::endl;
		std::cout<<"limit right      = "<<limitSecond<<std::endl;
		std::cout<<"ssigma left      = "<<sharpSigmas.first<<std::endl;
		std::cout<<"ssigma right     = "<<sharpSigmas.second<<std::endl;
		std::cout<<"bsigma scaling   = "<<broadSigmaScalingFactor<<std::endl;
		std::cout<<"bsigma left      = "<<broadSigmas.first<<std::endl;
		std::cout<<"bsigma right     = "<<broadSigmas.second<<std::endl;
		std::cout<<std::endl;
		std::cout<<"unscaled LH      = "<<_likelihoodFunction->Eval(deltaPhi)<<std::endl;
		std::cout<<"scaling factor   = "<<getScalingFactor(broadSigmaScalingFactor)<<std::endl;
		std::cout<<"likelihood       = "<<likelihood<<std::endl;
		std::cout<<"-------------------------------------------------------------"<<std::endl;
		if(print) return 0.;
	}
	if(likelihood < 0.) {
		likelihood = 0.;
	}
	return likelihood;
}

double antok::RpdHelperHelper::getScalingFactor(const double& xMass,
                                                const double& predictedProtonMom) const
{
	double broadSigmaScalingFactor = getBroadSigmasScalingFactor(xMass, predictedProtonMom);
	return getScalingFactor(broadSigmaScalingFactor);
}

double antok::RpdHelperHelper::getScalingFactor(const double& broadSigmasScalingFactor) const
{
	double sum = 2 * (1 + broadSigmasScalingFactor);
	return (1. / sum);
}

antok::RpdHelperHelper::RpdHelperHelper()
	: _likelihoodFunction(new TF1("antok_RpdHelperHelper_likelihoodFunction", "-1*(TMath::Erf((x-[0])/[2])-TMath::Erf((x+[1])/[3])+[4]*(TMath::Erf((x-[0])/[5])-TMath::Erf((x+[1])/[6])))")),
	  __correctedProtonPhi(),
	  __correctedVertexXY()
{

	RPD_SLAB_PHI_ANGLES.push_back(-0.196);
	RPD_SLAB_PHI_ANGLES.push_back(0.002);
	RPD_SLAB_PHI_ANGLES.push_back(0.194);

	std::vector<double> leftSideLower;
	leftSideLower.push_back(-0.05928);
	leftSideLower.push_back(-0.0204924);
	leftSideLower.push_back(-0.0791312);
	antok::FitParameters leftSideLowerParams(leftSideLower);
	std::vector<double> leftSideUpper;
	leftSideUpper.push_back(0.0642057);
	leftSideUpper.push_back(-0.00171496);
	leftSideUpper.push_back(-0.0122408);
	antok::FitParameters leftSideUpperParams(leftSideUpper);

	_leftSlabLimitsParameters = std::pair<antok::FitParameters,
	                                      antok::FitParameters>(leftSideLowerParams, leftSideUpperParams);

	std::vector<double> middleLower;
	middleLower.push_back(-0.130212);
	middleLower.push_back(-0.00156496);
	middleLower.push_back(-0.0124321);
	antok::FitParameters middleLowerParams(middleLower);
	std::vector<double> middleUpper;
	middleUpper.push_back(0.129598);
	middleUpper.push_back(0.00164405);
	middleUpper.push_back(-0.0124098);
	antok::FitParameters middleUpperParams(middleUpper);

	_middleSlabLimitsParameters = std::pair<antok::FitParameters,
	                                        antok::FitParameters>(middleLowerParams, middleUpperParams);

	std::vector<double> rightSideLower;
	rightSideLower.push_back(-0.0642652);
	rightSideLower.push_back(0.00158913);
	rightSideLower.push_back(-0.011811);
	antok::FitParameters rightSideLowerParams(rightSideLower);
	std::vector<double> rightSideUpper;
	rightSideUpper.push_back(0.0587655);
	rightSideUpper.push_back(0.0207241);
	rightSideUpper.push_back(-0.0791871);
	antok::FitParameters rightSideUpperParams(rightSideUpper);

	_rightSlabLimitsParameters = std::pair<antok::FitParameters,
	                                       antok::FitParameters>(rightSideLowerParams, rightSideUpperParams);

	std::vector<double> sharpSigmasLeft;
	sharpSigmasLeft.push_back(-0.0051137);
	sharpSigmasLeft.push_back(0.00091844);
	sharpSigmasLeft.push_back(0.0169681);
	antok::FitParameters sharpSigmasLeftParams(sharpSigmasLeft);

	std::vector<double> sharpSigmasRight;
	sharpSigmasRight.push_back(-0.00570525);
	sharpSigmasRight.push_back(0.000937381);
	sharpSigmasRight.push_back(0.0171201);
	antok::FitParameters sharpSigmasRightParams(sharpSigmasRight);

	_sharpSigmaParameters = std::pair<antok::FitParameters,
	                                  antok::FitParameters>(sharpSigmasLeftParams, sharpSigmasRightParams);

	std::vector<double> broadSigmasLeft;
	broadSigmasLeft.push_back(0.307881);
	broadSigmasLeft.push_back(-0.0189389);
	broadSigmasLeft.push_back(1.99952);
	broadSigmasLeft.push_back(-0.0123104);
	broadSigmasLeft.push_back(-2.88051);
	broadSigmasLeft.push_back(0.0494617);
	antok::FitParameters broadSigmasLeftParams(broadSigmasLeft);

	_broadSigmaParameters = std::pair<double, double>(0.215095, 0.214565);

	std::vector<double> broadSigmaContribution;
	broadSigmaContribution.push_back(-0.0337792);
	broadSigmaContribution.push_back(0.00578239);
	broadSigmaContribution.push_back(0.0325389);

	_broadSigmaScalingFactor = antok::FitParameters(broadSigmaContribution);

}

std::pair<double, double> antok::RpdHelperHelper::getLimits(const double& rpdPhi,
                                                            const TVector2& vertexXY) const
{

	std::pair<double, TVector2> transformedVariabes = rotateRpdAndVertex(rpdPhi, vertexXY);
	const double& correctedProtonPhi = transformedVariabes.first;
	const TVector2& correctedVertexXY = transformedVariabes.second;

	std::pair<double, double> retval = std::pair<double, double>(0., 0.);

	unsigned int slabType = getHitSlab(correctedProtonPhi);
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

std::pair<double, double> antok::RpdHelperHelper::getSharpSigmas(const double& m,
                                                                 const double& q) const
{
	std::pair<double, double> retval = std::pair<double, double>(0., 0.);
	retval.first = evaluateSharpSigmaFunction(m, q, _sharpSigmaParameters.first);
	retval.second = evaluateSharpSigmaFunction(m, q, _sharpSigmaParameters.second);
	return retval;
}

std::pair<double, double> antok::RpdHelperHelper::getBroadSigmas(const double& m,
                                                                 const double& q) const
{
	return _broadSigmaParameters;
}

double antok::RpdHelperHelper::getBroadSigmasScalingFactor(const double& m,
                                                           const double& q) const
{
	return evaluateSharpSigmaFunction(m, q, _broadSigmaScalingFactor);
}

std::vector<double> antok::RpdHelperHelper::getAllLimits() const {
	std::vector<double> retval(6, 0.);
	retval[0] = evaluatePlaneFunction(__correctedVertexXY.X(), __correctedVertexXY.Y(), _leftSlabLimitsParameters.first);
	retval[1] = evaluatePlaneFunction(__correctedVertexXY.X(), __correctedVertexXY.Y(), _leftSlabLimitsParameters.second);
	retval[2] = evaluatePlaneFunction(__correctedVertexXY.X(), __correctedVertexXY.Y(), _middleSlabLimitsParameters.first);
	retval[3] = evaluatePlaneFunction(__correctedVertexXY.X(), __correctedVertexXY.Y(), _middleSlabLimitsParameters.second);
	retval[4] = evaluatePlaneFunction(__correctedVertexXY.X(), __correctedVertexXY.Y(), _rightSlabLimitsParameters.first);
	retval[5] = evaluatePlaneFunction(__correctedVertexXY.X(), __correctedVertexXY.Y(), _rightSlabLimitsParameters.second);
	return retval;
}

double antok::RpdHelperHelper::evaluatePlaneFunction(const double& x,
                                                     const double& y,
                                                     const antok::FitParameters& fitParameters) const
{
	return fitParameters.getParameter(0) +
	       fitParameters.getParameter(1) * x +
	       fitParameters.getParameter(2) * y;
}

double antok::RpdHelperHelper::evaluateSharpSigmaFunction(const double& m,
                                                          const double& q,
                                                          const antok::FitParameters& fitParameters) const
{
	return fitParameters.getParameter(0) +
	       fitParameters.getParameter(1) * m*m +
	       fitParameters.getParameter(2) / q;
}

double antok::RpdHelperHelper::evaluateBroadSigmaFunction(const double& m,
                                                          const double& q,
                                                          const antok::FitParameters& fitParameters) const
{

//	[0]+[1]*(y-[2])^2+[3]*(x-[4])^2+[5]*x*y
	double fitParam2 = fitParameters.getParameter(2);
	double fitParam4 = fitParameters.getParameter(4);
	return fitParameters.getParameter(0) +
	       fitParameters.getParameter(1) * (m-fitParam2)*(m-fitParam2) +
	       fitParameters.getParameter(3) * (q-fitParam4)*(q-fitParam4) +
	       fitParameters.getParameter(5) * m * q;
}

std::pair<double, TVector2> antok::RpdHelperHelper::rotateRpdAndVertex(const double& rpdPhi,
                                                                       const TVector2 vertexXY) const
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

unsigned int antok::RpdHelperHelper::getHitSlab(const double& correctedProtonPhi) const
{
	const double RPD_SLAB_STEP = 0.025;
	unsigned int slabType = 0;
	for(; correctedProtonPhi > RPD_SLAB_PHI_ANGLES[slabType] + RPD_SLAB_STEP; ++slabType);
	return slabType;
}
