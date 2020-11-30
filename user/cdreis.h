#ifndef ANTOK_USER_CDREIS_H
#define ANTOK_USER_CDREIS_H

#include <yaml-cpp/yaml.h>

#include "data.h"

namespace antok {

	class Function;

	namespace user {

		namespace cdreis {

			// macro for shorter prototyping
			#define FUNCTION_PROTOTYPE(fctName) antok::Function* fctName (const YAML::Node& function, const std::vector<std::string>& quantityNames, const int index)
			FUNCTION_PROTOTYPE(getUserFunction);

			FUNCTION_PROTOTYPE(generateGetVector3VectorAttributes);
			FUNCTION_PROTOTYPE(generateGetVectorLorentzVectorAttributes);
			FUNCTION_PROTOTYPE(generateGetNominalMassDifferences);
			FUNCTION_PROTOTYPE(generateGetRecoilLorentzVec);
			FUNCTION_PROTOTYPE(generateGetECALCorrectedEnergy);
			FUNCTION_PROTOTYPE(generateGetECALCorrectedTiming);
			FUNCTION_PROTOTYPE(generateGetCleanedEcalClusters);
			FUNCTION_PROTOTYPE(generateGetECALVariables);
			FUNCTION_PROTOTYPE(generateGetPhotonLorentzVecs);
			FUNCTION_PROTOTYPE(generateGetPhotonPairParticles);
			FUNCTION_PROTOTYPE(generateGetPi0Pair);
			FUNCTION_PROTOTYPE(generateGetKinematicFittingMass);
			FUNCTION_PROTOTYPE(generateGetOmega);
			FUNCTION_PROTOTYPE(generateGetFittedOmegaMassVsPrecisionGoal);
			FUNCTION_PROTOTYPE(generateGetThreePionCombinationMass);

			const bool registerOutputVarTypes (antok::Data& data, const std::vector<std::string>& quantityNames, const std::vector<std::string>& outputVarTypes);
		}

	}

}

#endif  // ANTOK_USER_CDREIS_H
