#ifndef ANTOK_USER_CDREIS_H
#define ANTOK_USER_CDREIS_H

#include <yaml-cpp/yaml.h>

#include "data.h"
#include "generators_functions.h"

namespace antok {

	class Function;

	namespace user {

		namespace cdreis {

			// reuse macro from generators_functions.h
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

			bool registerOutputVarTypes(antok::Data& data, const std::vector<std::string>& quantityNames, const std::vector<std::string>& outputVarTypes);

		}

	}

}

#endif  // ANTOK_USER_CDREIS_H
