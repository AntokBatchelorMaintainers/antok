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

			FUNCTION_PROTOTYPE(generateGetSumOverVector);
			FUNCTION_PROTOTYPE(generateGetVector3VectorAttributes);
			FUNCTION_PROTOTYPE(generateGetVectorLorentzVectorAttributes);
			FUNCTION_PROTOTYPE(generateGetSumLorentzVectors);
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
			FUNCTION_PROTOTYPE(generateGetOmegaDalitzVariables);
			FUNCTION_PROTOTYPE(generateGetPionLVs);
			FUNCTION_PROTOTYPE(generateGetTwoPionCombinationLV);
			FUNCTION_PROTOTYPE(generateGetThreePionCombinationLV);
			FUNCTION_PROTOTYPE(generateGetFourPionCombinationLV);
			FUNCTION_PROTOTYPE(generateGetResolutions);
			FUNCTION_PROTOTYPE(generateGetPi0Resolutions);
			FUNCTION_PROTOTYPE(generateGetSelectedPhotonLVs);
			FUNCTION_PROTOTYPE(generateGetNeutralMeson);
			FUNCTION_PROTOTYPE(generateGetPiPiNeutralSystem);
			FUNCTION_PROTOTYPE(generateGetAngles3P);

			bool registerOutputVarTypes(antok::Data& data, const std::vector<std::string>& quantityNames, const std::vector<std::string>& outputVarTypes);

		}


		// stream operator for std::vector
		template<typename T>
		inline
		std::ostream&
		operator << (std::ostream&         out,
		             const std::vector<T>& vec)
		{
			const size_t vecSize = vec.size();
			if (vecSize == 0) {
				return out << "[]";
			}
			out << "[";
			for (unsigned int i = 0; i < (vecSize - 1); ++i)
				out << vec[i] << ", ";
			return out << vec[vecSize - 1] << "]";
		}

	}

}

#endif  // ANTOK_USER_CDREIS_H
