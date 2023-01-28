#ifndef ANTOK_USER_jbeckers_H
#define ANTOK_USER_jbeckers_H

#include<yaml-cpp/yaml.h>

#include "generators_functions.h"

namespace antok {

	class Function;

	namespace user {

		namespace jbeckers {

			// reuse macro from generators_functions.h
			FUNCTION_PROTOTYPE(getUserFunction);

			FUNCTION_PROTOTYPE(generateScale);

			FUNCTION_PROTOTYPE(generateIncrease);

			FUNCTION_PROTOTYPE(generateCompareIndices);

			FUNCTION_PROTOTYPE(generatePrintEvent);

			FUNCTION_PROTOTYPE(generateCalculateCollinearityAngle);

		}

	}

}

#endif  // ANTOK_USER_jbeckers_H
