#ifndef ANTOK_GENERATORS_CUTS_H
#define ANTOK_GENERATORS_CUTS_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Cut;

	namespace generators {

		antok::Cut* generateCut(const YAML::Node& cut);

	}

}

#endif

