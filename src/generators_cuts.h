#ifndef ANTOK_GENERATORS_CUTS_H
#define ANTOK_GENERATORS_CUTS_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Cut;

	namespace generators {

		antok::Cut* generateRangeCut(const YAML::Node& cut,
		                             const std::string& shortName,
		                             const std::string& longName,
		                             const std::string& abbreviation);
		antok::Cut* generateEqualityCut(const YAML::Node& cut,
		                                const std::string& shortName,
		                                const std::string& longName,
		                                const std::string& abbreviation);

	}

}

#endif

