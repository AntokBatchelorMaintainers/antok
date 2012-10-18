#ifndef ANTOK_GENERATORS_CUTS_H
#define ANTOK_GENERATORS_CUTS_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Cut;

	namespace generators {

		bool generateCut(const YAML::Node& cutEntry,
		                 antok::Cut*& antokCut,
		                 bool*& result);

	}

}

#endif

