#ifndef ANTOK_GENERATORS_PLOTS_H
#define ANTOK_GENERATORS_PLOTS_H

#include<yaml-cpp/yaml.h>

namespace antok {

	class Plot;

	namespace plotUtils {

		class GlobalPlotOptions;

	}

	namespace generators {

		antok::Plot* generate1DPlot(const YAML::Node& plot, antok::plotUtils::GlobalPlotOptions plotOptions);
		antok::Plot* generate2DPlot(const YAML::Node& plot, antok::plotUtils::GlobalPlotOptions plotOptions);

	}

}

#endif

