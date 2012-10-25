#ifndef ANTOK_PLOTTER_H
#define ANTOK_PLOTTER_H

#include<map>
#include<string>
#include<vector>

namespace YAML {
	class Node;
}

namespace antok {

	class Plot;

	class Plotter {

		friend class Initializer;

	  public:

		static Plotter* instance();

		void fill(long cutPattern);

		static bool handleAdditionalCuts(const YAML::Node& cuts, std::map<std::string, std::vector<long> >& map);

	  private:

		Plotter();

		static Plotter* _plotter;

		std::vector<antok::Plot*> _plots;

	};

	namespace plotUtils {

		struct GlobalPlotOptions {

			GlobalPlotOptions(const YAML::Node& optionNode);
			bool plotsForSequentialCuts;
			bool plotsWithSingleCutsOn;
			bool plotsWithSingleCutsOff;
			std::string statisticsHistInName;
			std::string statisticsHistOutName;
			std::map<std::string, std::vector<long> > cutMasks;

		  private:

			bool handleOnOffOption(std::string optionName, const YAML::Node& option, std::string location) const;

		};

	}

}

#endif
