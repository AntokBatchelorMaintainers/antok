#ifndef ANTOK_CUTTER_H
#define ANTOK_CUTTER_H

#include<map>
#include<string>
#include<vector>

#include<TH1D.h>

#include<cut.hpp>
#include<event.h>

class TTree;

namespace antok {

	class Cutter {

		friend class Initializer;

	  public:

		static Cutter* instance();

		bool cut();

		long getCutmaskForNames(std::vector<std::string> names) const;

		const long& getCutPattern() const { return _cutPattern; };

		long getAllCutsCutmaskForCutTrain(std::string cutTrainName) const;
		const std::vector<antok::Cut*>& getCutsForCutTrain(std::string cutTrainName) const;

		bool cutOnInCutmask(long mask, const antok::Cut* cut) const;

		std::string getAbbreviations(long cutPattern, std::string cutTrainName) const;

		const std::map<std::string, std::vector<long> >& getWaterfallCutmasks();



//		bool set_stats_histogram(TH1D* stats);
//		TH1D* get_stats_histogram() { return _statsHist; };

	  private:

		Cutter();
//		~Cutter();

		static Cutter* _cutter;

		long _cutPattern;

		std::map<std::string, std::map<std::string, antok::Cut*> > _cutTrainsMap;
		std::map<std::string, std::vector<antok::Cut*> > _cutTrainsCutOrderMap;
		std::map<std::string, antok::Cut*> _cutsMap;
		std::map<std::string, TTree*> _outTreeMap;

		std::map<std::string, std::vector<long> > _waterfallCutmasksCache;

		std::vector<std::pair<antok::Cut*, bool*> > _cuts;

//		TH1D* _statsHist;

	};

}

#endif

