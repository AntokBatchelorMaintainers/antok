#ifndef ANTOK_CUTTER_H
#define ANTOK_CUTTER_H

#include<map>
#include<string>
#include<vector>
#include<bitset>

class TTree;

namespace antok {

	using bitmask = std::bitset<64>;

	class Cut;

	class Cutter {

		friend class Initializer;

	  public:

		static Cutter* instance();

		bool cut();

		const bitmask& getCutPattern() const { return _cutPattern; };

		bool fillOutTrees() const;

		bitmask getCutmaskForNames(std::vector<std::string> names) const;

		bitmask getAllCutsCutmaskForCutTrain(std::string cutTrainName) const;
		const std::vector<antok::Cut*>& getCutsForCutTrain(std::string cutTrainName) const;

		bool cutOnInCutmask(bitmask mask, const antok::Cut* cut) const;

		std::string getAbbreviations(bitmask cutPattern, std::string cutTrainName) const;

		const std::map<std::string, std::vector<bitmask> >& getWaterfallCutmasks();
		const std::map<std::string, std::vector<bitmask> >& getCutmasksAllCutsOffSeparately();
		const std::map<std::string, std::vector<bitmask> >& getCutmasksAllCutsOnSeparately();

		bool cutInCutTrain(std::string cutName, std::string cutTrainName) const;

		const bool* getCutResult(antok::Cut*) const;

	  private:

		Cutter()
			: _cutPattern(0) { };

		static Cutter* _cutter;

		bitmask _cutPattern;

		std::map<std::string, std::map<std::string, antok::Cut*> > _cutTrainsMap;
		std::map<std::string, std::vector<antok::Cut*> > _cutTrainsCutOrderMap;
		std::map<std::string, antok::Cut*> _cutsMap;
		std::map<std::string, TTree*> _outTreeMap;

		std::map<std::string, std::vector<bitmask> > _waterfallCutmasksCache;
		std::map<std::string, std::vector<bitmask> > _singleOffCutmasksCache;
		std::map<std::string, std::vector<bitmask> > _singleOnCutmasksCache;

		std::vector<std::pair<TTree*, bitmask> > _treesToFill;

		std::vector<std::pair<antok::Cut*, bool*> > _cuts;

	};

}

#endif

