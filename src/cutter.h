#ifndef ANTOK_CUTTER_H
#define ANTOK_CUTTER_H

#include<map>
#include<string>
#include<vector>

class TTree;

namespace antok {

	class Cut;

	class Cutter {

		friend class Initializer;

	  public:

		static Cutter* instance();

		bool cut();

		const long& getCutPattern() const { return _cutPattern; };

		bool fillOutTrees() const;

		long getCutmaskForNames(std::vector<std::string> names) const;

		long getAllCutsCutmaskForCutTrain(std::string cutTrainName) const;
		const std::vector<antok::Cut*>& getCutsForCutTrain(std::string cutTrainName) const;

		bool cutOnInCutmask(long mask, const antok::Cut* cut) const;

		std::string getAbbreviations(long cutPattern, std::string cutTrainName) const;

		const std::map<std::string, std::vector<long> >& getWaterfallCutmasks();
		const std::map<std::string, std::vector<long> >& getCutmasksAllCutsOffSeparately();
		const std::map<std::string, std::vector<long> >& getCutmasksAllCutsOnSeparately();

		bool cutInCutTrain(std::string cutName, std::string cutTrainName) const;

		const bool* getCutResult(antok::Cut*) const;

	  private:

		Cutter() { };

		static Cutter* _cutter;

		long _cutPattern;

		std::map<std::string, std::map<std::string, antok::Cut*> > _cutTrainsMap;
		std::map<std::string, std::vector<antok::Cut*> > _cutTrainsCutOrderMap;
		std::map<std::string, antok::Cut*> _cutsMap;
		std::map<std::string, TTree*> _outTreeMap;

		std::map<std::string, std::vector<long> > _waterfallCutmasksCache;
		std::map<std::string, std::vector<long> > _singleOffCutmasksCache;
		std::map<std::string, std::vector<long> > _singleOnCutmasksCache;

		std::vector<std::pair<TTree*, long> > _treesToFill;

		std::vector<std::pair<antok::Cut*, bool*> > _cuts;

	};

}

#endif

