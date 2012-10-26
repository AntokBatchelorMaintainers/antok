#include<cutter.h>

#include<assert.h>
#include<sstream>

#include<TTree.h>

antok::Cutter* antok::Cutter::_cutter = 0;

antok::Cutter* antok::Cutter::instance() {
	if(_cutter == 0) {
		_cutter = new antok::Cutter();
//		_cutter->_statsHist = 0;
	}
	return _cutter;
}

antok::Cutter::Cutter() {

/*													// No cuts: 511
	_cuts.push_back(new TrigMask(0x1));				// 1	510
	_cuts.push_back(new VrtxZ(-29, -66));			// 2	508
	_cuts.push_back(new VrtxR(1.55));				// 4	504
	_cuts.push_back(new nRPDTracks(1));				// 8	496
	_cuts.push_back(new RPDProtMass(0.2));			// 16	480
	_cuts.push_back(new CedarKaon());				// 32	448
	_cuts.push_back(new TPrime(0.1, 1.));			// 64	384
	_cuts.push_back(new RPDPlanarity());			// 128	256
	_cuts.push_back(new Exclusivity(191.5, 3.78));	// 256	0
*/
}

bool antok::Cutter::cut() {

	bool success = true;
	_cutPattern = 0;
	for(unsigned int i = 0; i < _cuts.size(); ++i) {
		success = success and (*(_cuts[i].first))();
		bool result = (*(_cuts[i].second));
		if(result) {
			_cutPattern += (1<<i);
		}
	}
	return success;

/*	int cutmask = 0;
	bool cut_previously = false;
	for(unsigned int i = 0; i < _cuts.size(); ++i) {
		if((*(_cuts.at(i)))()) {
			cutmask += (1<<i);
			cut_previously = true;
		} else if (!cut_previously) {
			_statsHist->Fill(((_cuts.at(i))->get_longname()).c_str(), 1);
		}
	}
	return cutmask;
*/
};

long antok::Cutter::getAllCutsCutmaskForCutTrain(std::string cutTrainName) const {

	std::map<std::string, std::map<std::string, antok::Cut*> >::const_iterator cutTrainsMap_it = _cutTrainsMap.find(cutTrainName);
	assert(cutTrainsMap_it != _cutTrainsMap.end());
	const std::map<std::string, antok::Cut*>& cuts = cutTrainsMap_it->second;
	std::vector<std::string> names;
	for(std::map<std::string, antok::Cut*>::const_iterator cuts_it = cuts.begin(); cuts_it != cuts.end(); ++cuts_it) {
		names.push_back(cuts_it->first);
	}
	return getCutmaskForNames(names);

}

const std::vector<antok::Cut*>& antok::Cutter::getCutsForCutTrain(std::string cutTrainName) const {

	std::map<std::string, std::vector<antok::Cut*> >::const_iterator cutTrainsMap_it = _cutTrainsCutOrderMap.find(cutTrainName);
	assert(cutTrainsMap_it != _cutTrainsCutOrderMap.end());
	return cutTrainsMap_it->second;

};

bool antok::Cutter::cutOnInCutmask(long mask, const antok::Cut* cut) const {

	for(unsigned int i = 0; i < _cuts.size(); ++i) {
		const antok::Cut* innerCut = _cuts[i].first;
		if(cut == innerCut) {
			return ((mask>>i)&1);
		}
	}
	assert(false);

}

long antok::Cutter::getCutmaskForNames(std::vector<std::string> names) const {

	long cutmask = 0;
	for(unsigned int i = 0; i < names.size(); ++i) {
		std::map<std::string, antok::Cut*>::const_iterator cutsMap_it = _cutsMap.find(names[i]);
		assert(cutsMap_it != _cutsMap.end());
		const antok::Cut* cut = cutsMap_it->second;
		unsigned int index = 0;
		for( ; _cuts[index].first != cut && index < _cuts.size(); ++index);
		assert(index < _cuts.size());
		cutmask += 1<<index;
	}
	return cutmask;

};

std::string antok::Cutter::getAbbreviations(long cutPattern, std::string cutTrainName) const {

	const std::vector<antok::Cut*>& cuts = getCutsForCutTrain(cutTrainName);

	std::stringstream strStr;
	strStr<<"(";

	bool noCuts = true;
	bool first = true;
	for(unsigned int i = 0; i < cuts.size(); ++i) {
		if(cutOnInCutmask(cutPattern, cuts[i])) {
			noCuts = false;
			if(first) {
				first = false;
			} else {
				strStr<<"|";
			}
			strStr<<cuts[i]->getAbbreviation();
		}
	}
	if(noCuts) {
		return "(NoCuts)";
	}
	strStr<<")";
	return strStr.str();

/*
	unsigned int size = _cuts.size();
	assert(bitmask>>(size) == 0);
	std::stringstream sstr;
	sstr<<"(";
	if(bitmask == 0) {
		sstr<<"AllCuts";
	} else if (bitmask == ((1<<size)-1)) {
		sstr<<"NoCuts";
	} else {
		bool first = true;
		for(unsigned int i = 0; i < size; ++i) {
			if(!((bitmask>>i)&0x1)) {
				if(first) {
					first = false;
				} else {
					sstr<<"|";
				}
				sstr<<_cuts.at(i)->get_abbreviation();
			}
		}
	}
	sstr<<")";
	std::string retval = sstr.str();
	return retval;
*/
};

const std::map<std::string, std::vector<long> >& antok::Cutter::getWaterfallCutmasks() {

	if(_waterfallCutmasksCache.empty()) {

		for(std::map<std::string, std::vector<antok::Cut*> >::const_iterator cutTrainsCutOrder_it = _cutTrainsCutOrderMap.begin();
		    cutTrainsCutOrder_it != _cutTrainsCutOrderMap.end();
		    ++cutTrainsCutOrder_it)
		{

			const std::string& cutTrainName = cutTrainsCutOrder_it->first;
			const std::vector<antok::Cut*> cuts = cutTrainsCutOrder_it->second;

			std::vector<std::string> cutNames;
			_waterfallCutmasksCache[cutTrainName].push_back(0);

			for(unsigned int i = 0; i < cuts.size(); ++i) {
				cutNames.push_back(cuts[i]->getShortName());
				_waterfallCutmasksCache[cutTrainName].push_back(getCutmaskForNames(cutNames));
			}

		}

	}

	return _waterfallCutmasksCache;

}

const std::map<std::string, std::vector<long> >& antok::Cutter::getCutmasksAllCutsOffSeparately() {

	if(_singleOffCutmasksCache.empty()) {

		for(std::map<std::string, std::vector<antok::Cut*> >::const_iterator cutTrainsCutOrder_it = _cutTrainsCutOrderMap.begin();
		    cutTrainsCutOrder_it != _cutTrainsCutOrderMap.end();
		    ++cutTrainsCutOrder_it)
		{
			const std::string& cutTrainName = cutTrainsCutOrder_it->first;
			long cutmaskTemplate = getAllCutsCutmaskForCutTrain(cutTrainName);
			for(unsigned int i = 0; i < _cuts.size(); ++i) {
				if((cutmaskTemplate>>i)&1) {
					long cutmask = cutmaskTemplate - (1<<i);
					_singleOffCutmasksCache[cutTrainName].push_back(cutmask);
				}
			}
		}
	}

	return _singleOffCutmasksCache;

}

const std::map<std::string, std::vector<long> >& antok::Cutter::getCutmasksAllCutsOnSeparately() {

	if(_singleOnCutmasksCache.empty()) {

		for(std::map<std::string, std::vector<antok::Cut*> >::const_iterator cutTrainsCutOrder_it = _cutTrainsCutOrderMap.begin();
		    cutTrainsCutOrder_it != _cutTrainsCutOrderMap.end();
		    ++cutTrainsCutOrder_it)
		{
			const std::string& cutTrainName = cutTrainsCutOrder_it->first;
			long cutmaskTemplate = getAllCutsCutmaskForCutTrain(cutTrainName);
			for(unsigned int i = 0; i < _cuts.size(); ++i) {
				if((cutmaskTemplate>>i)&1) {
					long cutmask = (1<<i);
					_singleOnCutmasksCache[cutTrainName].push_back(cutmask);
				}
			}
		}

	}

	return _singleOnCutmasksCache;

}

bool antok::Cutter::cutInCutTrain(std::string cutName, std::string cutTrainName) const {

	std::map<std::string, std::map<std::string, antok::Cut*> >::const_iterator cutTrainsMap_it = _cutTrainsMap.find(cutTrainName);
	assert(cutTrainsMap_it != _cutTrainsMap.end());
	std::map<std::string, antok::Cut*>::const_iterator cuts_it = cutTrainsMap_it->second.find(cutName);
	return (cuts_it != cutTrainsMap_it->second.end());

}

bool antok::Cutter::fillOutTrees() const {

	bool success = true;
	for(unsigned int i = 0; i < _treesToFill.size(); ++i) {
		TTree* tree = _treesToFill[i].first;
		long mask = _treesToFill[i].second;
		if((mask&_cutPattern) == mask) {
			success = success and (tree->Fill() > 0);
		}
	}
	return success;

}

/*
bool antok::Cutter::set_stats_histogram(TH1D* stats) {

	if(_statsHist != 0) {
		return false;
	}
	_statsHist = stats;
	return true;

}

antok::Cutter::~Cutter() {

	for(unsigned int i = 0; i < _cuts.size(); ++i) {
		delete _cuts.at(i);
	}

};
*/
