#include<cutter.h>

#include<assert.h>
#include<sstream>

#include<TTree.h>

#include<cut.hpp>

antok::Cutter* antok::Cutter::_cutter = 0;

antok::Cutter* antok::Cutter::instance() {
	if(_cutter == 0) {
		_cutter = new antok::Cutter();
	}
	return _cutter;
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
			const std::vector<antok::Cut*>& cuts = cutTrainsCutOrder_it->second;
			for(unsigned int i = 0; i < cuts.size(); ++i) {
				for(unsigned int j = 0; j < _cuts.size(); ++j) {
					if(cuts[i] == _cuts[j].first) {
						long cutmask = cutmaskTemplate - (1<<j);
						_singleOffCutmasksCache[cutTrainName].push_back(cutmask);
						break;
					}
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
			const std::vector<antok::Cut*>& cuts = cutTrainsCutOrder_it->second;
			for(unsigned int i = 0; i < cuts.size(); ++i) {
				for(unsigned int j = 0; j < _cuts.size(); ++j) {
					if(cuts[i] == _cuts[j].first) {
						long cutmask = 1<<j;
						_singleOnCutmasksCache[cutTrainName].push_back(cutmask);
						break;
					}
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

const bool* antok::Cutter::getCutResult(antok::Cut* cut) const {

	for(unsigned int i = 0; i < _cuts.size(); ++i) {
		if(cut == _cuts[i].first) {
			return _cuts[i].second;
		}
	}
	assert(false);

}

