
#include<assert.h>
#include<sstream>

#include<cutter.h>

antok::Cutter* antok::Cutter::_cutter = 0;

antok::Cutter* antok::Cutter::instance() {
	if(_cutter == 0) {
		_cutter = new antok::Cutter();
		_cutter->_statsHist = 0;
	}
	return _cutter;
}

antok::Cutter::Cutter() {

													// No cuts: 511
	_cuts.push_back(new TrigMask(0x1));				// 1	510
	_cuts.push_back(new VrtxZ(-29, -66));			// 2	508
	_cuts.push_back(new VrtxR(1.55));				// 4	504
	_cuts.push_back(new nRPDTracks(1));				// 8	496
	_cuts.push_back(new RPDProtMass(0.2));			// 16	480
	_cuts.push_back(new CedarKaon());				// 32	448
	_cuts.push_back(new TPrime(0.1, 1.));			// 64	384
	_cuts.push_back(new RPDPlanarity());			// 128	256
	_cuts.push_back(new Exclusivity(191.5, 3.78));	// 256	0

}

int antok::Cutter::get_cutmask(const antok::Event& event) {

	int cutmask = 0;
	bool cut_previously = false;
	for(unsigned int i = 0; i < _cuts.size(); ++i) {
		if((*(_cuts.at(i)))(event)) {
			cutmask += (1<<i);
			cut_previously = true;
		} else if (!cut_previously) {
			_statsHist->Fill(((_cuts.at(i))->get_longname()).c_str(), 1);
		}
	}
	return cutmask;

};

std::string antok::Cutter::get_abbreviations(int bitmask) {

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

};

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
