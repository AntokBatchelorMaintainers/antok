
#include<assert.h>
#include<sstream>

#include<cutter.h>

hlib::Cutter* hlib::Cutter::_cutter = NULL;

hlib::Cutter* hlib::Cutter::instance() {
	if(_cutter == NULL) {
		_cutter = new hlib::Cutter();
	}
	return _cutter;
}

hlib::Cutter::Cutter() {

													// No cuts: 511
	_cuts.push_back(new TrigMask(0x1));				// 1	510
	_cuts.push_back(new VrtxZ(-28.4, -68.4));		// 2	508
	_cuts.push_back(new VrtxR(1.75));				// 4	504
	_cuts.push_back(new nRPDTracks(1));				// 8	496
	_cuts.push_back(new RPDProtMass(0.2));			// 16	480
	_cuts.push_back(new CedarKaon());				// 32	448
	_cuts.push_back(new TPrime(0.1, -1.));			// 64	384
	_cuts.push_back(new RPDPlanarity());			// 128	256
	_cuts.push_back(new Exclusivity(191., 3.28));	// 256	0

}

int hlib::Cutter::get_cutmask(const hlib::Event& event) {

	int cutmask = 0;
	for(unsigned int i = 0; i < _cuts.size(); ++i) {
		if((*(_cuts.at(i)))(event)) {
			cutmask += (1<<i);
		}
	}
	return cutmask;

};

std::string hlib::Cutter::get_abbreviations(int bitmask) {

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

hlib::Cutter::~Cutter() {

	for(unsigned int i = 0; i < _cuts.size(); ++i) {
		delete _cuts.at(i);
	}

};
