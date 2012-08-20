
#include<cutter.h>

hlib::Cutter::Cutter() {

	_cuts.push_back(new TrigMask(0x1));
	_cuts.push_back(new VrtxZ(-28.4, -68.4));
	_cuts.push_back(new VrtxR(1.75));
	_cuts.push_back(new nRPDTracks(1));
	_cuts.push_back(new RPDProtMass(0.2));
	_cuts.push_back(new CedarKaon());
	_cuts.push_back(new TPrime(0.1, -1.));
	_cuts.push_back(new RPDPlanarity());
	_cuts.push_back(new Exclusivity(191., 3.28));

}

int hlib::Cutter::get_cutmask(const hlib::Event& event) {

	int cutmask = 0;
	for(unsigned int i = 0; i < _cuts.size(); ++i) {
		if((*(_cuts.at(i)))(event)) {
			cutmask += (1<<i);
		}
	}
	return cutmask;

}
