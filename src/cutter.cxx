
#include<cutter.h>

hlib::Cutter::Cutter() {

	_cuts.push_back(new TrigMask());
	_cuts.push_back(new VrtxZ());
	_cuts.push_back(new VrtxR(1.75));

}
