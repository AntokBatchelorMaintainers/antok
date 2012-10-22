#include<event.h>

#include<basic_calcs.h>
#include<constants.h>
#include<functions.hpp>

antok::Event* antok::Event::_event = 0;

antok::Event* antok::Event::instance() {
	if(_event == 0) {
		_event = new antok::Event();
	}
	return _event;
}

bool antok::Event::update() {

	bool success = true;
	for(unsigned int i = 0; i < _functions.size(); ++i) {
		success = success and (*_functions[i])();
	}
	return success;

};

