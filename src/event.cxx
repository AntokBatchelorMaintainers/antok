#include "event.h"
#include "basic_calcs.h"
#include "constants.h"
#include "functions.hpp"


antok::Event* antok::Event::_event = nullptr;


antok::Event*
antok::Event::instance()
{
	if (_event == nullptr) {
		_event = new antok::Event();
	}
	return _event;
}


bool
antok::Event::update()
{
	bool success = true;
	for (auto& function : _functions) {
		if (not (*function)()) {
			std::cerr << "Error evaluating function '" << function->name() << "'." << std::endl;
			success = false;
		}
	}
	return success;
};

