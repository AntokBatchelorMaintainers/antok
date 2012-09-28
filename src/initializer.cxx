
#include<cutter.h>
#include<event.h>
#include<initializer.h>
#include<plotter.h>

antok::Initializer* antok::Initializer::_initializer = 0;

antok::Initializer* antok::Initializer::instance() {

	if(_initializer == 0) {
		_initializer = new antok::Initializer();
	};
	return _initializer;

};

antok::Initializer::Initializer()
	: _cutter(0),
	  _event(0),
	  _plotter(0)
{

};

bool antok::Initializer::readConfigFile(const std::string& filename) {

	return true;

};

antok::Cutter& antok::Initializer::get_cutter() {

	if(_cutter == 0) {
		_cutter = antok::Cutter::instance();
	}
	return *_cutter;

};

antok::Event& antok::Initializer::get_event() {

	if(_event == 0) {
		_event = antok::Event::instance();
	}
	return *_event;

};

antok::Plotter& antok::Initializer::get_plotter() {

	if(_plotter == 0) {
		_plotter = antok::Plotter::instance();
	}
	return *_plotter;

};

