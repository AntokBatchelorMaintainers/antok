#include<object_manager.h>

#include<iostream>

antok::ObjectManager* antok::ObjectManager::_objectManager = 0;

antok::ObjectManager* antok::ObjectManager::instance() {

	if(_objectManager == 0) {
		_objectManager = new antok::ObjectManager();
	}
	return _objectManager;

}

antok::ObjectManager::ObjectManager()
	: _cutter(0),
	  _data(0),
	  _event(0),
	  _plotter(0),
	  _inFile(0),
	  _inTree(0)
{

}

antok::Cutter& antok::ObjectManager::getCutter() {

	if(_cutter == 0) {
		std::cerr<<"Trying to get uninitialized Cutter."<<std::endl;
		throw 1;
	}
	return *_cutter;

}

antok::Data& antok::ObjectManager::getData() {

	if(_data == 0) {
		std::cerr<<"Trying to get uninitialized Data."<<std::endl;
		throw 1;
	}
	return *_data;

}

antok::Event& antok::ObjectManager::getEvent() {

	if(_event == 0) {
		std::cerr<<"Trying to get uninitialized Event."<<std::endl;
		throw 1;
	}
	return *_event;

}

antok::Plotter& antok::ObjectManager::getPlotter() {

	if(_plotter == 0) {
		std::cerr<<"Trying to get uninitialized Plotter."<<std::endl;
		throw 1;
	}
	return *_plotter;

}

bool antok::ObjectManager::setInFile(TFile* inFile) {

	if(inFile == 0) {
		return false;
	}
	_inFile = inFile;
	return true;

}

