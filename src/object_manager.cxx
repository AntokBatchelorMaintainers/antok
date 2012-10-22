#include<object_manager.h>

#include<iostream>

#include<TFile.h>
#include<TObject.h>

#include<cutter.h>
#include<event.h>
#include<plotter.h>

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
	  _outFile(0),
	  _inTree(0)
{

}

bool antok::ObjectManager::magic() {

	bool success = _event->update() and _cutter->cut();
	long cutPattern = _cutter->getCutPattern();
	_plotter->fill(cutPattern);
	return success;

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

bool antok::ObjectManager::setOutFile(TFile* outFile) {

	if(outFile == 0) {
		return false;
	}
	_outFile = outFile;
	return true;

}

bool antok::ObjectManager::registerObjectToWrite(TObject* object) {

	bool found = false;
	for(unsigned int i = 0; i < _objectsToWrite.size(); ++i) {
		if(object == _objectsToWrite[i]) {
			found = true;
		}
	}
	if(found) {
		return false;
	}
	_objectsToWrite.push_back(object);
	return true;

}

bool antok::ObjectManager::finish() {

	bool success = true;
	if(_outFile->Write() == 0) {
		success = false;
	}
	_outFile->Close();
	_inFile->Close();
	return success;

}

