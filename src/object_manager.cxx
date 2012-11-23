#include<object_manager.h>

#include<assert.h>
#include<iostream>
#include<sstream>

#include<TFile.h>
#include<TH1.h>
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

	bool success = _event->update() and _cutter->cut() and _cutter->fillOutTrees();
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

bool antok::ObjectManager::registerHistogramToCopy(TH1* histogram,
                                                   std::string cutTrainName,
                                                   std::string dirName,
                                                   std::string newName,
                                                   std::string newTitle,
                                                   unsigned int orderParameter)
{
	histogramCopyInformation histCopyInfo(histogram, newName, newTitle, orderParameter);
	std::stringstream strStr;
	strStr<<cutTrainName<<"/"<<dirName;
	_histogramsToCopy[strStr.str()].push_back(histCopyInfo);
	return true;
}

bool antok::ObjectManager::finish() {

	bool success = true;
	for(std::map<std::string, std::map<unsigned int, TH1*> >::const_iterator it = _histogramOrder.begin();
	    it != _histogramOrder.end();
	    ++it)
	{
		std::string dirName = it->first;
		TDirectory* dir = dynamic_cast<TDirectory*>(_outFile->Get(dirName.c_str()));
		assert(dir != 0);
		dir->cd();
		if(_histogramsToCopy.find(dirName) != _histogramsToCopy.end()) {
			std::vector<histogramCopyInformation>::const_iterator histogramsToCopy_it= _histogramsToCopy[dirName].begin();
			std::map<unsigned int, TH1*>::const_iterator orderOfExistingHistograms_it = it->second.begin();
			for(unsigned int i = 0;
			    orderOfExistingHistograms_it != it->second.end() || histogramsToCopy_it != _histogramsToCopy[dirName].end();
			    ++i)
			{
				if(orderOfExistingHistograms_it == it->second.end() || i == histogramsToCopy_it->orderParameter) {
					histogramCopyInformation info = (*histogramsToCopy_it);
					assert(i == info.orderParameter);
					TH1* copiedHist = dynamic_cast<TH1*>(info.histogram->Clone(info.newName.c_str()));
					assert(copiedHist != 0);
					copiedHist->SetTitle(info.newTitle.c_str());
					registerObjectToWrite(copiedHist);
					++histogramsToCopy_it;
				} else if(histogramsToCopy_it == _histogramsToCopy[dirName].end() || i == orderOfExistingHistograms_it->first) {
					TH1* hist = orderOfExistingHistograms_it->second;
					dir->Remove(hist);
					dir->Append(hist);
					++orderOfExistingHistograms_it;
				} else {
					assert(false);
				}
			}
		}
		_outFile->cd();
	}
	if(_outFile->Write() == 0) {
		success = false;
	}
	_outFile->Close();
	_inFile->Close();
	return success;

}

