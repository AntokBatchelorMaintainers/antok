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

antok::ObjectManager* antok::ObjectManager::_objectManager = nullptr;

antok::ObjectManager* antok::ObjectManager::instance() {

	if(_objectManager == nullptr) {
		_objectManager = new antok::ObjectManager();
	}
	return _objectManager;

}

antok::ObjectManager::ObjectManager()
	: _cutter(nullptr),
	  _data(nullptr),
	  _event(nullptr),
	  _plotter(nullptr),
	  _inFile(nullptr),
	  _outFile(nullptr),
	  _inTree(nullptr),
	  _histNameAppendix("")
{

}

bool antok::ObjectManager::magic() {

	bool success = _event->update() and _cutter->cut() and _cutter->fillOutTrees();
	const antok::bitmask& cutPattern = _cutter->getCutPattern();
	_plotter->fill(cutPattern);
	return success;

}

antok::Cutter& antok::ObjectManager::getCutter() {

	if(_cutter == nullptr) {
		std::cerr<<"Trying to get uninitialized Cutter."<<std::endl;
		throw 1;
	}
	return *_cutter;

}

antok::Data& antok::ObjectManager::getData() {

	if(_data == nullptr) {
		std::cerr<<"Trying to get uninitialized Data."<<std::endl;
		throw 1;
	}
	return *_data;

}

antok::Event& antok::ObjectManager::getEvent() {

	if(_event == nullptr) {
		std::cerr<<"Trying to get uninitialized Event."<<std::endl;
		throw 1;
	}
	return *_event;

}

antok::Plotter& antok::ObjectManager::getPlotter() {

	if(_plotter == nullptr) {
		std::cerr<<"Trying to get uninitialized Plotter."<<std::endl;
		throw 1;
	}
	return *_plotter;

}

bool antok::ObjectManager::setInFile(TFile* inFile) {

	if(inFile == nullptr) {
		return false;
	}
	_inFile = inFile;
	return true;

}
bool antok::ObjectManager::changeInFile(TFile* inFile){
	if( _inFile != nullptr ) _inFile->Close();
	return setInFile(inFile);
}

bool antok::ObjectManager::setOutFile(TFile* outFile) {

	if(outFile == nullptr) {
		return false;
	}
	_outFile = outFile;
	return true;

}

bool antok::ObjectManager::registerObjectToWrite(TDirectory* path, TObject* object) {

	if(_objectsToWrite.find(object) != _objectsToWrite.end()) {
		return false;
	}
	_objectsToWrite[object] = path;
	return true;

}

bool antok::ObjectManager::registerHistogramToCopy(TH1* histogram,
                                                   std::string path,
                                                   std::string newName,
                                                   std::string newTitle)
{
	histogramCopyInformation histCopyInfo(histogram, newName, newTitle);
	_histogramsToCopy[path].push_back(histCopyInfo);
	return true;
}

bool antok::ObjectManager::registerHistogramNameAppendix(const std::string& appendix) {
	_histNameAppendix = appendix;
	return true;
}

bool antok::ObjectManager::finish() {

	bool success = true;

	for(std::map<TObject*, TDirectory*>::const_iterator it = _objectsToWrite.begin(); it != _objectsToWrite.end(); ++it) {
		it->second->cd();
		if(it->first->Write() <= 0) {
			success = false;
		}
	}

	for(std::map<std::string, std::vector<histogramCopyInformation> >::const_iterator histsToCopy_it = _histogramsToCopy.begin();
	    histsToCopy_it != _histogramsToCopy.end();
	    ++histsToCopy_it)
	{
		std::string path = histsToCopy_it->first;
		std::vector<histogramCopyInformation> histsToCopy = histsToCopy_it->second;
		if(histsToCopy.size() == 0) {
			continue;
		}
		TDirectory* dir = _outFile->GetDirectory(path.c_str(), false);
		if(not dir) {
			std::string cutTrainDirName = path.substr(0, path.find_last_of('/'));
			std::string plotDirName = path.substr(path.find_last_of('/')+1);
			TDirectory* cutTrainDir = _outFile->GetDirectory(cutTrainDirName.c_str());
			assert(cutTrainDir != nullptr);
			dir = cutTrainDir->mkdir(plotDirName.c_str());
			assert(dir != nullptr);
		}
		dir->cd();
		for(unsigned int i = 0; i < histsToCopy.size(); ++i) {
			histogramCopyInformation info = histsToCopy[i];
			std::string histName = info.newName;
			if(_histNameAppendix != "") {
					std::stringstream strStr;
					strStr<<histName<<_histNameAppendix;
					histName = strStr.str();
			}
			TH1* copiedHist = dynamic_cast<TH1*>(info.histogram->Clone(histName.c_str()));
			assert(copiedHist != nullptr);
			copiedHist->SetTitle(info.newTitle.c_str());
			if(copiedHist->Write() <= 0) {
				success = false;
			}
		}
		dir->Close();
	}
	_outFile->cd();
	_outFile->Delete("tmptmptmp;*");

	_outFile->Close();
	_inFile->Close();
	return success;

}

