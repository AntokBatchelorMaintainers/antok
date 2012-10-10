#ifndef ANTOK_OBJECTMANAGER_H
#define ANTOK_OBJECTMANAGER_H

class TFile;
class TTree;

namespace antok {

	class Cutter;
	class Data;
	class Event;
	class Plotter;

	class ObjectManager {

		friend class Initializer;

	  public:

		static ObjectManager* instance();

		antok::Cutter& getCutter();
		antok::Data& getData();
		antok::Event& getEvent();
		antok::Plotter& getPlotter();

		TFile* getInFile() { return _inFile; };
		TFile* getOutFile() { return _outFile; };
		TTree* getInTree() { return _inTree; };

		bool setInFile(TFile* inFile);
		bool setOutFile(TFile* outFile);

	  private:

		ObjectManager();

		static ObjectManager* _objectManager;

		antok::Cutter* _cutter;
		antok::Data* _data;
		antok::Event* _event;
		antok::Plotter* _plotter;

		TFile* _inFile;
		TFile* _outFile;
		TTree* _inTree;

	};

}

#endif

