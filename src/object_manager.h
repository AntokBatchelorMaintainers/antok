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

	  public:

		static ObjectManager* instance();

		antok::Cutter& getCutter();
		antok::Data& getData();
		antok::Event& getEvent();
		antok::Plotter& getPlotter();

		TFile* getInFile() { return _inFile; };
		TTree* getInTree() { return _inTree; };

		bool setInFile(TFile* inFile);

	  private:

		static ObjectManager* _objectManager;

		ObjectManager();

		antok::Cutter* _cutter;
		antok::Data* _data;
		antok::Event* _event;
		antok::Plotter* _plotter;

		TFile* _inFile;
		TTree* _inTree;

		friend class Initializer;

	};

}

#endif

