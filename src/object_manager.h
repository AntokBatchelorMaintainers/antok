#ifndef ANTOK_OBJECTMANAGER_H
#define ANTOK_OBJECTMANAGER_H

#include<string>
#include<vector>

class TFile;
class TH1;
class TObject;
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

		bool magic();

		antok::Cutter& getCutter();
		antok::Data& getData();
		antok::Event& getEvent();
		antok::Plotter& getPlotter();

		TFile* getInFile() { return _inFile; };
		TFile* getOutFile() { return _outFile; };
		TTree* getInTree() { return _inTree; };

		bool setInFile(TFile* inFile);
		bool setOutFile(TFile* outFile);

		bool registerObjectToWrite(TObject* object);
		bool registerHistogramToCopy(TH1* histogram,
		                             std::string path,
		                             std::string newName,
		                             std::string newTitle);

		bool finish();

	  private:

		struct histogramCopyInformation {
			TH1* histogram;
			std::string path;
			std::string newName;
			std::string newTitle;
			histogramCopyInformation(TH1* hist, std::string inPath, std::string nName, std::string nTitle)
				: histogram(hist),
				  path(inPath),
				  newName(nName),
				  newTitle(nTitle) { };
		};

		ObjectManager();

		static ObjectManager* _objectManager;

		antok::Cutter* _cutter;
		antok::Data* _data;
		antok::Event* _event;
		antok::Plotter* _plotter;

		TFile* _inFile;
		TFile* _outFile;
		TTree* _inTree;
		std::vector<TObject*> _objectsToWrite;
		std::vector<histogramCopyInformation> _histogramsToCopy;

	};

}

#endif

