
#include<assert.h>
#include<sstream>

#include<TFile.h>

#include<cutter.h>
#include<plot.h>
#include<object_manager.h>

antok::Plot::Plot(std::map<std::string, std::vector<long> >& cutmasks, TH1* histTemplate, double* data1, double* data2)
	: _histTemplate(histTemplate),
	  _data1(data1),
	  _data2(data2)
{

	assert(_histTemplate != 0);
	assert(_data1 != 0);

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	antok::Cutter& cutter = objectManager->getCutter();

	for(std::map<std::string, std::vector<long> >::const_iterator cutmasks_it = cutmasks.begin(); cutmasks_it != cutmasks.end(); cutmasks_it++) {
		const std::string& cutTrainName = cutmasks_it->first;
		const std::vector<long>& masks = cutmasks_it->second;
		const std::vector<antok::Cut*>& cuts = cutter.getCutsForCutTrain(cutTrainName);

		for(unsigned int cutmask_i = 0; cutmask_i < masks.size(); ++cutmask_i) {

			long mask = masks[cutmask_i];
			std::stringstream strStr;
			strStr<<histTemplate->GetName()<<"_";

			for(unsigned int i = 0; i < cuts.size(); ++i) {
				if(cutter.cutOnInCutmask(mask, cuts[i])) {
					strStr<<"1";
				} else {
					strStr<<"0";
				}
			}

			TFile* outFile = objectManager->getOutFile();
			outFile->cd(cutTrainName.c_str());
			TH1* new_hist = dynamic_cast<TH1*>(histTemplate->Clone(strStr.str().c_str()));
			assert(new_hist != 0);
			strStr.str("");
			strStr<<histTemplate->GetTitle();
			strStr<<" "<<cutter.getAbbreviations(mask, cutTrainName);
			new_hist->SetTitle(strStr.str().c_str());
			assert(objectManager->registerObjectToWrite(new_hist));
			_histograms.push_back(std::pair<TH1*, long>(new_hist, mask));

		}

	}

};

void antok::Plot::fill(long cutPattern) {

	for(unsigned int i = 0; i < _histograms.size(); ++i) {
		TH1* hist = _histograms[i].first;
		long histMask = _histograms[i].second;
		if((histMask&cutPattern) == histMask) {
			if(_data2 == 0) {
				hist->Fill(*_data1);
			} else {
				hist->Fill(*_data1, *_data2);
			}
		}
	}

}

