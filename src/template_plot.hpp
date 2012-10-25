#ifndef ANTOK_TEMPLATE_PLOT_HPP
#define ANTOK_TEMPLATE_PLOT_HPP

#include<assert.h>
#include<map>
#include<vector>

#include<TH1.h>
#include<TFile.h>

#include<cut.hpp>
#include<cutter.h>
#include<plot.hpp>

namespace antok {

	template<typename T>
	class TemplatePlot : public Plot {

	  public:

		TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks, TH1* hist_template, T* data1, T* data2 = 0);

		void fill(long cutmask);

	  private:
		TH1* _histTemplate;
		std::vector<std::pair<TH1*, long> > _histograms;

		T* _data1;
		T* _data2;

	};

}
#include<iostream>
template<typename T>
antok::TemplatePlot<T>::TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks, TH1* histTemplate, T* data1, T* data2)
	: _histTemplate(histTemplate),
	  _data1(data1),
	  _data2(data2)
{

	assert(_histTemplate != 0);
	assert(_data1 != 0);

	antok::ObjectManager* objectManager = antok::ObjectManager::instance();
	antok::Cutter& cutter = objectManager->getCutter();
	TFile* outFile = objectManager->getOutFile();

	for(std::map<std::string, std::vector<long> >::const_iterator cutmasks_it = cutmasks.begin(); cutmasks_it != cutmasks.end(); ++cutmasks_it) {
		const std::string& cutTrainName = cutmasks_it->first;
		const std::vector<long>& masks = cutmasks_it->second;
		const std::vector<antok::Cut*>& cuts = cutter.getCutsForCutTrain(cutTrainName);

		histTemplate->SetDirectory(0);
		outFile->cd(cutTrainName.c_str());
		TDirectory* dir = TDirectory::CurrentDirectory();
		dir->mkdir(histTemplate->GetName());
		dir->cd(histTemplate->GetName());

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

	outFile->cd();
	histTemplate->Delete();

};

template<typename T>
void antok::TemplatePlot<T>::fill(long cutPattern) {

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

#endif
