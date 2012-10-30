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

		TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks,
		             TH1* hist_template,
		             T* data1,
		             T* data2 = 0);

		TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks,
                     TH1* hist_template,
                     std::vector<T*>* data1,
                     std::vector<T*>* data2 = 0);

		void fill(long cutmask);

	  private:

		void makePlot(std::map<std::string, std::vector<long> >& cutmasks, TH1* histTemplate);

		std::vector<std::pair<TH1*, long> > _histograms;

		unsigned int _mode;

		std::vector<T*> _vecData1;
		std::vector<T*> _vecData2;

		T* _data1;
		T* _data2;

	};

}

template<typename T>
antok::TemplatePlot<T>::TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks, TH1* histTemplate, T* data1, T* data2)
	: _data1(data1),
	  _data2(data2)
{

	assert(histTemplate != 0);
	assert(_data1 != 0);
	_mode = 0;
	makePlot(cutmasks, histTemplate);

};

template<typename T>
antok::TemplatePlot<T>::TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks,
                                     TH1* histTemplate,
                                     std::vector<T*>* vecData1,
                                     std::vector<T*>* vecData2)
{

	assert(histTemplate != 0);
	_vecData1 = *vecData1;
	if(vecData2 == 0) {
		_mode = 1;
	} else {
		assert(vecData1->size() == vecData2->size());
		_vecData2 = *vecData2;
		_mode = 2;
	}
	makePlot(cutmasks, histTemplate);

};


template<typename T>
void antok::TemplatePlot<T>::fill(long cutPattern) {

	for(unsigned int i = 0; i < _histograms.size(); ++i) {
		TH1* hist = _histograms[i].first;
		long histMask = _histograms[i].second;
		if((histMask&cutPattern) == histMask) {
			switch(_mode) {
				case 0:
					if(_data2 == 0) {
						hist->Fill(*_data1);
					} else {
						hist->Fill(*_data1, *_data2);
					}
					break;
				case 1:
					for(unsigned int j = 0; j < _vecData1.size(); ++j) {
						hist->Fill(*_vecData1[j]);
					}
					break;
				case 2:
					for(unsigned int j = 0; j < _vecData1.size(); ++j) {
						hist->Fill(*_vecData1[j], *_vecData2[j]);
					}
					break;
				default:
					throw 1;
			}
		}
	}

}

template<typename T>
void antok::TemplatePlot<T>::makePlot(std::map<std::string, std::vector<long> >& cutmasks, TH1* histTemplate)
{

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

}

#endif

