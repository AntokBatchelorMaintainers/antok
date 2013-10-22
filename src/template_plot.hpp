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

		TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks,
                     TH1* hist_template,
                     std::vector<T>* data1,
                     std::vector<T>* data2 = 0);

		TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks,
                     TH1* hist_template,
                     std::vector<std::vector<T>*>* data1,
                     std::vector<std::vector<T>*>* data2 = 0);

		void fill(long cutmask);

		~TemplatePlot() { };

	  private:

		void makePlot(std::map<std::string, std::vector<long> >& cutmasks, TH1* histTemplate);

		std::vector<std::pair<TH1*, long> > _histograms;

		std::map<long, TH1*> _cutmaskIndex;

		unsigned int _mode;

		std::vector<T*>* _vecData1;
		std::vector<T*>* _vecData2;

		std::vector<T>* _vecDataVector1;
		std::vector<T>* _vecDataVector2;

		std::vector<std::vector<T>*>* _multipleVecDataVectors1;
		std::vector<std::vector<T>*>* _multipleVecDataVectors2;

		T* _data1;
		T* _data2;

	};

}

template<typename T>
antok::TemplatePlot<T>::TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks,
                                     TH1* histTemplate,
                                     T* data1,
                                     T* data2)
	: Plot(),
	  _vecDataVector1(0),
	  _vecDataVector2(0),
	  _multipleVecDataVectors1(0),
	  _multipleVecDataVectors2(0),
	  _data1(data1),
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
	: Plot(),
	  _vecData1(vecData1),
	  _vecData2(vecData2),
	  _vecDataVector1(0),
	  _vecDataVector2(0),
	  _multipleVecDataVectors1(0),
	  _multipleVecDataVectors2(0),
	  _data1(0),
	  _data2(0)
{

	assert(histTemplate != 0);
	if(_vecData2 == 0) {
		_mode = 1;
	} else {
		assert(_vecData1->size() == _vecData2->size());
		_mode = 2;
	}
	makePlot(cutmasks, histTemplate);

};

template<typename T>
antok::TemplatePlot<T>::TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks,
                                     TH1* histTemplate,
                                     std::vector<T>* vecData1,
                                     std::vector<T>* vecData2)
	: Plot(),
	  _vecData1(0),
	  _vecData2(0),
	  _vecDataVector1(vecData1),
	  _vecDataVector2(vecData2),
	  _multipleVecDataVectors1(0),
	  _multipleVecDataVectors2(0),
	  _data1(0),
	  _data2(0)
{

	assert(histTemplate != 0);
	if(_vecDataVector2 == 0) {
		_mode = 3;
	} else {
		assert(_vecDataVector1->size() == _vecDataVector2->size());
		_mode = 4;
	}
	makePlot(cutmasks, histTemplate);

};

template<typename T>
antok::TemplatePlot<T>::TemplatePlot(std::map<std::string, std::vector<long> >& cutmasks,
                                     TH1* histTemplate,
                                     std::vector<std::vector<T>*>* vecData1,
                                     std::vector<std::vector<T>*>* vecData2)
	: Plot(),
	  _vecData1(0),
	  _vecData2(0),
	  _vecDataVector1(0),
	  _vecDataVector2(0),
	  _multipleVecDataVectors1(vecData1),
	  _multipleVecDataVectors2(vecData2),
	  _data1(0),
	  _data2(0)
{

	assert(histTemplate != 0);
	if(_multipleVecDataVectors2 == 0) {
		_mode = 5;
	} else {
		assert(_multipleVecDataVectors1->size() == _multipleVecDataVectors2->size());
		_mode = 6;
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
				case 0: // Pointers to single variables
					if(_data2 == 0) {
						hist->Fill(*_data1);
					} else {
						hist->Fill(*_data1, *_data2);
					}
					break;
				case 1: // Multiple values, 1 variable
					for(unsigned int j = 0; j < _vecData1->size(); ++j) {
						hist->Fill(*(*_vecData1)[j]);
					}
					break;
				case 2: // Multiple values, 2 variables
					for(unsigned int j = 0; j < _vecData1->size(); ++j) {
						hist->Fill(*(*_vecData1)[j], *(*_vecData2)[j]);
					}
					break;
				case 3: // Multiple values of varying size, 1 variable
					for(unsigned int j = 0; j < _vecDataVector1->size(); ++j) {
						hist->Fill((*_vecDataVector1)[j]);
					}
					break;
				case 4: // Multiple values of varying size, 2 variables
					if(_vecDataVector1->size() != _vecDataVector2->size()) {
						std::cerr<<"When filling plot with std::vectors as variables: "
						         <<"vectors have different size ("
						         <<_vecDataVector1->size()<<"!="<<_vecDataVector2->size()
						         <<"). Aborting..."<<std::endl;
						throw;
					}
					for(unsigned int j = 0; j < _vecDataVector1->size(); ++j) {
						hist->Fill((*_vecDataVector1)[j], (*_vecDataVector2)[j]);
					}
					break;
				case 5: // Multiple vectors with multiple values of varying size, 1 variable
					for(unsigned int j = 0; j < _multipleVecDataVectors1->size(); ++j) {
						for(unsigned int k = 0; k < (*_multipleVecDataVectors1)[j]->size(); ++k) {
							hist->Fill((*(*_multipleVecDataVectors1)[j])[k]);
						}
					}
					break;
				case 6: // Multiple vectors with multiple values of varying size, 2 variables
					for(unsigned int j = 0; j < _multipleVecDataVectors1->size(); ++j) {
						if((*_multipleVecDataVectors1)[j]->size() != (*_multipleVecDataVectors2)[j]->size()) {
							std::cerr<<"When filling plot with std::vectors as variables and "
							         <<"\"Indices\" used: sub-vectors for index "<<j<<" of different "
							         <<"size ("<<(*_multipleVecDataVectors1)[j]->size()<<"!="
							         <<(*_multipleVecDataVectors2)[j]->size()<<"). Aborting..."<<std::endl;
							throw;
						}
						for(unsigned int k = 0; k < (*_multipleVecDataVectors1)[j]->size(); ++k) {
							hist->Fill((*(*_multipleVecDataVectors1)[j])[k], (*(*_multipleVecDataVectors2)[j])[k]);
						}
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
		std::stringstream strStr;
		strStr<<"tmptmptmp/"<<cutTrainName;
		outFile->cd(strStr.str().c_str());
		TDirectory* dir = TDirectory::CurrentDirectory();
		dir->mkdir(histTemplate->GetName());
		dir->cd(histTemplate->GetName());

		for(unsigned int cutmask_i = 0; cutmask_i < masks.size(); ++cutmask_i) {

			long mask = masks[cutmask_i];
			strStr.str("");
			strStr<<histTemplate->GetName()<<"_";

			for(unsigned int i = 0; i < cuts.size(); ++i) {
				if(cutter.cutOnInCutmask(mask, cuts[i])) {
					strStr<<"1";
				} else {
					strStr<<"0";
				}
			}
			std::string histName = strStr.str();

			strStr.str("");
			strStr<<histTemplate->GetTitle();
			strStr<<" "<<cutter.getAbbreviations(mask, cutTrainName);
			std::string histTitle = strStr.str();

			strStr.str("");
			strStr<<cutTrainName<<"/"<<histTemplate->GetName();
			std::string path = strStr.str();

			TH1* hist = 0;
			if(_cutmaskIndex.find(mask) == _cutmaskIndex.end()) {
				hist = dynamic_cast<TH1*>(histTemplate->Clone(histName.c_str()));
				assert(hist != 0);
				hist->SetTitle(histTitle.c_str());
				_histograms.push_back(std::pair<TH1*, long>(hist, mask));
				_cutmaskIndex[mask] = hist;
			} else {
				hist = _cutmaskIndex.find(mask)->second;
			}
			assert(objectManager->registerHistogramToCopy(hist,
			                                              path,
			                                              histName,
			                                              histTitle));
		}

	}

	outFile->cd();
	histTemplate->Delete();

}

#endif

