#ifndef ANTOK_TEMPLATE_PLOT_HPP
#define ANTOK_TEMPLATE_PLOT_HPP

#include<assert.h>
#include<iostream>
#include<map>
#include<vector>

#include<TH1.h>
#include<TFile.h>

#include<cut.hpp>
#include<cutter.h>
#include<plot.hpp>

namespace antok {

	/**
	 * Template class for 1D, 2D, or 3D plots for different data types.
	 * Different modes are possible:
	 *  1D, 2D, 3D plots:
	 *     - mSiscalar:  All single scalar variables
	 *     - mMulscalar: All multiple scalar variables
	 *     - mSivector:  All single vector variables
	 *     - mMulvector: All multiple vector variables
	 *  2D plots:
	 *     - mSiscalarVsSivector:   Single scalar   vs single vector variables
	 *     - mSivectorVsSiscalar:   Single vector   vs single scalar variables
	 *     - mMulscalarVsMulvector: Multiple scalar vs multiple vector variables
	 *     - mMulvectorVsMulscalar: Multiple vector vs multiple scalar variables
	 *
	 */
	template<typename T1, typename T2 = T1, typename T3 = T2>
	class TemplatePlot : public Plot {
		enum Modes {
			mSiscalar, mMulscalar, mSivector, mMulvector,
			mSiscalarVsSivector, mSivectorVsSiscalar,
			mMulscalarVsMulvector, mMulvectorVsMulscalar
		};

	public:

		TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
		             TH1 *hist_template,
		             T1 *data1,
		             T2 *data2 = nullptr,
		             T3 *data3 = nullptr);

		TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
		             TH1 *hist_template,
		             std::vector<T1> *data1,
		             std::vector<T2> *data2 = nullptr,
		             std::vector<T3> *data3 = nullptr);

		TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
		             TH1 *hist_template,
		             std::vector<T1 *> *data1,
		             std::vector<T2 *> *data2 = nullptr,
		             std::vector<T3 *> *data3 = nullptr);

		TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
		             TH1 *hist_template,
		             std::vector<std::vector<T1> *> *data1,
		             std::vector<std::vector<T2> *> *data2 = nullptr,
		             std::vector<std::vector<T3> *> *data3 = nullptr);

		TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
		             TH1 *hist_template,
		             T1 *data1,
		             std::vector<T2> *data2);

		TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
		             TH1 *hist_template,
		             std::vector<T1> *data1,
		             T2 *data2);

		TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
		             TH1 *hist_template,
		             std::vector<T1 *> *data1,
		             std::vector<std::vector<T2> *> *data2);

		TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
		             TH1 *hist_template,
		             std::vector<std::vector<T1> *> *data1,
		             std::vector<T2 *> *data2);


		virtual void fill(long cutmask);

		~TemplatePlot() {};

	private:

		void makePlot(std::map<std::string, std::vector<long> > &cutmasks, TH1 *histTemplate);

		void fill(TH1 *hist, const T1 d1);
		void fill(TH1 *hist, const T1 d1, const T2 d2);
		void fill(TH1 *hist, const T1 d1, const T2 d2, const T3 d3);
		void fill(TH1 *hist, const T1 *d1) { fill(hist, *d1); }
		void fill(TH1 *hist, const T1 *d1, const T2 *d2) { fill(hist, *d1, *d2); }
		void fill(TH1 *hist, const T1 *d1, const T2 *d2, const T3 *d3) { fill(hist, *d1, *d2, *d3); }

		std::vector<std::pair<TH1 *, long> > _histograms;

		std::map<long, TH1 *> _cutmaskIndex;

		Modes _mode;

		std::vector<T1 *> *_multipleData1;
		std::vector<T2 *> *_multipleData2;
		std::vector<T3 *> *_multipleData3;

		std::vector<T1> *_dataVector1;
		std::vector<T2> *_dataVector2;
		std::vector<T3> *_dataVector3;

		std::vector<std::vector<T1> *> *_multipleDataVector1;
		std::vector<std::vector<T2> *> *_multipleDataVector2;
		std::vector<std::vector<T3> *> *_multipleDataVector3;

		T1 *_data1;
		T2 *_data2;
		T3 *_data3;

	};
}

template<typename T1, typename T2, typename T3>
antok::TemplatePlot<T1, T2, T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
                                              TH1 *histTemplate,
                                              T1 *data1,
                                              T2 *data2,
                                              T3 *data3)
		: Plot(),
		  _mode(mSiscalar),
		  _dataVector1(nullptr),
		  _dataVector2(nullptr),
		  _dataVector3(nullptr),
		  _multipleDataVector1(nullptr),
		  _multipleDataVector2(nullptr),
		  _multipleDataVector3(nullptr),
		  _data1(data1),
		  _data2(data2),
		  _data3(data3) {

	assert(histTemplate != nullptr);
	assert(_data1 != nullptr);
	makePlot(cutmasks, histTemplate);

};

template<typename T1, typename T2, typename T3>
antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
                                            TH1 *histTemplate,
                                            std::vector<T1*> *vecData1,
                                            std::vector<T2*> *vecData2,
                                            std::vector<T3*> *vecData3)
		: Plot(),
		  _mode(mMulscalar),
		  _multipleData1(vecData1),
		  _multipleData2(vecData2),
		  _multipleData3(vecData3),
		  _dataVector1(nullptr),
		  _dataVector2(nullptr),
		  _dataVector3(nullptr),
		  _multipleDataVector1(nullptr),
		  _multipleDataVector2(nullptr),
		  _multipleDataVector3(nullptr),
		  _data1(nullptr),
		  _data2(nullptr),
		  _data3(nullptr) {

	assert(histTemplate != nullptr);
	assert(_multipleData1 != nullptr);
	makePlot(cutmasks, histTemplate);

};

template<typename T1, typename T2, typename T3>
antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
                                            TH1 *histTemplate,
                                            std::vector<T1> *vecData1,
                                            std::vector<T2> *vecData2,
                                            std::vector<T3> *vecData3)
		: Plot(),
		  _mode(mSivector),
		  _multipleData1(nullptr),
		  _multipleData2(nullptr),
		  _multipleData3(nullptr),
		  _dataVector1(vecData1),
		  _dataVector2(vecData2),
		  _dataVector3(vecData3),
		  _multipleDataVector1(nullptr),
		  _multipleDataVector2(nullptr),
		  _multipleDataVector3(nullptr),
		  _data1(nullptr),
		  _data2(nullptr),
		  _data3(nullptr) {

	assert(histTemplate != nullptr);
	assert(_dataVector1 != nullptr);
	makePlot(cutmasks, histTemplate);

};

template<typename T1, typename T2, typename T3>
antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
                                            TH1 *histTemplate,
                                            std::vector<std::vector<T1>*> *vecData1,
                                            std::vector<std::vector<T2>*> *vecData2,
                                            std::vector<std::vector<T3>*> *vecData3)
		: Plot(),
		  _mode(mMulvector),
		  _multipleData1(nullptr),
		  _multipleData2(nullptr),
		  _multipleData3(nullptr),
		  _dataVector1(nullptr),
		  _dataVector2(nullptr),
		  _dataVector3(nullptr),
		  _multipleDataVector1(vecData1),
		  _multipleDataVector2(vecData2),
		  _multipleDataVector3(vecData3),
		  _data1(nullptr),
		  _data2(nullptr),
		  _data3(nullptr) {

	assert(histTemplate != 0);
	assert(_multipleDataVector1 != nullptr);
	makePlot(cutmasks, histTemplate);

};

template<typename T1, typename T2, typename T3>
antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
                                            TH1 *histTemplate,
                                            T1* data1,
                                            std::vector<T2> *vecData2)
		: Plot(),
		  _mode(mSiscalarVsSivector),
		  _multipleData1(nullptr),
		  _multipleData2(nullptr),
		  _multipleData3(nullptr),
		  _dataVector1(nullptr),
		  _dataVector2(vecData2),
		  _dataVector3(nullptr),
		  _multipleDataVector1(nullptr),
		  _multipleDataVector2(nullptr),
		  _multipleDataVector3(nullptr),
		  _data1(data1),
		  _data2(nullptr),
		  _data3(nullptr) {

	assert(histTemplate != nullptr);
	assert(_data1 != nullptr);
	assert(_dataVector2 != nullptr);
	makePlot(cutmasks, histTemplate);

};

template<typename T1, typename T2, typename T3>
antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
                                            TH1 *histTemplate,
                                            std::vector<T1> *vecData1,
                                            T2* data2)
		: Plot(),
		  _mode(mSivectorVsSiscalar),
		  _multipleData1(nullptr),
		  _multipleData2(nullptr),
		  _multipleData3(nullptr),
		  _dataVector1(vecData1),
		  _dataVector2(nullptr),
		  _dataVector3(nullptr),
		  _multipleDataVector1(nullptr),
		  _multipleDataVector2(nullptr),
		  _multipleDataVector3(nullptr),
		  _data1(nullptr),
		  _data2(data2),
		  _data3(nullptr) {

	assert(histTemplate != nullptr);
	assert(_dataVector1 != nullptr);
	assert(_data2 != nullptr);
	makePlot(cutmasks, histTemplate);

};

template<typename T1, typename T2, typename T3>
antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
                                            TH1 *histTemplate,
                                            std::vector<T1*> *vecData1,
                                            std::vector<std::vector<T2>*> *vecData2)
		: Plot(),
		  _mode(mMulscalarVsMulvector),
		  _multipleData1(vecData1),
		  _multipleData2(nullptr),
		  _multipleData3(nullptr),
		  _dataVector1(nullptr),
		  _dataVector2(nullptr),
		  _dataVector3(nullptr),
		  _multipleDataVector1(nullptr),
		  _multipleDataVector2(vecData2),
		  _multipleDataVector3(nullptr),
		  _data1(nullptr),
		  _data2(nullptr),
		  _data3(nullptr) {

	assert(histTemplate != nullptr);
	assert(_dataVector1 != nullptr);
	assert(_multipleDataVector2 != nullptr);
	assert(_dataVector1->size() != _multipleDataVector2->size());
	makePlot(cutmasks, histTemplate);

};

template<typename T1, typename T2, typename T3>
antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
                                            TH1 *histTemplate,
                                            std::vector<std::vector<T1>*> *vecData1,
                                            std::vector<T2*> *vecData2)
		: Plot(),
		  _mode(mMulvectorVsMulscalar),
		  _multipleData1(nullptr),
		  _multipleData2(vecData2),
		  _multipleData3(nullptr),
		  _dataVector1(nullptr),
		  _dataVector2(nullptr),
		  _dataVector3(nullptr),
		  _multipleDataVector1(vecData1),
		  _multipleDataVector2(nullptr),
		  _multipleDataVector3(nullptr),
		  _data1(nullptr),
		  _data2(nullptr),
		  _data3(nullptr) {

	assert(histTemplate != nullptr);
	assert(_multipleDataVector1 != nullptr);
	assert(_dataVector2 != nullptr);
	assert(_multipleDataVector1->size() != _dataVector2->size());
	makePlot(cutmasks, histTemplate);

};

template<typename T1, typename T2, typename T3>
void antok::TemplatePlot<T1,T2,T3>::fill(TH1* hist, const T1 d1) {
	hist->Fill(static_cast<Double_t>(d1));
}
template<typename T1, typename T2, typename T3>
void antok::TemplatePlot<T1,T2,T3>::fill(TH1* hist, const T1 d1, const T2 d2) {
	hist->Fill(static_cast<Double_t>(d1), static_cast<Double_t>(d2));
}
template<typename T1, typename T2, typename T3>
void antok::TemplatePlot<T1,T2,T3>::fill(TH1* hist, const T1 d1, const T2 d2, const T3 d3) {
	static_cast<TH3D*>(hist)->Fill(static_cast<Double_t>(d1), static_cast<Double_t>(d2), static_cast<Double_t>(d3));
}

template<typename T1, typename T2, typename T3>
void antok::TemplatePlot<T1,T2,T3>::fill(long cutPattern) {

	for (unsigned int i = 0; i < _histograms.size(); ++i) {
		TH1 *hist = _histograms[i].first;
		long histMask = _histograms[i].second;
		if ((histMask & cutPattern) == histMask) {
			switch (_mode) {
				case mSiscalar:
					if (_data2 == nullptr) {
						fill(hist,_data1);
					} else if (_data3 == nullptr) {
						fill(hist,_data1, _data2);
					} else {
						fill(hist,_data1, _data2, _data3);
					}
					break;
				case mMulscalar:
					if( _multipleData2 == nullptr ) {
						for (unsigned int j = 0; j < _multipleData1->size(); ++j) {
							fill(hist,(*_multipleData1)[j]);
						}
					} else if( _multipleData3 == nullptr ) {
						assert(_multipleData1->size() == _multipleData2->size());
						for (unsigned int j = 0; j < _multipleData1->size(); ++j) {
							fill(hist,(*_multipleData1)[j], (*_multipleData2)[j]);
						}
					} else {
						assert(_multipleData1->size() == _multipleData2->size());
						assert(_multipleData1->size() == _multipleData3->size());
						for (unsigned int j = 0; j < _multipleData1->size(); ++j) {
							fill(hist,(*_multipleData1)[j], (*_multipleData2)[j], (*_multipleData3)[j]);
						}
					}
					break;
				case mSivector:
					if( _dataVector2 == nullptr ) {
						for (unsigned int j = 0; j < _dataVector1->size(); ++j) {
							fill(hist,(*_dataVector1)[j]);
						}
					} else if( _dataVector3 == nullptr ) {
						assert(_dataVector1->size() == _dataVector2->size());
						for (unsigned int j = 0; j < _dataVector1->size(); ++j) {
							fill(hist,(*_dataVector1)[j], (*_dataVector2)[j]);
						}
					} else {
						assert(_dataVector1->size() == _dataVector2->size());
						assert(_dataVector1->size() == _dataVector3->size());
						for (unsigned int j = 0; j < _dataVector1->size(); ++j) {
							fill(hist,(*_dataVector1)[j], (*_dataVector2)[j], (*_dataVector3)[j]);
						}
					}
					break;
				case mMulvector:
					if( _multipleDataVector2 == nullptr ) {
						for ( unsigned int j = 0; j < _multipleDataVector1->size(); ++j ) {
							for(unsigned int k = 0; k < (*_multipleDataVector1)[j]->size(); ++k ) {
								fill(hist,(*(*_multipleDataVector1)[j])[k]);
							}
						}
					} else if( _multipleDataVector3 == nullptr ) {
						assert(_multipleDataVector1->size() == _multipleDataVector2->size());
						for ( unsigned int j = 0; j < _multipleDataVector1->size(); ++j ) {
							assert( (*_multipleDataVector1)[j]->size() == (*_multipleDataVector2)[j]->size() );
							for(unsigned int k = 0; k < (*_multipleDataVector1)[j]->size(); ++k ) {
								fill(hist,(*(*_multipleDataVector1)[j])[k],(*(*_multipleDataVector2)[j])[k]);
							}
						}
					} else {
						assert(_multipleDataVector1->size() == _multipleDataVector2->size());
						assert(_multipleDataVector1->size() == _multipleDataVector3->size());
						for ( unsigned int j = 0; j < _multipleDataVector1->size(); ++j ) {
							assert( (*_multipleDataVector1)[j]->size() == (*_multipleDataVector2)[j]->size() );
							assert( (*_multipleDataVector1)[j]->size() == (*_multipleDataVector3)[j]->size() );
							for(unsigned int k = 0; k < (*_multipleDataVector1)[j]->size(); ++k ) {
								fill(hist,(*(*_multipleDataVector1)[j])[k],(*(*_multipleDataVector2)[j])[k],(*(*_multipleDataVector3)[j])[k]);
							}
						}
					}
					break;
				case mSiscalarVsSivector:
					for( unsigned int j = 0; j < _dataVector2->size(); ++j ) {
						fill(hist,*_data1, (*_dataVector2)[j]);
					}
					break;
				case mSivectorVsSiscalar:
					for( unsigned int j = 0; j < _dataVector1->size(); ++j ) {
						fill(hist,(*_dataVector1)[j], *_data2);
					}
					break;
				case mMulscalarVsMulvector:
					for( unsigned int j = 0; j < _dataVector1->size(); ++j ) {
						for( unsigned int k = 0; k < (*_multipleDataVector2)[j]->size(); ++k ) {
							fill(hist,*(*_multipleData1)[j], (*(*_multipleDataVector2)[j])[k]);
						}
					}
					break;
				case mMulvectorVsMulscalar:
					for( unsigned int j = 0; j < _dataVector2->size(); ++j ) {
						for( unsigned int k = 0; k < (*_multipleDataVector1)[j]->size(); ++k ) {
							fill(hist, (*(*_multipleDataVector1)[j])[k],  *(*_multipleData2)[j]);
						}
					}
					break;

				default:
					throw 1;
			}
		}
	}

}

template<typename T1, typename T2, typename T3>
void
antok::TemplatePlot<T1, T2, T3>::makePlot(std::map<std::string, std::vector<long> > &cutmasks, TH1 *histTemplate) {

	antok::ObjectManager *objectManager = antok::ObjectManager::instance();
	antok::Cutter &cutter = objectManager->getCutter();
	TFile *outFile = objectManager->getOutFile();

	for (std::map<std::string, std::vector<long> >::const_iterator cutmasks_it = cutmasks.begin(); cutmasks_it != cutmasks.end(); ++cutmasks_it) {
		const std::string &cutTrainName = cutmasks_it->first;
		const std::vector<long> &masks = cutmasks_it->second;
		const std::vector<antok::Cut *> &cuts = cutter.getCutsForCutTrain(cutTrainName);

		histTemplate->SetDirectory(0);
		std::stringstream strStr;
		strStr << "tmptmptmp/" << cutTrainName;
		outFile->cd(strStr.str().c_str());
		TDirectory *dir = TDirectory::CurrentDirectory();
		dir->mkdir(histTemplate->GetName());
		dir->cd(histTemplate->GetName());

		for (unsigned int cutmask_i = 0; cutmask_i < masks.size(); ++cutmask_i) {

			long mask = masks[cutmask_i];
			strStr.str("");
			strStr << histTemplate->GetName() << "_";

			for (unsigned int i = 0; i < cuts.size(); ++i) {
				if (cutter.cutOnInCutmask(mask, cuts[i])) {
					strStr << "1";
				} else {
					strStr << "0";
				}
			}
			std::string histName = strStr.str();

			strStr.str("");
			strStr << histTemplate->GetTitle();
			strStr << " " << cutter.getAbbreviations(mask, cutTrainName);
			std::string histTitle = strStr.str();

			strStr.str("");
			strStr << cutTrainName << "/" << histTemplate->GetName();
			std::string path = strStr.str();

			TH1 *hist = 0;
			if (_cutmaskIndex.find(mask) == _cutmaskIndex.end()) {
				hist = dynamic_cast<TH1 *>(histTemplate->Clone(histName.c_str()));
				assert(hist != 0);
				hist->SetTitle(histTitle.c_str());
				_histograms.push_back(std::pair<TH1 *, long>(hist, mask));
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
