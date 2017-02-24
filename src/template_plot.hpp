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

	template<typename T1, typename T2 = T1, typename T3 = T2>
	class TemplatePlot : public Plot {

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
		             std::vector<T1*> *data1,
		             std::vector<std::vector<T2>*> *data2);

		TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
		             TH1 *hist_template,
		             std::vector<std::vector<T1>*> *data1,
		             std::vector<T2*> *data2);


		virtual void fill(long cutmask);

		~TemplatePlot() {};

	private:

		void makePlot(std::map<std::string, std::vector<long> > &cutmasks, TH1 *histTemplate);

		std::vector<std::pair<TH1 *, long> > _histograms;

		std::map<long, TH1 *> _cutmaskIndex;

		unsigned int _mode;

		std::vector<T1 *> *_vecData1;
		std::vector<T2 *> *_vecData2;
		std::vector<T3 *> *_vecData3;

		std::vector<T1> *_vecDataVector1;
		std::vector<T2> *_vecDataVector2;
		std::vector<T3> *_vecDataVector3;

		std::vector<std::vector<T1> *> *_multipleVecDataVectors1;
		std::vector<std::vector<T2> *> *_multipleVecDataVectors2;
		std::vector<std::vector<T3> *> *_multipleVecDataVectors3;

		T1 *_data1;
		T2 *_data2;
		T3 *_data3;

	};

	template<typename T1, typename T2, typename T3>
	antok::TemplatePlot<T1, T2, T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
	                                              TH1 *histTemplate,
	                                              T1 *data1,
	                                              T2 *data2,
	                                              T3 *data3)
			: Plot(),
			  _mode(0),
			  _vecDataVector1(nullptr),
			  _vecDataVector2(nullptr),
			  _vecDataVector3(nullptr),
			  _multipleVecDataVectors1(nullptr),
			  _multipleVecDataVectors2(nullptr),
			  _multipleVecDataVectors3(nullptr),
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
			  _mode(1),
			  _vecData1(vecData1),
			  _vecData2(vecData2),
			  _vecData3(vecData3),
			  _vecDataVector1(nullptr),
			  _vecDataVector2(nullptr),
			  _vecDataVector3(nullptr),
			  _multipleVecDataVectors1(nullptr),
			  _multipleVecDataVectors2(nullptr),
			  _multipleVecDataVectors3(nullptr),
			  _data1(nullptr),
			  _data2(nullptr),
			  _data3(nullptr) {

		assert(histTemplate != nullptr);
		assert(_vecData1 != nullptr);
		makePlot(cutmasks, histTemplate);

	};

	template<typename T1, typename T2, typename T3>
	antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
	                                            TH1 *histTemplate,
	                                            std::vector<T1> *vecData1,
	                                            std::vector<T2> *vecData2,
	                                            std::vector<T3> *vecData3)
			: Plot(),
			  _mode(2),
			  _vecData1(nullptr),
			  _vecData2(nullptr),
			  _vecData3(nullptr),
			  _vecDataVector1(vecData1),
			  _vecDataVector2(vecData2),
			  _vecDataVector3(vecData3),
			  _multipleVecDataVectors1(nullptr),
			  _multipleVecDataVectors2(nullptr),
			  _multipleVecDataVectors3(nullptr),
			  _data1(nullptr),
			  _data2(nullptr),
			  _data3(nullptr) {

		assert(histTemplate != nullptr);
		assert(_vecDataVector1 != nullptr);
		makePlot(cutmasks, histTemplate);

	};

	template<typename T1, typename T2, typename T3>
	antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
	                                            TH1 *histTemplate,
	                                            std::vector<std::vector<T1>*> *vecData1,
	                                            std::vector<std::vector<T2>*> *vecData2,
	                                            std::vector<std::vector<T3>*> *vecData3)
			: Plot(),
			  _mode(3),
			  _vecData1(nullptr),
			  _vecData2(nullptr),
			  _vecData3(nullptr),
			  _vecDataVector1(nullptr),
			  _vecDataVector2(nullptr),
			  _vecDataVector3(nullptr),
			  _multipleVecDataVectors1(vecData1),
			  _multipleVecDataVectors2(vecData2),
			  _multipleVecDataVectors3(vecData3),
			  _data1(nullptr),
			  _data2(nullptr),
			  _data3(nullptr) {

		assert(histTemplate != 0);
		assert(_multipleVecDataVectors1 != nullptr);
		makePlot(cutmasks, histTemplate);

	};

	template<typename T1, typename T2, typename T3>
	antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
	                                            TH1 *histTemplate,
	                                            T1* data1,
	                                            std::vector<T2> *vecData2)

			: Plot(),
			  _mode(4),
			  _vecData1(nullptr),
			  _vecData2(nullptr),
			  _vecData3(nullptr),
			  _vecDataVector1(nullptr),
			  _vecDataVector2(vecData2),
			  _vecDataVector3(nullptr),
			  _multipleVecDataVectors1(nullptr),
			  _multipleVecDataVectors2(nullptr),
			  _multipleVecDataVectors3(nullptr),
			  _data1(data1),
			  _data2(nullptr),
			  _data3(nullptr) {

		assert(histTemplate != nullptr);
		assert(_data1 != nullptr);
		assert(_vecDataVector2 != nullptr);
		makePlot(cutmasks, histTemplate);

	};

	template<typename T1, typename T2, typename T3>
	antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
	                                            TH1 *histTemplate,
	                                            std::vector<T1> *vecData1,
	                                            T2* data2)

			: Plot(),
			  _mode(5),
			  _vecData1(nullptr),
			  _vecData2(nullptr),
			  _vecData3(nullptr),
			  _vecDataVector1(vecData1),
			  _vecDataVector2(nullptr),
			  _vecDataVector3(nullptr),
			  _multipleVecDataVectors1(nullptr),
			  _multipleVecDataVectors2(nullptr),
			  _multipleVecDataVectors3(nullptr),
			  _data1(nullptr),
			  _data2(data2),
			  _data3(nullptr) {

		assert(histTemplate != nullptr);
		assert(_vecDataVector1 != nullptr);
		assert(_data2 != nullptr);
		makePlot(cutmasks, histTemplate);

	};

	template<typename T1, typename T2, typename T3>
	antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
	                                            TH1 *histTemplate,
	                                            std::vector<T1*> *vecData1,
	                                            std::vector<std::vector<T2>*> *vecData2)

			: Plot(),
			  _mode(6),
			  _vecData1(vecData1),
			  _vecData2(nullptr),
			  _vecData3(nullptr),
			  _vecDataVector1(nullptr),
			  _vecDataVector2(nullptr),
			  _vecDataVector3(nullptr),
			  _multipleVecDataVectors1(nullptr),
			  _multipleVecDataVectors2(vecData2),
			  _multipleVecDataVectors3(nullptr),
			  _data1(nullptr),
			  _data2(nullptr),
			  _data3(nullptr) {

		assert(histTemplate != nullptr);
		assert(_vecDataVector1 != nullptr);
		assert(_multipleVecDataVectors2 != nullptr);
		assert(_vecDataVector1->size() != _multipleVecDataVectors2->size());
		makePlot(cutmasks, histTemplate);

	};

	template<typename T1, typename T2, typename T3>
	antok::TemplatePlot<T1,T2,T3>::TemplatePlot(std::map<std::string, std::vector<long> > &cutmasks,
	                                            TH1 *histTemplate,
	                                            std::vector<std::vector<T1>*> *vecData1,
	                                            std::vector<T2*> *vecData2)


			: Plot(),
			  _mode(7),
			  _vecData1(nullptr),
			  _vecData2(vecData2),
			  _vecData3(nullptr),
			  _vecDataVector1(nullptr),
			  _vecDataVector2(nullptr),
			  _vecDataVector3(nullptr),
			  _multipleVecDataVectors1(vecData1),
			  _multipleVecDataVectors2(nullptr),
			  _multipleVecDataVectors3(nullptr),
			  _data1(nullptr),
			  _data2(nullptr),
			  _data3(nullptr) {

		assert(histTemplate != nullptr);
		assert(_multipleVecDataVectors1 != nullptr);
		assert(_vecDataVector2 != nullptr);
		assert(_multipleVecDataVectors1->size() != _vecDataVector2->size());
		makePlot(cutmasks, histTemplate);

	};

	template<typename T1, typename T2, typename T3>
	void antok::TemplatePlot<T1,T2,T3>::fill(long cutPattern) {

		for (unsigned int i = 0; i < _histograms.size(); ++i) {
			TH1 *hist = _histograms[i].first;
			long histMask = _histograms[i].second;
			if ((histMask & cutPattern) == histMask) {
				switch (_mode) {
					case 0:
						if (_data2 == nullptr) {
							hist->Fill(*_data1);
						} else if (_data3 == nullptr) {
							hist->Fill(*_data1, *_data2);
						} else {
							static_cast<TH3D *>(hist)->Fill(*_data1, *_data2, *_data3);
						}
						break;
					case 1:
						if( _vecData2 == nullptr ) {
							for (unsigned int j = 0; j < _vecData1->size(); ++j) {
								hist->Fill(*(*_vecData1)[j]);
							}
						} else if( _vecData3 == nullptr ) {
							assert(_vecData1->size() == _vecData2->size());
							for (unsigned int j = 0; j < _vecData1->size(); ++j) {
								hist->Fill(*(*_vecData1)[j], *(*_vecData2)[j]);
							}
						} else {
							assert(_vecData1->size() == _vecData2->size());
							assert(_vecData1->size() == _vecData3->size());
							for (unsigned int j = 0; j < _vecData1->size(); ++j) {
								for (unsigned int k = 0; k < _vecData2->size(); ++k) {
									static_cast<TH3D *>(hist)->Fill(*(*_vecData1)[j], *(*_vecData2)[j], *(*_vecData3)[j]);
								}
							}
						}
						break;
					case 2:
						if( _vecDataVector2 == nullptr ) {
							for (unsigned int j = 0; j < _vecData1->size(); ++j) {
								hist->Fill((*_vecDataVector1)[j]);
							}
						} else if( _vecDataVector3 == nullptr ) {
							assert(_vecDataVector1->size() == _vecDataVector2->size());
							for (unsigned int j = 0; j < _vecDataVector1->size(); ++j) {
								hist->Fill((*_vecDataVector1)[j], (*_vecDataVector2)[j]);
							}
						} else {
							assert(_vecDataVector1->size() == _vecDataVector2->size());
							assert(_vecDataVector1->size() == _vecDataVector3->size());
							for (unsigned int j = 0; j < _vecDataVector1->size(); ++j) {
								for (unsigned int k = 0; k < _vecDataVector2->size(); ++k) {
									static_cast<TH3D *>(hist)->Fill((*_vecDataVector1)[j], (*_vecDataVector2)[j], (*_vecDataVector3)[j]);
								}
							}
						}
						break;
					case 3:
						if( _multipleVecDataVectors2 == nullptr ) {
							for ( unsigned int j = 0; j < _multipleVecDataVectors1->size(); ++j ) {
								for(unsigned int k = 0; k < (*_multipleVecDataVectors1)[j]->size(); ++k ) {
									hist->Fill((*(*_multipleVecDataVectors1)[j])[k]);
								}
							}
						} else if( _multipleVecDataVectors3 == nullptr ) {
							assert(_multipleVecDataVectors1->size() == _multipleVecDataVectors2->size());
							for ( unsigned int j = 0; j < _multipleVecDataVectors1->size(); ++j ) {
								assert( (*_multipleVecDataVectors1)[j]->size() == (*_multipleVecDataVectors2)[j]->size() );
								for(unsigned int k = 0; k < (*_multipleVecDataVectors1)[j]->size(); ++k ) {
									hist->Fill((*(*_multipleVecDataVectors1)[j])[k],(*(*_multipleVecDataVectors2)[j])[k]);
								}
							}
						} else {
							assert(_multipleVecDataVectors1->size() == _multipleVecDataVectors2->size());
							assert(_multipleVecDataVectors1->size() == _multipleVecDataVectors3->size());
							for ( unsigned int j = 0; j < _multipleVecDataVectors1->size(); ++j ) {
								assert( (*_multipleVecDataVectors1)[j]->size() == (*_multipleVecDataVectors2)[j]->size() );
								assert( (*_multipleVecDataVectors1)[j]->size() == (*_multipleVecDataVectors3)[j]->size() );
								for(unsigned int k = 0; k < (*_multipleVecDataVectors1)[j]->size(); ++k ) {
									static_cast<TH3D *>(hist)->Fill((*(*_multipleVecDataVectors1)[j])[k],(*(*_multipleVecDataVectors2)[j])[k],(*(*_multipleVecDataVectors3)[j])[k]);
								}
							}
						}
						break;
					case 4:
						for( unsigned int j = 0; j < _vecDataVector2->size(); ++j ) {
							hist->Fill( *_data1, (*_vecDataVector2)[j] );
						}
						break;
					case 5:
						for( unsigned int j = 0; j < _vecDataVector1->size(); ++j ) {
							hist->Fill((*_vecDataVector1)[j], *_data2);
						}
						break;
					case 6:
						for( unsigned int j = 0; j < _vecDataVector1->size(); ++j ) {
							for( unsigned int k = 0; k < (*_multipleVecDataVectors2)[j]->size(); ++k ) {
								hist->Fill( *(*_vecData1)[j], (*(*_multipleVecDataVectors2)[j])[k]);
							}
						}
						break;
					case 7:
						for( unsigned int j = 0; j < _vecDataVector2->size(); ++j ) {
							for( unsigned int k = 0; k < (*_multipleVecDataVectors1)[j]->size(); ++k ) {
								hist->Fill((*(*_multipleVecDataVectors1)[j])[k],  *(*_vecData2)[j]);
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

		for (std::map<std::string, std::vector<long> >::const_iterator cutmasks_it = cutmasks.begin();
		     cutmasks_it != cutmasks.end(); ++cutmasks_it) {
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
}

#endif
