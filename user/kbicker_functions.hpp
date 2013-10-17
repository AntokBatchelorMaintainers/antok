#ifndef ANTOK_USER_KBICKER_FUNCTIONS_HPP
#define ANTOK_USER_KBICKER_FUNCTIONS_HPP

#include<iostream>
#include<vector>

#include<TLorentzVector.h>

#include<kbicker.h>

namespace antok {

	namespace user {

		namespace kbicker {

			namespace functions {

				class GetRpdPhi : public Function
				{

				  public:

					GetRpdPhi(TLorentzVector* beamLorentzVecAddr,
							  TLorentzVector* rpdProtonLorentzVecAddr,
							  TLorentzVector* xLorentzVecAddr,
							  double* rpdDeltaPhiAddr,
							  double* rpdDeltaPhiResAddr,
							  int method,
							  double* protonPhiAddr = 0,
							  double* xPhiAddr = 0,
							  TVector3* vertex = 0
					)
						: _beamLorentzVecAddr(beamLorentzVecAddr),
						  _rpdProtonLorentzVecAddr(rpdProtonLorentzVecAddr),
						  _xLorentzVecAddr(xLorentzVecAddr),
						  _vertex(vertex),
						  _rpdDeltaPhiAddr(rpdDeltaPhiAddr),
						  _rpdDeltaPhiResAddr(rpdDeltaPhiResAddr),
						  _protonPhiAddr(protonPhiAddr),
						  _xPhiAddr(xPhiAddr),
						  _method(method) { }

					virtual ~GetRpdPhi() { }

					bool operator() () {
						switch(_method)
						{
							case 0:
								antok::user::kbicker::getRPDDeltaPhiResProjection((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr));
								return true;
							case 1:
								if(_protonPhiAddr == 0 and _xPhiAddr == 0) {
									antok::user::kbicker::getRPDDeltaPhiResRotation((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr));
								} else {
									antok::user::kbicker::getRPDDeltaPhiResRotation((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr), (*_protonPhiAddr), (*_xPhiAddr));
								}
								return true;
							case 2:
								if(not _vertex) {
									std::cerr<<"Found null pointer for vertex in GetRpdPhi while using prediction method."<<std::endl;
									return false;
								}
								if(_protonPhiAddr == 0 and _xPhiAddr == 0) {
									antok::user::kbicker::getRPDDeltaPhiResPrediction((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_vertex), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr));
								} else {
									antok::user::kbicker::getRPDDeltaPhiResPrediction((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_vertex), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr), (*_protonPhiAddr), (*_xPhiAddr));
								}
								return true;
						}
						return false;
					};

				  private:

					TLorentzVector* _beamLorentzVecAddr;
					TLorentzVector* _rpdProtonLorentzVecAddr;
					TLorentzVector* _xLorentzVecAddr;
					TVector3* _vertex;
					double* _rpdDeltaPhiAddr;
					double* _rpdDeltaPhiResAddr;
					double* _protonPhiAddr;
					double* _xPhiAddr;
					int _method;

				};

				class GetRPDExpectedHitsParameters : public Function
				{

				  public:

					GetRPDExpectedHitsParameters(TLorentzVector* pBeamAddr,
												 TLorentzVector* pXAddr,
												 TVector3* vertexAddr,
												 double* xOffsetAddr,
												 double* yOffsetAddr,
												 double* xAngleAddr,
												 double* yAngleAddr,
												 double* rpdPhiRingAAddr,
												 double* rpdPhiRingBAddr,
												 double* rpdZRingAAddr,
												 double* rpdZRingBAddr)
						: _pBeamAddr(pBeamAddr),
						  _pXAddr(pXAddr),
						  _vertexAddr(vertexAddr),
						  _xOffsetAddr(xOffsetAddr),
						  _yOffsetAddr(yOffsetAddr),
						  _xAngleAddr(xAngleAddr),
						  _yAngleAddr(yAngleAddr),
						  _rpdPhiRingAAddr(rpdPhiRingAAddr),
						  _rpdPhiRingBAddr(rpdPhiRingBAddr),
						  _rpdZRingAAddr(rpdZRingAAddr),
						  _rpdZRingBAddr(rpdZRingBAddr) { }

					virtual ~GetRPDExpectedHitsParameters() { }

					bool operator() () {
						antok::user::kbicker::getRPDExpectedHitsParameters(*_pBeamAddr,
						                                                   *_pXAddr,
						                                                   *_vertexAddr,
						                                                   *_xOffsetAddr,
						                                                   *_yOffsetAddr,
						                                                   *_xAngleAddr,
						                                                   *_yAngleAddr,
						                                                   *_rpdPhiRingAAddr,
						                                                   *_rpdPhiRingBAddr,
						                                                   *_rpdZRingAAddr,
						                                                   *_rpdZRingBAddr);
						return true;
					}

				  private:

					TLorentzVector* _pBeamAddr;
					TLorentzVector* _pXAddr;
					TVector3* _vertexAddr;
					double* _xOffsetAddr;
					double* _yOffsetAddr;
					double* _xAngleAddr;
					double* _yAngleAddr;
					double* _rpdPhiRingAAddr;
					double* _rpdPhiRingBAddr;
					double* _rpdZRingAAddr;
					double* _rpdZRingBAddr;

				};

				class GetCutOnExtraTracks : public Function
				{

				  public:

					GetCutOnExtraTracks(std::vector<double>* trackTimesAddr,
					                    std::vector<double>* trackTimeSigmasAddr,
					                    std::vector<double>* trackNHitsAddr,
					                    std::vector<double>* trackZFirstAddr,
					                    std::vector<double>* trackZLastAddr,
					                    std::vector<double>* trackQAddr,
					                    double* cutAddr,
					                    int* nExtraTracksAddr,
					                    double* xTracksMeanTime)
						: _nParticles(antok::Constants::nParticles()),
						  _trackTimesAddr(trackTimesAddr),
						  _trackTimeSigmasAddr(trackTimeSigmasAddr),
						  _trackNHitsAddr(trackNHitsAddr),
						  _trackZFirstAddr(trackZFirstAddr),
						  _trackZLastAddr(trackZLastAddr),
						  _trackQAddr(trackQAddr),
						  _cutAddr(cutAddr),
						  _nExtraTracksAddr(nExtraTracksAddr),
						  _xTracksMeanTime(xTracksMeanTime) { }

					virtual ~GetCutOnExtraTracks() { }

					bool operator() () {
						// Do the cut
						bool cut  = antok::user::kbicker::extraTracksCut(*_trackTimesAddr,
						                                                 *_trackTimeSigmasAddr,
						                                                 *_trackNHitsAddr,
						                                                 *_trackZFirstAddr,
						                                                 *_trackZLastAddr,
						                                                 *_trackQAddr);
						(*_cutAddr) = cut ? 1. : 0.;
						// Number of extra tracks
						(*_nExtraTracksAddr) = (int)(_trackTimesAddr->size() - _nParticles);
						// Mean time of track from X
						(*_xTracksMeanTime) = 0.;
						double weightSum = 0.;
						for(unsigned int i = 0; i < _nParticles; ++i) {
							double weight = 1. / ((*_trackTimeSigmasAddr)[i] * (*_trackTimeSigmasAddr)[i]);
							(*_xTracksMeanTime) += (*_trackTimesAddr)[i] * weight;
							weightSum += weight;
						}
						(*_xTracksMeanTime) /= weightSum;
						return true;
					}

					  private:

						const unsigned int& _nParticles;

						std::vector<double>* _trackTimesAddr;
						std::vector<double>* _trackTimeSigmasAddr;
						std::vector<double>* _trackNHitsAddr;
						std::vector<double>* _trackZFirstAddr;
						std::vector<double>* _trackZLastAddr;
						std::vector<double>* _trackQAddr;
						double* _cutAddr;
						int* _nExtraTracksAddr;
						double* _xTracksMeanTime;

				};

			}

		}

	}

}

#endif
