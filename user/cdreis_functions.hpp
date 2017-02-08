#ifndef ANTOK_USER_CDREIS_FUNCTIONS_HPP
#define ANTOK_USER_CDREIS_FUNCTIONS_HPP

#include<iostream>
#include<vector>


#include<TLorentzVector.h>
#include<TRotation.h>
#include<TLorentzRotation.h>

namespace antok {

	namespace user {

		namespace cdreis {

			namespace functions {

				class GetRecoilLorentzVec : public Function {
				public:
					GetRecoilLorentzVec(TLorentzVector *BeamLorentzVec,
					                    TLorentzVector *XLorentzVec,
					                    double *RecoilMass,
					                    TLorentzVector *RecoilLorentzVec)
							: _BeamLorentzVec(BeamLorentzVec),
							  _XLorentzVec(XLorentzVec),
							  _RecoilMass(RecoilMass),
							  _RecoilLorentzVec(RecoilLorentzVec) {}

					virtual ~GetRecoilLorentzVec() {}

					bool operator()() {
						_RecoilLorentzVec->SetVectM(((_BeamLorentzVec->Vect()) - (_XLorentzVec->Vect())), *_RecoilMass);
						return true;
					}

				private:
					TLorentzVector *_BeamLorentzVec;
					TLorentzVector *_XLorentzVec;
					double *_RecoilMass;
					TLorentzVector *_RecoilLorentzVec;
				};


				class GetPhotonLorentzVecs : public Function {
				public:
					GetPhotonLorentzVecs(std::vector<double> *VectorXAddr,
					                     std::vector<double> *VectorYAddr,
					                     std::vector<double> *VectorZAddr,
					                     std::vector<double> *VectorEAddr,
					                     std::vector<double> *VectorTimeAddr,
					                     double *xPVAddr,
					                     double *yPVAddr,
					                     double *zPVAddr,
					                     std::vector<TLorentzVector> *resultVec,
					                     std::vector<int> *resultECALIndex)
							: _VectorXAddr(VectorXAddr),
							  _VectorYAddr(VectorYAddr),
							  _VectorZAddr(VectorZAddr),
							  _VectorEAddr(VectorEAddr),
							  _VectorTimeAddr(VectorTimeAddr),
							  _xPVAddr(xPVAddr),
							  _yPVAddr(yPVAddr),
							  _zPVAddr(zPVAddr),
							  _resultVec(resultVec),
							  _resultECALIndex(resultECALIndex) {}

					virtual ~GetPhotonLorentzVecs() {}

					bool operator()() {
						(*_resultVec).reserve(_VectorXAddr->size());
						(*_resultVec).clear();
						(*_resultECALIndex).reserve(_VectorXAddr->size());
						(*_resultECALIndex).clear();

						for (unsigned int i = 0; i < _VectorXAddr->size(); i++) {
							TVector3 v3(((*_VectorXAddr)[i] - (*_xPVAddr)), ((*_VectorYAddr)[i] - (*_yPVAddr)),
							            ((*_VectorZAddr)[i] - (*_zPVAddr)));
							v3.SetMag((*_VectorEAddr)[i]);
							TLorentzVector lv;
							lv.SetXYZT(v3.X(), v3.Y(), v3.Z(), (*_VectorEAddr)[i]);
							(*_resultVec).push_back(lv);
							// ECAL 1
							if ((*_VectorZAddr)[i] < 2500) {
								(*_resultECALIndex).push_back(1);
							}
								// ECAL 2
							else {
								(*_resultECALIndex).push_back(2);
							}
						}
						return true;
					}

				private:
					std::vector<double> *_VectorXAddr;
					std::vector<double> *_VectorYAddr;
					std::vector<double> *_VectorZAddr;
					std::vector<double> *_VectorEAddr;
					std::vector<double> *_VectorTimeAddr;
					double *_xPVAddr;
					double *_yPVAddr;
					double *_zPVAddr;
					std::vector<TLorentzVector> *_resultVec;
					std::vector<int> *_resultECALIndex;
				};


				class GetVectorLorentzVectorAttributes : public Function {
				public:
					GetVectorLorentzVectorAttributes(std::vector<TLorentzVector> *VectorLV,
					                                 std::vector<double> *resultVecX,
					                                 std::vector<double> *resultVecY,
					                                 std::vector<double> *resultVecZ,
					                                 std::vector<double> *resultVecE,
					                                 std::vector<double> *resultVecTheta,
					                                 std::vector<double> *resultVecPhi,
					                                 std::vector<double> *resultVecMag)
							: _VectorLV(VectorLV),
							  _resultVecX(resultVecX),
							  _resultVecY(resultVecY),
							  _resultVecZ(resultVecZ),
							  _resultVecE(resultVecE),
							  _resultVecTheta(resultVecTheta),
							  _resultVecPhi(resultVecPhi),
							  _resultVecMag(resultVecMag) {}

					virtual ~GetVectorLorentzVectorAttributes() {}

					bool operator()() {
						(*_resultVecX).reserve(_VectorLV->size());
						(*_resultVecX).clear();
						(*_resultVecY).reserve(_VectorLV->size());
						(*_resultVecY).clear();
						(*_resultVecZ).reserve(_VectorLV->size());
						(*_resultVecZ).clear();
						(*_resultVecE).reserve(_VectorLV->size());
						(*_resultVecE).clear();
						(*_resultVecTheta).reserve(_VectorLV->size());
						(*_resultVecTheta).clear();
						(*_resultVecPhi).reserve(_VectorLV->size());
						(*_resultVecPhi).clear();
						(*_resultVecMag).reserve(_VectorLV->size());
						(*_resultVecMag).clear();
						for (unsigned int i = 0; i < _VectorLV->size(); i++) {
							_resultVecX->push_back((*_VectorLV)[i].X());
							_resultVecY->push_back((*_VectorLV)[i].Y());
							_resultVecZ->push_back((*_VectorLV)[i].Z());
							_resultVecE->push_back((*_VectorLV)[i].E());
							_resultVecTheta->push_back((*_VectorLV)[i].Theta());
							_resultVecPhi->push_back((*_VectorLV)[i].Phi());
							_resultVecMag->push_back((*_VectorLV)[i].Mag());
						}
						return true;
					}

				private:
					std::vector<TLorentzVector> *_VectorLV;
					std::vector<double> *_resultVecX;
					std::vector<double> *_resultVecY;
					std::vector<double> *_resultVecZ;
					std::vector<double> *_resultVecE;
					std::vector<double> *_resultVecTheta;
					std::vector<double> *_resultVecPhi;
					std::vector<double> *_resultVecMag;
				};

				class GetPi0s : public Function {
				public:
					GetPi0s(std::vector<TLorentzVector> *VectorLV,
					        std::vector<int> *ECALIndex,
					        double *mass,
					        std::vector<TLorentzVector> *resultVecLV)
							: _VectorLV(VectorLV),
							  _ECALIndex(ECALIndex),
							  _mass(mass),
							  _resultVecLV(resultVecLV) {}

					virtual ~GetPi0s() {}

					bool operator()() {
						(*_resultVecLV).reserve(_VectorLV->size() * (_VectorLV->size()) * 100);
						(*_resultVecLV).clear();
						if (_VectorLV->size() < 2) {
							return true;
						}
						for (unsigned int i = 0; i < _VectorLV->size(); i++) {
							TLorentzVector comparator = _VectorLV->back();
							int ECALIndexComparator = (*_ECALIndex)[i];
							double massResolution;
							_VectorLV->pop_back();
							_ECALIndex->pop_back();
							for (unsigned int j = 0; j < _VectorLV->size(); j++) {
								TLorentzVector lv = (*_VectorLV)[j];
								TLorentzVector pi0Candidate = comparator + lv;
								int ECALIndexLV = (*_ECALIndex)[j];
								// ECAL 1
								if (ECALIndexComparator == 1 && ECALIndexLV == 1) {
									massResolution = 3. * 0.028819;
								}
								// ECAL 2
								else if (ECALIndexComparator == 2 && ECALIndexLV == 2) {
									massResolution = 3. * 0.0151949;
								}
								// ECAL 1 + 2 mixed
								else {
									massResolution = 3. * 0.0197151;
								}
								if (std::fabs(pi0Candidate.Mag() - *_mass) < massResolution) {

									_resultVecLV->push_back(pi0Candidate);
								}
							}
						}
						return true;
					}

				private:
					std::vector<TLorentzVector> *_VectorLV;
					std::vector<int> *_ECALIndex;
					double *_mass;
					std::vector<TLorentzVector> *_resultVecLV;
				};

				class GetCleanedEcalClusters : public Function {
				public:
					GetCleanedEcalClusters(std::vector<double> *VectorXAddr,
					                       std::vector<double> *VectorYAddr,
					                       std::vector<double> *VectorZAddr,
					                       std::vector<double> *VectorEAddr,
					                       std::vector<double> *VectorTAddr,
					                       std::vector<double> *resultVectorX,
					                       std::vector<double> *resultVectorY,
					                       std::vector<double> *resultVectorZ,
					                       std::vector<double> *resultVectorE,
					                       std::vector<double> *resultVectorT,
					                       std::vector<int> *resultVectorIndex)
							: _VectorXAddr(VectorXAddr),
							  _VectorYAddr(VectorYAddr),
							  _VectorZAddr(VectorZAddr),
							  _VectorEAddr(VectorEAddr),
							  _VectorTAddr(VectorTAddr),
							  _resultVectorXAddr(resultVectorX),
							  _resultVectorYAddr(resultVectorY),
							  _resultVectorZAddr(resultVectorZ),
							  _resultVectorEAddr(resultVectorE),
							  _resultVectorTAddr(resultVectorT),
							  _resultVectorIndexAddr(resultVectorIndex) {}

					virtual ~GetCleanedEcalClusters() {}

					bool operator()() {
						(*_resultVectorXAddr).reserve(_VectorXAddr->size());
						(*_resultVectorXAddr).clear();
						(*_resultVectorYAddr).reserve(_VectorXAddr->size());
						(*_resultVectorYAddr).clear();
						(*_resultVectorZAddr).reserve(_VectorXAddr->size());
						(*_resultVectorZAddr).clear();
						(*_resultVectorEAddr).reserve(_VectorXAddr->size());
						(*_resultVectorEAddr).clear();
						(*_resultVectorTAddr).reserve(_VectorXAddr->size());
						(*_resultVectorTAddr).clear();
						(*_resultVectorIndexAddr).reserve(_VectorXAddr->size());
						(*_resultVectorIndexAddr).clear();

						// Energy threshold and timing
						for (unsigned int i = 0; i < _VectorXAddr->size(); i++) {
							// ECAL 1
							if ((*_VectorZAddr)[i] < 2500) {
								if ((*_VectorEAddr)[i] < 0.6 || fabs((*_VectorTAddr)[i]) > 3.75) {
									continue;
								}
								(*_resultVectorIndexAddr).push_back(1);
							}
								// ECAL 2
							else {
								if ((*_VectorEAddr)[i] < 1.2 || fabs((*_VectorTAddr)[i]) > 3.00) {
									continue;
								}
								(*_resultVectorIndexAddr).push_back(2);
							}

							(*_resultVectorXAddr).push_back((*_VectorXAddr)[i]);
							(*_resultVectorYAddr).push_back((*_VectorYAddr)[i]);
							(*_resultVectorZAddr).push_back((*_VectorZAddr)[i]);
							(*_resultVectorEAddr).push_back((*_VectorEAddr)[i]);
							(*_resultVectorTAddr).push_back((*_VectorTAddr)[i]);
						}
						return true;
					}

				private:
					std::vector<double> *_VectorXAddr;
					std::vector<double> *_VectorYAddr;
					std::vector<double> *_VectorZAddr;
					std::vector<double> *_VectorEAddr;
					std::vector<double> *_VectorTAddr;
					std::vector<double> *_resultVectorXAddr;
					std::vector<double> *_resultVectorYAddr;
					std::vector<double> *_resultVectorZAddr;
					std::vector<double> *_resultVectorEAddr;
					std::vector<double> *_resultVectorTAddr;
					std::vector<int> *_resultVectorIndexAddr;
				};

				class GetPi0Pair : public Function {
				public:
					GetPi0Pair(std::vector<TLorentzVector> *VectorPhotonLV,
					           std::vector<TLorentzVector> *resultVecLV,
					           TLorentzVector *resultVecLV0,
					           TLorentzVector *resultVecLV1,
					           int *resultGoodPair)
							: _VectorPhotonLV(VectorPhotonLV),
							  _resultVecLV(resultVecLV),
							  _resultVecLV0(resultVecLV0),
							  _resultVecLV1(resultVecLV1),
							  _resultGoodPair(resultGoodPair) {}

					virtual ~GetPi0Pair() {}

					bool operator()() {
						(*_resultGoodPair) = 0;
						(*_resultVecLV).reserve(2);
						(*_resultVecLV).clear();
						if (_VectorPhotonLV->size() < 2) {
							return true;
						}

						unsigned int numberCandidatePairs = 0;
						TLorentzVector pi0Candidate0(0, 0, 0, 0);
						TLorentzVector pi0Candidate1(0, 0, 0, 0);
						for (unsigned int i = 0; i < _VectorPhotonLV->size(); i++) {
							if (numberCandidatePairs > 1) {
								break;
							}
							for (unsigned int j = i + 1; j < _VectorPhotonLV->size(); j++) {
								if (numberCandidatePairs > 1) {
									break;
								}
								TLorentzVector photon0 = (*_VectorPhotonLV)[i];
								TLorentzVector photon1 = (*_VectorPhotonLV)[j];
								pi0Candidate0 = photon0 + photon1;
								if (std::fabs(pi0Candidate0.Mag() - 0.13957018) > 0.06) {
									continue;
								}

								for (unsigned int m = i + 1; m < _VectorPhotonLV->size(); m++) {
									if (numberCandidatePairs > 1) {
										break;
									}
									for (unsigned int n = m + 1; n < _VectorPhotonLV->size(); n++) {
										if (numberCandidatePairs > 1) {
											break;
										}
										if (m == j || n == j) {
											continue;
										}
										TLorentzVector photon2 = (*_VectorPhotonLV)[m];
										TLorentzVector photon3 = (*_VectorPhotonLV)[n];
										pi0Candidate1 = photon2 + photon3;
										if (std::fabs(pi0Candidate1.Mag() - 0.13957018) > 0.06) {
											continue;
										}
										numberCandidatePairs++;
									}
								}
							}
						}
						if (numberCandidatePairs == 1) {
							(*_resultVecLV).push_back(pi0Candidate0);
							(*_resultVecLV).push_back(pi0Candidate1);
						}

						if (_resultVecLV->size() == 2) {
							(*_resultGoodPair) = 1;
							TLorentzVector both = pi0Candidate0 + pi0Candidate1;
						}

						return true;
					}

				private:
					std::vector<TLorentzVector> *_VectorPhotonLV;
					std::vector<TLorentzVector> *_resultVecLV;
					TLorentzVector *_resultVecLV0;
					TLorentzVector *_resultVecLV1;
					int *_resultGoodPair;
				};

			}

		}

	}

}

#endif