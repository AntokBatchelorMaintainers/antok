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


				class GetPhotonLorentzVecs : public Function
				{
				public:
					GetPhotonLorentzVecs(std::vector<double>         *VectorXAddr,
					                     std::vector<double>         *VectorYAddr,
					                     std::vector<double>         *VectorZAddr,
					                     std::vector<double>         *VectorEAddr,
					                     std::vector<double>         *VectorTimeAddr,
					                     double                      *xPVAddr,
					                     double                      *yPVAddr,
					                     double                      *zPVAddr,
					                     double                      *RangeECAL1,
					                     std::vector<TLorentzVector> *resultVec,
					                     std::vector<int>            *resultECALIndex)
							: _VectorXAddr    ( VectorXAddr     ),
							  _VectorYAddr    ( VectorYAddr     ),
							  _VectorZAddr    ( VectorZAddr     ),
							  _VectorEAddr    ( VectorEAddr     ),
							  _VectorTimeAddr ( VectorTimeAddr  ),
							  _xPVAddr        ( xPVAddr         ),
							  _yPVAddr        ( yPVAddr         ),
							  _zPVAddr        ( zPVAddr         ),
							  _RangeECAL1     ( RangeECAL1      ),
							  _resultVec      ( resultVec       ),
							  _resultECALIndex( resultECALIndex ) {}

					virtual ~GetPhotonLorentzVecs() {}

					bool operator()()
					{
						_resultVec->reserve(_VectorXAddr->size());
						_resultVec->clear();
						_resultECALIndex->reserve(_VectorXAddr->size());
						_resultECALIndex->clear();

						for (unsigned int i = 0; i < _VectorXAddr->size(); i++)
						{
							TVector3 v3(((*_VectorXAddr)[i] - (*_xPVAddr)), ((*_VectorYAddr)[i] - (*_yPVAddr)),
							            ((*_VectorZAddr)[i] - (*_zPVAddr)));
							v3.SetMag((*_VectorEAddr)[i]);
							TLorentzVector lv;
							lv.SetXYZT(v3.X(), v3.Y(), v3.Z(), (*_VectorEAddr)[i]);
							_resultVec->push_back(lv);
							if ((*_VectorZAddr)[i] < *_RangeECAL1 )
							{
								_resultECALIndex->push_back(1);
							}
							else
							{
								_resultECALIndex->push_back(2);
							}
						}
						return true;
					}

				private:
					std::vector<double>         *_VectorXAddr;
					std::vector<double>         *_VectorYAddr;
					std::vector<double>         *_VectorZAddr;
					std::vector<double>         *_VectorEAddr;
					std::vector<double>         *_VectorTimeAddr;
					double                      *_xPVAddr;
					double                      *_yPVAddr;
					double                      *_zPVAddr;
					double                      *_RangeECAL1;
					std::vector<TLorentzVector> *_resultVec;
					std::vector<int>            *_resultECALIndex;
				};


				class GetVectorLorentzVectorAttributes : public Function
				{
				public:
					GetVectorLorentzVectorAttributes( std::vector<TLorentzVector> *VectorLV,
					                                  std::vector<double>         *resultVecX,
					                                  std::vector<double>         *resultVecY,
					                                  std::vector<double>         *resultVecZ,
					                                  std::vector<double>         *resultVecE,
					                                  std::vector<double>         *resultVecTheta,
					                                  std::vector<double>         *resultVecPhi,
					                                  std::vector<double>         *resultVecMag )
							: _VectorLV      ( VectorLV       ),
							  _resultVecX    ( resultVecX     ),
							  _resultVecY    ( resultVecY     ),
							  _resultVecZ    ( resultVecZ     ),
							  _resultVecE    ( resultVecE     ),
							  _resultVecTheta( resultVecTheta ),
							  _resultVecPhi  ( resultVecPhi   ),
							  _resultVecMag  ( resultVecMag   ) {}

					virtual ~GetVectorLorentzVectorAttributes() {}

					bool operator()()
					{
						_resultVecX->reserve(_VectorLV->size());
						_resultVecX->clear();
						_resultVecY->reserve(_VectorLV->size());
						_resultVecY->clear();
						_resultVecZ->reserve(_VectorLV->size());
						_resultVecZ->clear();
						_resultVecE->reserve(_VectorLV->size());
						_resultVecE->clear();
						_resultVecTheta->reserve(_VectorLV->size());
						_resultVecTheta->clear();
						_resultVecPhi->reserve(_VectorLV->size());
						_resultVecPhi->clear();
						_resultVecMag->reserve(_VectorLV->size());
						_resultVecMag->clear();
						for( unsigned int i = 0; i < _VectorLV->size(); i++)
						{
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
					std::vector<double>         *_resultVecX;
					std::vector<double>         *_resultVecY;
					std::vector<double>         *_resultVecZ;
					std::vector<double>         *_resultVecE;
					std::vector<double>         *_resultVecTheta;
					std::vector<double>         *_resultVecPhi;
					std::vector<double>         *_resultVecMag;
				};

				class GetCleanedEcalClusters : public Function
				{
				public:
					GetCleanedEcalClusters(std::vector<double> *VectorXAddr,
					                       std::vector<double> *VectorYAddr,
					                       std::vector<double> *VectorZAddr,
					                       std::vector<double> *VectorEAddr,
					                       std::vector<double> *VectorTAddr,
					                       double              *RangeECAL1,
					                       double              *ThresholdEnergyECAL1,
					                       double              *ThresholdTimingECAL1,
					                       double              *ThresholdEnergyECAL2,
					                       double              *ThresholdTimingECAL2,
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
							  _RangeECAL1(RangeECAL1),
							  _ThresholdEnergyECAL1(ThresholdEnergyECAL1),
							  _ThresholdTimingECAL1(ThresholdTimingECAL1),
							  _ThresholdEnergyECAL2(ThresholdEnergyECAL2),
							  _ThresholdTimingECAL2(ThresholdTimingECAL2),
							  _resultVectorXAddr(resultVectorX),
							  _resultVectorYAddr(resultVectorY),
							  _resultVectorZAddr(resultVectorZ),
							  _resultVectorEAddr(resultVectorE),
							  _resultVectorTAddr(resultVectorT),
							  _resultVectorIndexAddr(resultVectorIndex) {}

					virtual ~GetCleanedEcalClusters() {}

					bool operator()()
					{
						_resultVectorXAddr->reserve(_VectorXAddr->size());
						_resultVectorXAddr->clear();
						_resultVectorYAddr->reserve(_VectorXAddr->size());
						_resultVectorYAddr->clear();
						_resultVectorZAddr->reserve(_VectorXAddr->size());
						_resultVectorZAddr->clear();
						_resultVectorEAddr->reserve(_VectorXAddr->size());
						_resultVectorEAddr->clear();
						_resultVectorTAddr->reserve(_VectorXAddr->size());
						_resultVectorTAddr->clear();
						_resultVectorIndexAddr->reserve(_VectorXAddr->size());
						_resultVectorIndexAddr->clear();

						for (unsigned int i = 0; i < _VectorXAddr->size(); i++)
						{
							if ((*_VectorZAddr)[i] < *_RangeECAL1)
							{
								if ((*_VectorEAddr)[i] < *_ThresholdEnergyECAL1 || fabs((*_VectorTAddr)[i]) > *_ThresholdTimingECAL1)
								{
									continue;
								}
								_resultVectorIndexAddr->push_back(1);
							}
							else
							{
								if ((*_VectorEAddr)[i] < *_ThresholdEnergyECAL2 || fabs((*_VectorTAddr)[i]) > *_ThresholdTimingECAL2)
								{
									continue;
								}
								_resultVectorIndexAddr->push_back(2);
							}

							_resultVectorXAddr->push_back((*_VectorXAddr)[i]);
							_resultVectorYAddr->push_back((*_VectorYAddr)[i]);
							_resultVectorZAddr->push_back((*_VectorZAddr)[i]);
							_resultVectorEAddr->push_back((*_VectorEAddr)[i]);
							_resultVectorTAddr->push_back((*_VectorTAddr)[i]);
						}
						return true;
					}

				private:
					std::vector<double> *_VectorXAddr;
					std::vector<double> *_VectorYAddr;
					std::vector<double> *_VectorZAddr;
					std::vector<double> *_VectorEAddr;
					std::vector<double> *_VectorTAddr;
					double              *_RangeECAL1;
					double              *_ThresholdEnergyECAL1;
					double              *_ThresholdTimingECAL1;
					double              *_ThresholdEnergyECAL2;
					double              *_ThresholdTimingECAL2;
					std::vector<double> *_resultVectorXAddr;
					std::vector<double> *_resultVectorYAddr;
					std::vector<double> *_resultVectorZAddr;
					std::vector<double> *_resultVectorEAddr;
					std::vector<double> *_resultVectorTAddr;
					std::vector<int>    *_resultVectorIndexAddr;
				};

				class GetPi0Pair : public Function
				{
				public:
					GetPi0Pair(std::vector<TLorentzVector> *VectorPhotonLV,
					           std::vector<int>            *ECALIndex,
					           double                      *Mass,
					           double                      *ResolutionECAL,
					           double                      *ResolutionECAL1,
					           double                      *ResolutionECAL2,
					           std::vector<TLorentzVector> *resultVecLV,
					           TLorentzVector              *resultVecLV0,
					           TLorentzVector              *resultVecLV1,
					           int                         *resultGoodPair)
							: _VectorPhotonLV( VectorPhotonLV   ),
							  _Mass           ( Mass            ),
							  _ResolutionECAL ( ResolutionECAL  ),
							  _ResolutionECAL1( ResolutionECAL1 ),
							  _ResolutionECAL2( ResolutionECAL2 ),
							  _ECALIndex      ( ECALIndex       ),
							  _resultVecLV    ( resultVecLV     ),
							  _resultVecLV0   ( resultVecLV0    ),
							  _resultVecLV1   ( resultVecLV1    ),
							  _resultGoodPair ( resultGoodPair  ) {}

					virtual ~GetPi0Pair() {}

					bool operator()()
					{
						(*_resultGoodPair) = 0;
						_resultVecLV->reserve(2);
						_resultVecLV->clear();
						if( _VectorPhotonLV->size() < 2 )
						{
							return true;
						}

						unsigned int numberCandidatePairs = 0;
						double resolution;
						TLorentzVector pi0Candidate0(0, 0, 0, 0);
						TLorentzVector pi0Candidate1(0, 0, 0, 0);
						for( unsigned int i = 0; i < _VectorPhotonLV->size(); i++ )
						{
							if( numberCandidatePairs > 1 )
							{
								break;
							}
							for( unsigned int j = i + 1; j < _VectorPhotonLV->size(); j++)
							{
								if(numberCandidatePairs > 1 )
								{
									break;
								}
								pi0Candidate0 = (*_VectorPhotonLV)[i] + (*_VectorPhotonLV)[j];
								if( (*_ECALIndex)[i] == 1 && (*_ECALIndex)[j] == 1 )
								{
									resolution = *_ResolutionECAL1;
								}
								else if( (*_ECALIndex)[i] == 2 && (*_ECALIndex)[j] == 2 )
								{
									resolution = *_ResolutionECAL2;
								}
								else if( (*_ECALIndex)[i] != (*_ECALIndex)[j] )
								{
									resolution = *_ResolutionECAL;
								}
								if (std::fabs(pi0Candidate0.Mag() - *_Mass) > resolution)
								{
									continue;
								}
								for (unsigned int m = i + 1; m < _VectorPhotonLV->size(); m++)
								{
									if (numberCandidatePairs > 1)
									{
										break;
									}
									for (unsigned int n = m + 1; n < _VectorPhotonLV->size(); n++)
									{
										if (numberCandidatePairs > 1)
										{
											break;
										}
										if (m == j || n == j)
										{
											continue;
										}
										pi0Candidate1 = (*_VectorPhotonLV)[m] + (*_VectorPhotonLV)[n];
										if( (*_ECALIndex)[m] == 1 && (*_ECALIndex)[n] == 1 )
										{
											resolution = *_ResolutionECAL1;
										}
										else if( (*_ECALIndex)[m] == 2 && (*_ECALIndex)[n] == 2 )
										{
											resolution = *_ResolutionECAL2;
										}
										else if( (*_ECALIndex)[m] != (*_ECALIndex)[n] )
										{
											resolution = *_ResolutionECAL;
										}
										if (std::fabs(pi0Candidate1.Mag() - *_Mass) > resolution)
										{
											continue;
										}
										if( numberCandidatePairs == 0 )
										{
											(*_resultVecLV).push_back(pi0Candidate0);
											(*_resultVecLV).push_back(pi0Candidate1);
											(*_resultVecLV0) = pi0Candidate0;
											(*_resultVecLV1) = pi0Candidate1;
										}
										numberCandidatePairs++;
									}
								}
							}
						}

						if ( numberCandidatePairs == 1)
						{
							(*_resultGoodPair) = 1;
						}

						return true;
					}

				private:
					std::vector<TLorentzVector> *_VectorPhotonLV;
					std::vector<int>            *_ECALIndex;
					double                      *_Mass;
					double                      *_ResolutionECAL;
					double                      *_ResolutionECAL1;
					double                      *_ResolutionECAL2;
					std::vector<TLorentzVector> *_resultVecLV;
					TLorentzVector              *_resultVecLV0;
					TLorentzVector              *_resultVecLV1;
					int                         *_resultGoodPair;
				};

				class GetOmega : public Function
				{
				public:
					GetOmega( TLorentzVector* Pi0_OAddr,
					          TLorentzVector* Pi0_1Addr,
					          TLorentzVector* Scattered0Addr,
					          TLorentzVector* Scattered1Addr,
					          TLorentzVector* Scattered2Addr,
					          int*            Charge0Addr,
					          int*            Charge1Addr,
					          int*            Charge2Addr,
					          double*         Mass,
					          double*         ResolutionOmega,
					          TLorentzVector* resultOmega,
					          int*            resultAccepted,
					          TLorentzVector* resultPi0,
					          TLorentzVector* resultPiMinus )
							: _Pi0_OAddr      ( Pi0_OAddr       ),
							  _Pi0_1Addr      ( Pi0_1Addr       ),
							  _Scattered0Addr ( Scattered0Addr  ),
							  _Scattered1Addr ( Scattered1Addr  ),
							  _Scattered2Addr ( Scattered2Addr  ),
							  _Charge0Addr    ( Charge0Addr     ),
							  _Charge1Addr    ( Charge1Addr     ),
							  _Charge2Addr    ( Charge2Addr     ),
							  _Mass           ( Mass            ),
							  _ResolutionOmega( ResolutionOmega ),
							  _resultOmega    ( resultOmega     ),
							  _resultAccepted ( resultAccepted  ),
							  _resultPi0      ( resultPi0       ),
							  _resultPiMinus  ( resultPiMinus   )
					{}

					virtual ~GetOmega() {}

					bool operator() ()
					{
						(*_resultAccepted) = 0;

						std::vector<TLorentzVector*> pi0s;
						pi0s.resize(2);
						pi0s.clear();
						pi0s.push_back(_Pi0_OAddr);
						pi0s.push_back(_Pi0_1Addr);

						std::vector<TLorentzVector*> chargedLV;
						chargedLV.resize(3);
						chargedLV.clear();
						chargedLV.push_back(_Scattered0Addr);
						chargedLV.push_back(_Scattered1Addr);
						chargedLV.push_back(_Scattered2Addr);

						std::vector<int*> charge;
						charge.resize(3);
						charge.clear();
						charge.push_back(_Charge0Addr);
						charge.push_back(_Charge1Addr);
						charge.push_back(_Charge2Addr);

						unsigned int numberCandidates = 0;

						// Loop over available pi0s
						for( unsigned int i = 0; i < pi0s.size(); i++ )
						{
							// Find suitable pi+/pi- pair
							for(unsigned int j = 0; j < chargedLV.size(); j++ )
							{
								for(unsigned int k = j + 1; k < chargedLV.size(); k++ )
								{
									// Check if charge is consistent with zero
									if( (*charge[j]) + (*charge[k]) == 0 )
									{
										// Check if mass fits
										const TLorentzVector temp = (*pi0s[i]) + (*chargedLV[j]) + (*chargedLV[k]);
										if( std::fabs( temp.Mag() - *_Mass ) < *_ResolutionOmega )
										{
											// Count candidates
											numberCandidates++;
											(*_resultOmega) = (*pi0s[i]) + (*chargedLV[j]) + (*chargedLV[k]);
											for( unsigned int l = 0; l < pi0s.size(); l++ )
											{
												if( i != l ) (*_resultPi0) = (*pi0s[l]);
											}
											for( unsigned int m = 0; m < chargedLV.size(); m++ )
											{
												if( j != m && k != m ) (*_resultPiMinus) = (*chargedLV[m]);
											}
										}
									}
								}
							}
						}

						if( numberCandidates ==  1 )
						{
							(*_resultAccepted) = 1;
						}

						return true;
					}

				private:
					TLorentzVector* _Pi0_OAddr;
					TLorentzVector* _Pi0_1Addr;
					TLorentzVector* _Scattered0Addr;
					TLorentzVector* _Scattered1Addr;
					TLorentzVector* _Scattered2Addr;
					int*            _Charge0Addr;
					int*            _Charge1Addr;
					int*            _Charge2Addr;
					double*         _Mass;
					double*         _ResolutionOmega;
					TLorentzVector* _resultOmega;
					int*            _resultAccepted;
					TLorentzVector* _resultPi0;
					TLorentzVector* _resultPiMinus;
				};

				class GetECALCorrectedEnergy : public Function
				{
				public:
					GetECALCorrectedEnergy( std::vector<double>                     *EnergyAddr,
					                        std::vector<double>                     *ClusterZ,
					                        double                                  *RangeECAL1,
					                        int                                     *RunNumberAddr,
					                        std::map<int, std::pair<double,double>>  correction,
					                        std::vector<double>                     *resultEnergy
					)
						: _EnergyAddr     ( EnergyAddr      ),
						  _ClusterZ       ( ClusterZ        ),
						  _RangeECAL1     ( RangeECAL1      ),
						  _RunNumberAddr  ( RunNumberAddr   ),
						  _correction     ( correction      ),
						  _resultEnergy   ( resultEnergy    ) {}

					virtual ~GetECALCorrectedEnergy() {}

					bool operator() ()
					{
						std::vector<double> correctedEnergies;
						correctedEnergies.resize( _EnergyAddr->size() );
						double correction;
						for( unsigned int i = 0; i < _EnergyAddr->size(); ++i )
						{
							if( _correction.count((*_RunNumberAddr)) == 0 )
							{
								correction = 1;
							}
							else
							{
								if( (*_ClusterZ)[i] < *_RangeECAL1 )
								{
									correction = _correction[(*_RunNumberAddr)].first;
								}
								else
								{
									correction = _correction[(*_RunNumberAddr)].second;
								}
							}
							correctedEnergies[i] = (*_EnergyAddr)[i] * correction;
						}
						(*_resultEnergy) = correctedEnergies;

						return true;
					}

				private:
					std::vector<double>                     *_EnergyAddr;
					std::vector<double>                     *_ClusterZ;
					double                                  *_RangeECAL1;
					int                                     *_RunNumberAddr;
					std::map<int, std::pair<double,double>>  _correction;
					std::vector<double>                     *_resultEnergy;
				};

				class GetECALCorrectedTiming : public Function
				{
				public:
					GetECALCorrectedTiming( std::vector<double>                       *Timing,
					                        std::vector<double>                       *Energy,
					                        std::vector<double>                       *ClusterZ,
					                        double                                    *RangeECAL1,
					                        std::map<std::string,std::vector<double>>  CalibrationValues,
					                        std::vector<double>                       *resultTiming
					)
							: _Timing           ( Timing            ),
							  _Energy           ( Energy            ),
							  _ClusterZ         ( ClusterZ          ),
							  _RangeECAL1       ( RangeECAL1        ),
							  _CalibrationValues( CalibrationValues ),
							  _resultTiming     ( resultTiming      ) {}

					virtual ~GetECALCorrectedTiming() {}

					bool operator() ()
					{
						std::vector<double> values;
						double correction;
						_resultTiming->resize( _Timing->size(), double() );
						for( unsigned int i = 0; i < _Timing->size(); ++i )
						{
							if( (*_ClusterZ)[i] < *_RangeECAL1 )
							{
								values     = _CalibrationValues["ECAL1"];
								correction = values[0] + values[1] /  (*_Energy)[i]                  - values[2] *  (*_Energy)[i]
								                       - values[3] / ((*_Energy)[i] * (*_Energy)[i]) + values[4] * ((*_Energy)[i] * (*_Energy)[i]);
							}
							else
							{
								values     = _CalibrationValues["ECAL2"];
								correction = values[0] + values[1] /  (*_Energy)[i]          + values[2] * ( *_Energy)[i]
								                       - values[3] / ((*_Energy)[i] * (*_Energy)[i])                 - values[4] * ((*_Energy)[i] * (*_Energy)[i])
								                       + values[5] / ((*_Energy)[i] * (*_Energy)[i] * (*_Energy)[i]) + values[6] * ((*_Energy)[i] * (*_Energy)[i] * (*_Energy)[i]);
							}
							(*_resultTiming)[i] = (*_Timing)[i] - correction;
						}
						return true;
					}

				private:
					std::vector<double>                       *_Timing;
					std::vector<double>                       *_Energy;
					std::vector<double>                       *_ClusterZ;
					double                                    *_RangeECAL1;
					std::map<std::string,std::vector<double>>  _CalibrationValues;
					std::vector<double>                       *_resultTiming;
				};

				class GetPhotonPairParticles : public Function
				{
					public:
						GetPhotonPairParticles( std::vector<TLorentzVector> *Photons,
						                        double                      *Mass,
						                        double                      *ECALResolution,
						                        double                      *ECAL1Resolution,
						                        double                      *ECAL2Resolution,
						                        std::vector<int>            *ECALIndex,
						                        std::vector<TLorentzVector> *resultParticles,
						                        int                         *resultHasParticles
						)
							: _Photons           ( Photons            ),
							  _Mass              ( Mass               ),
							  _ECALResolution    ( ECALResolution     ),
							  _ECAL1Resolution   ( ECAL1Resolution    ),
							  _ECAL2Resolution   ( ECAL2Resolution    ),
							  _ECALIndex         ( ECALIndex          ),
							  _resultParticles   ( resultParticles    ),
							  _resultHasParticles( resultHasParticles ) {}

					virtual ~GetPhotonPairParticles() {}

					bool operator()()
					{
						_resultParticles->reserve( _Photons->size() * (_Photons->size()) );
						_resultParticles->clear();
						if (_Photons->size() < 2)
						{
							return true;
						}
						for( unsigned int i = 0; i < _Photons->size(); i++)
						{
							for( unsigned int j = i+1; j < _Photons->size(); j++ )
							{
								const TLorentzVector candidate = (*_Photons)[i] + (*_Photons)[j];
								double massResolution;
								if( (*_ECALIndex)[i] == 1 && (*_ECALIndex)[j] == 1 )
								{
									massResolution = *_ECAL1Resolution;
								}
								else if( (*_ECALIndex)[i] == 2 && (*_ECALIndex)[j] == 2 )
								{
									massResolution = *_ECAL2Resolution;
								}
								else
								{
									massResolution = *_ECALResolution;
								}
								if (std::fabs(candidate.Mag() - *_Mass) < massResolution)
								{
									_resultParticles->push_back(candidate);
								}
							}
						}
						if( _resultParticles->size() > 1 )
						{
							(*_resultHasParticles) = 1;
						}
						else
						{
							(*_resultHasParticles) = 0;
						}
						return true;
					}

				private:
					std::vector<TLorentzVector> *_Photons;
					double                      *_Mass;
					double                      *_ECALResolution;
					double                      *_ECAL1Resolution;
					double                      *_ECAL2Resolution;
					std::vector<int>            *_ECALIndex;
					std::vector<TLorentzVector> *_resultParticles;
					int                         *_resultHasParticles;
				};

			}

		}

	}

}

#endif
