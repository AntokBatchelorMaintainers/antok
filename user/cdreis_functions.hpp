#ifndef ANTOK_USER_CDREIS_FUNCTIONS_HPP
#define ANTOK_USER_CDREIS_FUNCTIONS_HPP

#include <iostream>
#include <limits>
#include <vector>

#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TRotation.h"

#include "functions.hpp"
#include "neutral_fit.h"

namespace antok {

	namespace user {

		namespace cdreis {

			namespace functions {

				class GetRecoilLorentzVec : public Function {

				public:

					GetRecoilLorentzVec(const TLorentzVector& BeamLV,          // Lorentz vector of beam particle
					                    const TLorentzVector& XLV,             // Lorentz vector of X intermediate state
					                    const double&         RecoilMass,      // mass of recoil particle
					                    TLorentzVector&       ResultRecoilLV)  // Lorentz vector of recoil particle
						: _BeamLV        (BeamLV),
						  _XLV           (XLV),
						  _RecoilMass    (RecoilMass),
						  _ResultRecoilLV(ResultRecoilLV)
					{ }

					virtual ~GetRecoilLorentzVec() { }

					bool
					operator() ()
					{
						_ResultRecoilLV.SetVectM(_BeamLV.Vect() - _XLV.Vect(), _RecoilMass);
						return true;
					}

				private:

					const TLorentzVector& _BeamLV;
					const TLorentzVector& _XLV;
					const double          _RecoilMass;  // constant parameter, needs to be copied
					TLorentzVector&       _ResultRecoilLV;

				};


				class GetPhotonLorentzVecs : public Function
				{

				public:

					GetPhotonLorentzVecs(const std::vector<double>&   xPositions,         // x positions of ECAL photon clusters
					                     const std::vector<double>&   yPositions,         // y positions of ECAL photon clusters
					                     const std::vector<double>&   zPositions,         // z positions of ECAL photon clusters
					                     const std::vector<double>&   Energies,           // energies of ECAL photon clusters
					                     const std::vector<double>&   Times,              // times of ECAL photon clusters
					                     const double&                xPV,                // x positions of primary vertex
					                     const double&                yPV,                // y positions of primary vertex
					                     const double&                zPV,                // z positions of primary vertex
					                     const double&                zECAL1,             // z position that distinguishes between ECAL 1 and 2
					                     std::vector<TLorentzVector>& ResultLVs,          // photon Lorentz vectors
					                     std::vector<int>&            ResultECALIndices)  // indices of the ECAL that measured the photons
						: _xPositions       (xPositions),
						  _yPositions       (yPositions),
						  _zPositions       (zPositions),
						  _Energies         (Energies),
						  _Times            (Times),
						  _xPV              (xPV),
						  _yPV              (yPV),
						  _zPV              (zPV),
						  _zECAL1           (zECAL1),
						  _ResultLVs        (ResultLVs),
						  _ResultECALIndices(ResultECALIndices)
					{ }

					virtual ~GetPhotonLorentzVecs() { }

					bool
					operator() ()
					{
						const size_t nmbPhotons = _xPositions.size();
						if (   (_yPositions.size() != nmbPhotons)
						    or (_zPositions.size() != nmbPhotons)
						    or (_Energies.size()   != nmbPhotons)
						    or (_Times.size()      != nmbPhotons)) {
							return false;
						}

						_ResultLVs.resize(nmbPhotons);
						_ResultECALIndices.resize(nmbPhotons);
						for (size_t i = 0; i < nmbPhotons; ++i) {
							TVector3 mom(_xPositions[i] - _xPV, _yPositions[i] - _yPV, _zPositions[i] - _zPV);
							mom.SetMag(_Energies[i]);
							_ResultLVs[i] = TLorentzVector(mom, _Energies[i]);
							if (_zPositions[i] < _zECAL1) {
								_ResultECALIndices[i] = 1;
							} else {
								_ResultECALIndices[i] = 2;
							}
						}
						return true;
					}

				private:

					const std::vector<double>&   _xPositions;
					const std::vector<double>&   _yPositions;
					const std::vector<double>&   _zPositions;
					const std::vector<double>&   _Energies;
					const std::vector<double>&   _Times;
					const double&                _xPV;
					const double&                _yPV;
					const double&                _zPV;
					const double                 _zECAL1;  // constant parameter, needs to be copied
					std::vector<TLorentzVector>& _ResultLVs;
					std::vector<int>&            _ResultECALIndices;

				};


				class GetCleanedEcalClusters : public Function
				{

				public:

					GetCleanedEcalClusters(const std::vector<double>& xPositions,            // x positions of ECAL photon clusters
					                       const std::vector<double>& yPositions,            // y positions of ECAL photon clusters
					                       const std::vector<double>& zPositions,            // z positions of ECAL photon clusters
					                       const std::vector<double>& Energies,              // energies of ECAL photon clusters
					                       const std::vector<double>& Times,                 // times of ECAL photon clusters
					                       const double&              zECAL1,                // z position that distinguishes between ECAL 1 and 2
					                       const double&              ThresholdEnergyECAL1,  // energy threshold applied to ECAL1 clusters
					                       const double&              WindowTimingECAL1,     // window applied on ECAL1 cluster times
					                       const double&              ThresholdEnergyECAL2,  // energy threshold applied to ECAL2 clusters
					                       const double&              WindowTimingECAL2,     // window applied on ECAL1 cluster times
					                       std::vector<double>&       ResultXPositions,      // x positions of ECAL photon clusters
					                       std::vector<double>&       ResultYPositions,      // y positions of ECAL photon clusters
					                       std::vector<double>&       ResultZPositions,      // z positions of ECAL photon clusters
					                       std::vector<double>&       ResultEnergies,        // energies of ECAL photon clusters
					                       std::vector<double>&       ResultTimes,           // times of ECAL photon clusters
					                       std::vector<int>&          ResultECALIndices)     // indices of the ECAL that measured the photons
						: _xPositions          (xPositions),
						  _yPositions          (yPositions),
						  _zPositions          (zPositions),
						  _Energies            (Energies),
						  _Times               (Times),
						  _zECAL1              (zECAL1),
						  _ThresholdEnergyECAL1(ThresholdEnergyECAL1),
						  _WindowTimingECAL1   (WindowTimingECAL1),
						  _ThresholdEnergyECAL2(ThresholdEnergyECAL2),
						  _WindowTimingECAL2   (WindowTimingECAL2),
						  _ResultXPositions    (ResultXPositions),
						  _ResultYPositions    (ResultYPositions),
						  _ResultZPositions    (ResultZPositions),
						  _ResultEnergies      (ResultEnergies),
						  _ResultTimes         (ResultTimes),
						  _ResultECALIndices   (ResultECALIndices)
					{ }

					virtual ~GetCleanedEcalClusters() { }

					bool
					operator() ()
					{
						const size_t nmbPhotons = _xPositions.size();
						if (   (_yPositions.size() != nmbPhotons)
						    or (_zPositions.size() != nmbPhotons)
						    or (_Energies.size()   != nmbPhotons)
						    or (_Times.size()      != nmbPhotons)) {
							return false;
						}

						_ResultXPositions.reserve(nmbPhotons);
						_ResultXPositions.clear();
						_ResultYPositions.reserve(nmbPhotons);
						_ResultYPositions.clear();
						_ResultZPositions.reserve(nmbPhotons);
						_ResultZPositions.clear();
						_ResultEnergies.reserve(nmbPhotons);
						_ResultEnergies.clear();
						_ResultTimes.reserve(nmbPhotons);
						_ResultTimes.clear();
						_ResultECALIndices.reserve(nmbPhotons);
						_ResultECALIndices.clear();

						for (size_t i = 0; i < nmbPhotons; ++i) {
							if (_zPositions[i] < _zECAL1) {
								if ((_Energies[i] < _ThresholdEnergyECAL1) or (fabs(_Times[i]) > _WindowTimingECAL1)) {
									continue;
								}
								_ResultECALIndices.push_back(1);
							} else {
								if ((_Energies[i] < _ThresholdEnergyECAL2) or (fabs(_Times[i]) > _WindowTimingECAL2)) {
									continue;
								}
								_ResultECALIndices.push_back(2);
							}
							_ResultXPositions.push_back(_xPositions[i]);
							_ResultYPositions.push_back(_yPositions[i]);
							_ResultZPositions.push_back(_zPositions[i]);
							_ResultEnergies.push_back  (_Energies[i]);
							_ResultTimes.push_back     (_Times[i]);
						}
						return true;
					}

				private:

					const std::vector<double>& _xPositions;
					const std::vector<double>& _yPositions;
					const std::vector<double>& _zPositions;
					const std::vector<double>& _Energies;
					const std::vector<double>& _Times;
					const double               _zECAL1;                // constant parameter, needs to be copied
					const double               _ThresholdEnergyECAL1;  // constant parameter, needs to be copied
					const double               _WindowTimingECAL1;     // constant parameter, needs to be copied
					const double               _ThresholdEnergyECAL2;  // constant parameter, needs to be copied
					const double               _WindowTimingECAL2;     // constant parameter, needs to be copied
					std::vector<double>&       _ResultXPositions;
					std::vector<double>&       _ResultYPositions;
					std::vector<double>&       _ResultZPositions;
					std::vector<double>&       _ResultEnergies;
					std::vector<double>&       _ResultTimes;
					std::vector<int>&          _ResultECALIndices;

				};


				class GetVectorLorentzVectorAttributes : public Function
				{

				public:

					GetVectorLorentzVectorAttributes(const std::vector<TLorentzVector>& LVs,                // Lorentz vectors
					                                 std::vector<double>&               ResultXComponents,  // x components of Lorentz vectors
					                                 std::vector<double>&               ResultYComponents,  // y components of Lorentz vectors
					                                 std::vector<double>&               ResultZComponents,  // z components of Lorentz vectors
					                                 std::vector<double>&               ResultEnergies,     // energies of Lorentz vectors
					                                 std::vector<double>&               ResultThetas,       // polar angles of Lorentz vectors
					                                 std::vector<double>&               ResultPhis,         // azimuthal angles of Lorentz vectors
					                                 std::vector<double>&               ResultMags)         // magnitudes of Lorentz vectors
						: _LVs              (LVs),
						  _ResultXComponents(ResultXComponents),
						  _ResultYComponents(ResultYComponents),
						  _ResultZComponents(ResultZComponents),
						  _ResultEnergies   (ResultEnergies),
						  _ResultThetas     (ResultThetas),
						  _ResultPhis       (ResultPhis),
						  _ResultMags       (ResultMags)
					{ }

					virtual ~GetVectorLorentzVectorAttributes() { }

					bool operator() ()
					{
						const size_t nmbPhotons = _LVs.size();
						_ResultXComponents.resize(nmbPhotons);
						_ResultYComponents.resize(nmbPhotons);
						_ResultZComponents.resize(nmbPhotons);
						_ResultEnergies.resize   (nmbPhotons);
						_ResultThetas.resize     (nmbPhotons);
						_ResultPhis.resize       (nmbPhotons);
						_ResultMags.resize       (nmbPhotons);
						for (size_t i = 0; i < nmbPhotons; ++i) {
							_ResultXComponents[i] = _LVs[i].X();
							_ResultYComponents[i] = _LVs[i].Y();
							_ResultZComponents[i] = _LVs[i].Z();
							_ResultEnergies   [i] = _LVs[i].E();
							_ResultThetas     [i] = _LVs[i].Theta();
							_ResultPhis       [i] = _LVs[i].Phi();
							_ResultMags       [i] = _LVs[i].Mag();
						}
						return true;
					}

				private:

					const std::vector<TLorentzVector>& _LVs;
					std::vector<double>&               _ResultXComponents;
					std::vector<double>&               _ResultYComponents;
					std::vector<double>&               _ResultZComponents;
					std::vector<double>&               _ResultEnergies;
					std::vector<double>&               _ResultThetas;
					std::vector<double>&               _ResultPhis;
					std::vector<double>&               _ResultMags;

				};


				class GetPi0Pair : public Function
				{

				public:

					GetPi0Pair(const std::vector<TLorentzVector>& PhotonLVs,            // Lorentz vectors of photons
					           const std::vector<int>&            ECALIndices,          // indices of the ECAL that measured the photons
					           const double&                      Pi0Mass,              // nominal pi^0 mass
					           const double&                      MassWindowECALMixed,  // m(gamma gamma) cut applied around Pi0Mass when the photon pair is in ECAL1 and 2
					           const double&                      MassWindowECAL1,      // m(gamma gamma) cut applied around Pi0Mass when both photons are in ECAL1
					           const double&                      MassWindowECAL2,      // m(gamma gamma) cut applied around Pi0Mass when both photons are in ECAL2
					           std::vector<TLorentzVector>&       ResultPi0PairLVs,     // Lorentz vectors of the two pi^0 in the first found pair
					           TLorentzVector&                    ResultPi0LV_0,        // Lorentz vectors of the 1st pi^0 in the first found pair
					           TLorentzVector&                    ResultPi0LV_1,        // Lorentz vectors of the 2nd pi^0 in the first found pair
					           int&                               ResultGoodPi0Pair,    // 1 if exactly one pi^0 pair was found; else 0
					           std::vector<int>&                  ResultECALIndices)    // indices of the ECAL that measured the selected photons
						: _PhotonLVs          (PhotonLVs),
						  _ECALIndices        (ECALIndices),
						  _Pi0Mass            (Pi0Mass),
						  _MassWindowECALMixed(MassWindowECALMixed),
						  _MassWindowECAL1    (MassWindowECAL1),
						  _MassWindowECAL2    (MassWindowECAL2),
						  _ResultPi0PairLVs   (ResultPi0PairLVs),
						  _ResultPi0LV_0      (ResultPi0LV_0),
						  _ResultPi0LV_1      (ResultPi0LV_1),
						  _ResultGoodPi0Pair  (ResultGoodPi0Pair),
						  _ResultECALIndices  (ResultECALIndices)
					{ }

					virtual ~GetPi0Pair() { }

					bool
					operator() ()
					{
						const size_t nmbPhotons = _PhotonLVs.size();
						if (_ECALIndices.size() != nmbPhotons) {
							return false;
						}

						_ResultPi0PairLVs.reserve(2);
						_ResultPi0PairLVs.clear();
						_ResultPi0LV_0.SetPxPyPzE(0, 0, 0, 0);
						_ResultPi0LV_1.SetPxPyPzE(0, 0, 0, 0);
						_ResultGoodPi0Pair = 0;
						_ResultECALIndices = {-1, -1, -1, -1};
						if (nmbPhotons < 2) {
							return true;
						}

						// searches for pairs of pi^0 candidates
						// search is aborted if more than one pair is found
						// always the first found pair is returned
						size_t nmbCandidatePairs = 0;
						for (size_t i = 0; i < nmbPhotons; ++i) {
							if (nmbCandidatePairs > 1) {
								break;
							}
							for (size_t j = i + 1; j < nmbPhotons; ++j) {
								if (nmbCandidatePairs > 1) {
									break;
								}
								const TLorentzVector pi0Candidate0 = _PhotonLVs[i] + _PhotonLVs[j];
								if (std::fabs(pi0Candidate0.M() - _Pi0Mass) > getECALMassWindow(_ECALIndices[i] == 1, _ECALIndices[j])) {
									continue;
								}
								for (size_t m = i + 1; m < nmbPhotons; ++m) {
									if (nmbCandidatePairs > 1) {
										break;
									}
									for (size_t n = m + 1; n < nmbPhotons; ++n) {
										if (nmbCandidatePairs > 1) {
											break;
										}
										if ((m == j) or (n == j)) {
											continue;
										}
										const TLorentzVector pi0Candidate1 = _PhotonLVs[m] + _PhotonLVs[n];
										if (std::fabs(pi0Candidate1.M() - _Pi0Mass) > getECALMassWindow(_ECALIndices[m] == 1, _ECALIndices[n])) {
											continue;
										}
										if (nmbCandidatePairs == 0) {
											_ResultPi0PairLVs.push_back(pi0Candidate0);
											_ResultPi0PairLVs.push_back(pi0Candidate1);
											_ResultPi0LV_0 = pi0Candidate0;
											_ResultPi0LV_1 = pi0Candidate1;
											_ResultECALIndices = {(int)i, (int)j, (int)m, (int)n};
										}
										nmbCandidatePairs++;
									}
								}
							}
						}
						if (nmbCandidatePairs == 1) {
							_ResultGoodPi0Pair = 1;
						}

						// std::cout << "!!!HERE nmbPhotons = " << nmbPhotons << ", nmbCandidatePairs = " << nmbCandidatePairs << std::endl;
						return true;
					}

				private:

					double getECALMassWindow(const int ECALIndex1,
					                         const int ECALIndex2) const
					{
						if        ((ECALIndex1 == 1) and (ECALIndex2 == 1)) {
							return _MassWindowECAL1;
						} else if ((ECALIndex1 == 2) and (ECALIndex2 == 2)) {
							return _MassWindowECAL2;
						} else if (ECALIndex1 != ECALIndex2) {
							return _MassWindowECALMixed;
						}
						//TODO maybe throwing an exception would be more meaningful here?
						return std::numeric_limits<double>::quiet_NaN();
					}

					const std::vector<TLorentzVector>& _PhotonLVs;
					const std::vector<int>&            _ECALIndices;
					const double                       _Pi0Mass;              // constant parameter, needs to be copied
					const double                       _MassWindowECALMixed;  // constant parameter, needs to be copied
					const double                       _MassWindowECAL1;      // constant parameter, needs to be copied
					const double                       _MassWindowECAL2;      // constant parameter, needs to be copied
					std::vector<TLorentzVector>&       _ResultPi0PairLVs;
					TLorentzVector&                    _ResultPi0LV_0;
					TLorentzVector&                    _ResultPi0LV_1;
					int&                               _ResultGoodPi0Pair;
					std::vector<int>&                  _ResultECALIndices;

				};


				class GetOmega : public Function
				{

				public:

					GetOmega(const TLorentzVector& Pi0LV_O,                 // Lorentz vector of 1st pi^0
					         const TLorentzVector& Pi0LV_1,                 // Lorentz vector of 2nd pi^0
					         const TLorentzVector& ChargedPartLV_0,         // Lorentz vector of 1st charged particle
					         const TLorentzVector& ChargedPartLV_1,         // Lorentz vector of 2nd charged particle
					         const TLorentzVector& ChargedPartLV_2,         // Lorentz vector of 3rd charged particle
					         const int&            Charge_0,                // charge of 1st charged particle
					         const int&            Charge_1,                // charge of 2nd charged particle
					         const int&            Charge_2,                // charge of 3rd charged particle
					         const double&         OmegaMass,               // nominal omega mass
					         const double&         MassWindowOmega,         // cut around OmegaMass applied on m(pi^- pi^0 pi^+)
					         TLorentzVector&       ResultOmegaLV,           // Lorentz vector of last found omega candidate
					         int&                  ResultAccepted,          // 1 if there is exactly one omega candidate, 0 otherwise
					         TLorentzVector&       ResultNotUsedPi0LV,      // Lorentz vector of the pi^0 that is not part of the omega
					         TLorentzVector&       ResultNotUsedPiMinusLV)  // Lorentz vector of the pi^- that is not part of the omega
						: _Pi0LV_O               (Pi0LV_O),
						  _Pi0LV_1               (Pi0LV_1),
						  _ChargedPartLV_0       (ChargedPartLV_0),
						  _ChargedPartLV_1       (ChargedPartLV_1),
						  _ChargedPartLV_2       (ChargedPartLV_2),
						  _Charge_0              (Charge_0),
						  _Charge_1              (Charge_1),
						  _Charge_2              (Charge_2),
						  _OmegaMass             (OmegaMass),
						  _MassWindowOmega       (MassWindowOmega),
						  _ResultOmegaLV         (ResultOmegaLV),
						  _ResultAccepted        (ResultAccepted),
						  _ResultNotUsedPi0LV    (ResultNotUsedPi0LV),
						  _ResultNotUsedPiMinusLV(ResultNotUsedPiMinusLV)
					{ }

					virtual ~GetOmega() { }

					bool
					operator() ()
					{
						// searches for omega(782) -> pi^- pi^0 pi^+ candidates in pi^- pi^- pi^0 pi^0 pi^+ final state
						// if several omega candidates are found, the last found one is returned
						const std::vector<const TLorentzVector*> pi0s       = {&_Pi0LV_O, &_Pi0LV_1};
						const std::vector<const TLorentzVector*> chargedLVs = {&_ChargedPartLV_0, &_ChargedPartLV_1, &_ChargedPartLV_2};
						const std::vector<const int*>            charges    = {&_Charge_0,        &_Charge_1,        &_Charge_2};
						size_t numberCandidates = 0;
						_ResultAccepted = 0;
						// Loop over available pi^0s
						for (size_t i = 0; i < pi0s.size(); ++i) {
							// Loop over charged particles
							for (size_t j = 0; j < chargedLVs.size(); ++j) {
								for (size_t k = j + 1; k < chargedLVs.size(); ++k) {
									// Find pi+ pi- pair
									if (*charges[j] + *charges[k] == 0) {
										// Check if mass fits omega(782) nominal mass
										const TLorentzVector candidate = *pi0s[i] + *chargedLVs[j] + *chargedLVs[k];
										if (std::fabs(candidate.M() - _OmegaMass) < _MassWindowOmega) {
											_ResultOmegaLV = candidate;
											numberCandidates++;  // Count omega candidates
											// find pi^0 that is not part of the omega candidate
											for (size_t l = 0; l < pi0s.size(); l++) {
												if (i != l) {
													_ResultNotUsedPi0LV = *pi0s[l];
												}
											}
											// find pi^- that is not part of the omega candidate
											for (size_t m = 0; m < chargedLVs.size(); ++m) {
												if ((m != j) and (m != k)) {
													_ResultNotUsedPiMinusLV = *chargedLVs[m];
												}
											}
										}
									}
								}
							}
						}
						if (numberCandidates ==  1) {
							_ResultAccepted = 1;
						}

						return true;
					}

				private:

					const TLorentzVector& _Pi0LV_O;
					const TLorentzVector& _Pi0LV_1;
					const TLorentzVector& _ChargedPartLV_0;
					const TLorentzVector& _ChargedPartLV_1;
					const TLorentzVector& _ChargedPartLV_2;
					const int&            _Charge_0;
					const int&            _Charge_1;
					const int&            _Charge_2;
					const double          _OmegaMass;        // constant parameter, needs to be copied
					const double          _MassWindowOmega;  // constant parameter, needs to be copied
					TLorentzVector&       _ResultOmegaLV;
					int&                  _ResultAccepted;
					TLorentzVector&       _ResultNotUsedPi0LV;
					TLorentzVector&       _ResultNotUsedPiMinusLV;

				};


				class GetECALCorrectedEnergy : public Function
				{

				public:

					GetECALCorrectedEnergy(const std::vector<double>&                      Energies,                 // energies of ECAL clusters to be corrected
					                       const std::vector<double>&                      zPositions,               // z positions of ECAL clusters
					                       const double&                                   zECAL1,                   // z position that distinguishes between ECAL 1 and 2
					                       const int&                                      RunNumber,                // run number of the event
					                       const std::map<int, std::pair<double, double>>& Corrections,              // energy-correction factors for ECAL1 and 2 by run number
					                       std::vector<double>&                            ResultCorrectedEnergies)  // corrected energies of ECAL clusters
						: _Energies               (Energies),
						  _zPositions             (zPositions),
						  _zECAL1                 (zECAL1),
						  _RunNumber              (RunNumber),
						  _Corrections            (Corrections),
						  _ResultCorrectedEnergies(ResultCorrectedEnergies)
					{ }

					virtual ~GetECALCorrectedEnergy() { }

					bool operator() ()
					{
						const size_t nmbClusters = _Energies.size();
						if (_zPositions.size() != nmbClusters) {
							return false;
						}

						_ResultCorrectedEnergies.resize(nmbClusters);
						for (size_t i = 0; i < nmbClusters; ++i) {
							double correction = 1;
							auto it = _Corrections.find(_RunNumber);
							if (it != _Corrections.end()) {
								if (_zPositions[i] < _zECAL1) {
									correction = it->second.first;
								} else {
									correction = it->second.second;
								}
							}
							_ResultCorrectedEnergies[i] = correction * _Energies[i];
						}
						return true;
					}

				private:

					const std::vector<double>&                     _Energies;
					const std::vector<double>&                     _zPositions;
					const double                                   _zECAL1;       // constant parameter, needs to be copied
					const int&                                     _RunNumber;
					const std::map<int, std::pair<double, double>> _Corrections;  // constant parameter, needs to be copied
					std::vector<double>&                           _ResultCorrectedEnergies;

				};


				class GetECALCorrectedTiming : public Function
				{

				public:

					GetECALCorrectedTiming(const std::vector<double>&                        Times,                 // times of ECAL clusters to be corrected
					                       const std::vector<double>&                        Energies,              // energies of ECAL clusters
					                       const std::vector<double>&                        zPositions,            // z positions of ECAL clusters
					                       const double&                                     zECAL1,                // z position that distinguishes between ECAL 1 and 2
					                       const std::map<std::string, std::vector<double>>& CalibrationCoeffs,     // calibration coefficients used to correct times
					                       std::vector<double>&                              ResultCorrectedTimes)  // corrected times of ECAL clusters
						: _Times               (Times),
						  _Energies            (Energies),
						  _zPositions          (zPositions),
						  _zECAL1              (zECAL1),
						  _CalibrationCoeffs   (CalibrationCoeffs),
						  _ResultCorrectedTimes(ResultCorrectedTimes)
					{ }

					virtual ~GetECALCorrectedTiming() { }

					bool
					operator() ()
					{
						const size_t nmbClusters = _Times.size();
						if ((_Energies.size() != nmbClusters) and (_zPositions.size() != nmbClusters)) {
							return false;
						}

						_ResultCorrectedTimes.resize(nmbClusters);
						for (size_t i = 0; i < nmbClusters; ++i) {
							double       correction = 0;
							const double energy     = _Energies[i];
							const double energy2    = energy * energy;
							if (_zPositions[i] < _zECAL1) {
								const std::vector<double>& coefficients = _CalibrationCoeffs.at("ECAL1");
								correction = coefficients[0] + coefficients[1] / energy  - coefficients[2] * energy
								                             - coefficients[3] / energy2 + coefficients[4] * energy2;
							} else {
							  const double               energy3      = energy * energy2;
								const std::vector<double>& coefficients = _CalibrationCoeffs.at("ECAL2");
								correction = coefficients[0] + coefficients[1] / energy  + coefficients[2] * energy
								                             - coefficients[3] / energy2 - coefficients[4] * energy2
								                             + coefficients[5] / energy3 + coefficients[6] * energy3;
							}
							_ResultCorrectedTimes[i] = _Times[i] - correction;
						}
						return true;
					}

				private:

					const std::vector<double>&                        _Times;
					const std::vector<double>&                        _Energies;
					const std::vector<double>&                        _zPositions;
					const double                                      _zECAL1;             // constant parameter, needs to be copied
					const std::map<std::string, std::vector<double>>  _CalibrationCoeffs;  // constant parameter, needs to be copied
					std::vector<double>&                              _ResultCorrectedTimes;

				};


				class GetPhotonPairParticles : public Function
				{

				public:

					//TODO is ResultHasParticles really needed? size of ResultParticleLVs contains the same info
					GetPhotonPairParticles(const std::vector<TLorentzVector>& PhotonLVs,            // Lorentz vectors of photons
					                       const double&                      NominalMass,          // nominal mass of photon pair
					                       const double&                      MassWindowECALMixed,  // m(gamma gamma) cut applied around NominalMass when the photon pair is in ECAL1 and 2
					                       const double&                      MassWindowECAL1,      // m(gamma gamma) cut applied around NominalMass when both photons are in ECAL1
					                       const double&                      MassWindowECAL2,      // m(gamma gamma) cut applied around NominalMass when both photons are in ECAL2
					                       const std::vector<int>&            ECALIndices,          // indices of the ECAL that measured the photons
					                       std::vector<TLorentzVector>&       ResultParticleLVs,    // Lorentz vectors of all particles reconstructed from photon pairs
					                       int&                               ResultHasParticles)   // 1 if at least one photon pair was found; else 0
						: _PhotonLVs          (PhotonLVs),
						  _NominalMass        (NominalMass),
						  _MassWindowECALMixed(MassWindowECALMixed),
						  _MassWindowECAL1    (MassWindowECAL1),
						  _MassWindowECAL2    (MassWindowECAL2),
						  _ECALIndices        (ECALIndices),
						  _ResultParticles    (ResultParticleLVs),
						  _ResultHasParticles (ResultHasParticles)
					{ }

					virtual ~GetPhotonPairParticles() { }

					bool
					operator() ()
					{
						const size_t nmbPhotons = _PhotonLVs.size();
						if (_ECALIndices.size() != nmbPhotons) {
							return false;
						}

						_ResultParticles.reserve(nmbPhotons * nmbPhotons);
						_ResultParticles.clear();
						_ResultHasParticles = 0;
						if (nmbPhotons < 2) {
							return true;
						}
						for (size_t i = 0; i < nmbPhotons; ++i) {
							for (size_t j = i + 1; j < nmbPhotons; ++j) {
								const TLorentzVector candidate = _PhotonLVs[i] + _PhotonLVs[j];
								if (std::fabs(candidate.M() - _NominalMass) < getECALResolution(_ECALIndices[i] == 1, _ECALIndices[j])) {
									_ResultParticles.push_back(candidate);
								}
							}
						}
						if (_ResultParticles.size() > 1) {
							_ResultHasParticles = 1;
						}
						return true;
					}

				private:

					double getECALResolution(const int ECALIndex1,
					                         const int ECALIndex2) const
					{
						if        ((ECALIndex1 == 1) and (ECALIndex2 == 1)) {
							return _MassWindowECAL1;
						} else if ((ECALIndex1 == 2) and (ECALIndex2 == 2)) {
							return _MassWindowECAL2;
						} else if (ECALIndex1 != ECALIndex2) {
							return _MassWindowECALMixed;
						}
						return std::numeric_limits<double>::quiet_NaN();
					}

					const std::vector<TLorentzVector>& _PhotonLVs;
					const double                       _NominalMass;          // constant parameter, needs to be copied
					const double                       _MassWindowECALMixed;  // constant parameter, needs to be copied
					const double                       _MassWindowECAL1;      // constant parameter, needs to be copied
					const double                       _MassWindowECAL2;      // constant parameter, needs to be copied
					const std::vector<int>&            _ECALIndices;
					std::vector<TLorentzVector>&       _ResultParticles;
					int&                               _ResultHasParticles;

				};


				class GetKinematicFittingMass : public Function
				{

				public:

					GetKinematicFittingMass(const TVector3&              VertexPosition,
					                        const std::vector<TVector3>& ClusterPositions,
					                        const std::vector<TVector3>& ClusterPositionErrors,
					                        const std::vector<double>&   ClusterEnergies,
					                        const std::vector<double>&   ClusterEnergyErrors,
					                        const std::vector<int>&      ClusterIndices,
					                        const double&                Mass,
					                        const double&                Window,
					                        const int&                   EnergyErrorType,
					                        std::vector<TLorentzVector>& ResultLorentzVectors,
					                        std::vector<double>&         ResultChi2s,
					                        std::vector<double>&         ResultPullsX0,
					                        std::vector<double>&         ResultPullsY0,
					                        std::vector<double>&         ResultPullsE0,
					                        std::vector<double>&         ResultPullsX1,
					                        std::vector<double>&         ResultPullsY1,
					                        std::vector<double>&         ResultPullsE1,
					                        std::vector<double>&         ResultCLs,
					                        int&                         ResultSuccess)
						: _VertexPosition       (VertexPosition),
						  _ClusterPositions     (ClusterPositions),
						  _ClusterPositionErrors(ClusterPositionErrors),
						  _ClusterEnergies      (ClusterEnergies),
						  _ClusterEnergyErrors  (ClusterEnergyErrors),
						  _ClusterIndices       (ClusterIndices),
						  _Mass                 (Mass),
						  _Window               (Window),
						  _EnergyErrorType      (EnergyErrorType),
						  _ResultLorentzVectors (ResultLorentzVectors),
						  _ResultChi2s          (ResultChi2s),
						  _ResultPullsX0        (ResultPullsX0),
						  _ResultPullsY0        (ResultPullsY0),
						  _ResultPullsE0        (ResultPullsE0),
						  _ResultPullsX1        (ResultPullsX1),
						  _ResultPullsY1        (ResultPullsY1),
						  _ResultPullsE1        (ResultPullsE1),
						  _ResultCLs            (ResultCLs),
						  _ResultSuccess        (ResultSuccess)
					{ }

					virtual ~GetKinematicFittingMass() { }

					bool
					operator() ()
					{
						if (   (_ClusterIndices.size()        != 4)
						    or (_ClusterPositions.size()      != _ClusterEnergies.size())
						    or (_ClusterPositionErrors.size() != _ClusterEnergyErrors.size())) {
							return false;
						}

						_ResultLorentzVectors.reserve(2);
						_ResultLorentzVectors.clear();
						_ResultChi2s.reserve(2);
						_ResultChi2s.clear();
						_ResultPullsX0.reserve(2);
						_ResultPullsX0.clear();
						_ResultPullsY0.reserve(2);
						_ResultPullsY0.clear();
						_ResultPullsE0.reserve(2);
						_ResultPullsE0.clear();
						_ResultPullsX1.reserve(2);
						_ResultPullsX1.clear();
						_ResultPullsY1.reserve(2);
						_ResultPullsY1.clear();
						_ResultPullsE1.reserve(2);
						_ResultPullsE1.clear();
						_ResultCLs.reserve(2);
						_ResultCLs.clear();
						_ResultSuccess = 0;

						if (   (_ClusterIndices[0] == -1)
						    or (_ClusterIndices[1] == -1)
						    or (_ClusterIndices[2] == -1)
						    or (_ClusterIndices[3] == -1)) {
							return true;
						}
						const double maxPosError = 1e3;
						for (size_t i = 0; i < 4; ++i) {
							if (   (_ClusterPositionErrors[_ClusterIndices[i]].X() > maxPosError)
							    or (_ClusterPositionErrors[_ClusterIndices[i]].Y() > maxPosError)
							    or (_ClusterPositionErrors[_ClusterIndices[i]].Z() > maxPosError)) {
								return true;
							}
						}

						antok::NeutralFit neutralFit0(_VertexPosition,
						                              _ClusterPositions     [_ClusterIndices[0]],
						                              _ClusterPositions     [_ClusterIndices[1]],
						                              _ClusterPositionErrors[_ClusterIndices[0]],
						                              _ClusterPositionErrors[_ClusterIndices[1]],
						                              _ClusterEnergies      [_ClusterIndices[0]],
						                              _ClusterEnergies      [_ClusterIndices[1]],
						                              _ClusterEnergyErrors  [_ClusterIndices[0]],
						                              _ClusterEnergyErrors  [_ClusterIndices[1]],
						                              _Mass,
						                              _Window,
						                              _EnergyErrorType);
						const bool success0 = neutralFit0.doFit();
						if (success0) {
							_ResultLorentzVectors.push_back(neutralFit0.getLVSum());
							_ResultPullsX0.push_back       (neutralFit0.getPulls()[0]);
							_ResultPullsY0.push_back       (neutralFit0.getPulls()[1]);
							_ResultPullsE0.push_back       (neutralFit0.getPulls()[2]);
							_ResultPullsX1.push_back       (neutralFit0.getPulls()[3]);
							_ResultPullsY1.push_back       (neutralFit0.getPulls()[4]);
							_ResultPullsE1.push_back       (neutralFit0.getPulls()[5]);
							_ResultChi2s.push_back         (neutralFit0.getChi2());
							_ResultCLs.push_back           (neutralFit0.getCL());
						}

						antok::NeutralFit neutralFit1(_VertexPosition,
						                              _ClusterPositions     [_ClusterIndices[2]],
						                              _ClusterPositions     [_ClusterIndices[3]],
						                              _ClusterPositionErrors[_ClusterIndices[2]],
						                              _ClusterPositionErrors[_ClusterIndices[3]],
						                              _ClusterEnergies      [_ClusterIndices[2]],
						                              _ClusterEnergies      [_ClusterIndices[3]],
						                              _ClusterEnergyErrors  [_ClusterIndices[2]],
						                              _ClusterEnergyErrors  [_ClusterIndices[3]],
						                              _Mass,
						                              _Window,
						                              _EnergyErrorType);
						const bool success1 = neutralFit1.doFit();
						if (success1) {
							_ResultLorentzVectors.push_back(neutralFit1.getLVSum());
							_ResultPullsX0.push_back       (neutralFit1.getPulls()[0]);
							_ResultPullsY0.push_back       (neutralFit1.getPulls()[1]);
							_ResultPullsE0.push_back       (neutralFit1.getPulls()[2]);
							_ResultPullsX1.push_back       (neutralFit1.getPulls()[3]);
							_ResultPullsY1.push_back       (neutralFit1.getPulls()[4]);
							_ResultPullsE1.push_back       (neutralFit1.getPulls()[5]);
							_ResultChi2s.push_back         (neutralFit1.getChi2());
							_ResultCLs.push_back           (neutralFit1.getCL());
						}

						if (success0 and success1) {
							_ResultSuccess = 1;
						}
						return true;
					}

				private:

					const TVector3&              _VertexPosition;
					const std::vector<TVector3>& _ClusterPositions;
					const std::vector<TVector3>& _ClusterPositionErrors;
					const std::vector<double>&   _ClusterEnergies;
					const std::vector<double>&   _ClusterEnergyErrors;
					const std::vector<int>&      _ClusterIndices;
					const double                 _Mass;             // constant parameter, needs to be copied
					const double                 _Window;           // constant parameter, needs to be copied
					const int                    _EnergyErrorType;  // constant parameter, needs to be copied
					std::vector<TLorentzVector>& _ResultLorentzVectors;
					std::vector<double>&         _ResultChi2s;
					std::vector<double>&         _ResultPullsX0;
					std::vector<double>&         _ResultPullsY0;
					std::vector<double>&         _ResultPullsE0;
					std::vector<double>&         _ResultPullsX1;
					std::vector<double>&         _ResultPullsY1;
					std::vector<double>&         _ResultPullsE1;
					std::vector<double>&         _ResultCLs;
					int&                         _ResultSuccess;

				};


				class GetThreePionCombinationMass : public Function
				{

				public:

					GetThreePionCombinationMass(const TLorentzVector& Pi0LV_0,          // Lorentz vector of 1st pi^0
					                            const TLorentzVector& Pi0LV_1,          // Lorentz vector of 2nd pi^0
					                            const TLorentzVector& ChargedPartLV_0,  // Lorentz vector of 1st charged particle
					                            const TLorentzVector& ChargedPartLV_1,  // Lorentz vector of 2nd charged particle
					                            const TLorentzVector& ChargedPartLV_2,  // Lorentz vector of 3rd charged particle
					                            const int&            Charge_0,         // charge of 1st charged particle
					                            const int&            Charge_1,         // charge of 2nd charged particle
					                            const int&            Charge_2,         // charge of 3rd charged particle
					                            const int&            UseMassSquared,   // switch between mass squared (= 1) and mass (all other values)
					                            std::vector<double>&  Result)           // mass or mass squared (see above)
						: _Pi0LV_0        (Pi0LV_0),
						  _Pi0LV_1        (Pi0LV_1),
						  _ChargedPartLV_0(ChargedPartLV_0),
						  _ChargedPartLV_1(ChargedPartLV_1),
						  _ChargedPartLV_2(ChargedPartLV_2),
						  _Charge_0       (Charge_0),
						  _Charge_1       (Charge_1),
						  _Charge_2       (Charge_2),
						  _UseMassSquared (UseMassSquared),
						  _Result         (Result)
					{ }

					virtual ~GetThreePionCombinationMass() { }

					bool
					operator() ()
					{
						_Result.reserve(4);
						_Result.clear();
						const std::vector<const TLorentzVector*> Pi0LVs         = {&_Pi0LV_0, &_Pi0LV_1};
						const std::vector<const TLorentzVector*> ChargedPartLVs = {&_ChargedPartLV_0, &_ChargedPartLV_1, &_ChargedPartLV_2};
						const std::vector<int>                   Charges        = {_Charge_0, _Charge_1, _Charge_2};

						for (size_t i = 0; i < Pi0LVs.size(); ++i) {
							const TLorentzVector* Pi0LV = Pi0LVs[i];
							for (size_t j = 0; j < ChargedPartLVs.size(); ++j) {
								const int             chargeFirst  = Charges[j];
								const TLorentzVector* chargedFirst = ChargedPartLVs[j];
								for (size_t k = j+1; k < ChargedPartLVs.size(); ++k) {
									const int             chargeSecond  = Charges[k];
									const TLorentzVector* chargedSecond = ChargedPartLVs[k];
									if (chargeFirst == chargeSecond) {
										continue;
									} else {
										const TLorentzVector sum = *Pi0LV + *chargedFirst + *chargedSecond;
										if (_UseMassSquared == 1) {
											_Result.push_back(sum.M2());
										} else {
											_Result.push_back(sum.M());
										}
									}
								}
							}
						}

						return true;
					}

				private:

					const TLorentzVector& _Pi0LV_0;
					const TLorentzVector& _Pi0LV_1;
					const TLorentzVector& _ChargedPartLV_0;
					const TLorentzVector& _ChargedPartLV_1;
					const TLorentzVector& _ChargedPartLV_2;
					const int&            _Charge_0;
					const int&            _Charge_1;
					const int&            _Charge_2;
					const int             _UseMassSquared;  // constant parameter, needs to be copied
					std::vector<double>&  _Result;

				};


			}  // functions namespace

		}  // cdreis namespace

	}  // user namespace

}  // antok namespace

#endif  // ANTOK_USER_CDREIS_FUNCTIONS_HPP
