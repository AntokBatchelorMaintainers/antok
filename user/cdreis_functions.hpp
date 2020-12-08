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

				class GetVector3VectorAttributes : public Function
				{

				public:

					GetVector3VectorAttributes(const std::vector<TVector3>& Vector3s,           // TVector3 vectors
					                                 std::vector<double>&   ResultXComponents,  // x components of TVector3 vectors
					                                 std::vector<double>&   ResultYComponents,  // y components of TVector3 vectors
					                                 std::vector<double>&   ResultZComponents)  // z components of TVector3 vectors
						: _Vector3s         (Vector3s),
						  _ResultXComponents(ResultXComponents),
						  _ResultYComponents(ResultYComponents),
						  _ResultZComponents(ResultZComponents)
					{ }

					virtual ~GetVector3VectorAttributes() { }

					bool
					operator() ()
					{
						const size_t vecSize = _Vector3s.size();
						_ResultXComponents.resize(vecSize);
						_ResultYComponents.resize(vecSize);
						_ResultZComponents.resize(vecSize);

						// loop over TVector3 vector and split into vector for each component
						for (size_t i = 0; i < vecSize; ++i) {
							_ResultXComponents[i] = _Vector3s[i].X();
							_ResultYComponents[i] = _Vector3s[i].Y();
							_ResultZComponents[i] = _Vector3s[i].Z();
						}
						return true;
					}

				private:

					const std::vector<TVector3>& _Vector3s;
					std::vector<double>&         _ResultXComponents;
					std::vector<double>&         _ResultYComponents;
					std::vector<double>&         _ResultZComponents;

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
						const size_t vecSize = _LVs.size();
						_ResultXComponents.resize(vecSize);
						_ResultYComponents.resize(vecSize);
						_ResultZComponents.resize(vecSize);
						_ResultEnergies.resize   (vecSize);
						_ResultThetas.resize     (vecSize);
						_ResultPhis.resize       (vecSize);
						_ResultMags.resize       (vecSize);
						for (size_t i = 0; i < vecSize; ++i) {
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


				class getNominalMassDifferences : public Function
				{

				public:

					getNominalMassDifferences(const std::vector<TLorentzVector>& VectorLV,
					                          const double                       NominalMass,
					                          std::vector<double>&               MassDifferences)
						: _VectorLV       (VectorLV),
						  _NominalMass    (NominalMass),
						  _MassDifferences(MassDifferences)
					{ }

					virtual ~getNominalMassDifferences() { }

					bool
					operator() ()
					{
						const size_t sizeVec = _VectorLV.size();
						_MassDifferences.resize(sizeVec);
						for (size_t i = 0; i < sizeVec; ++i) {
							_MassDifferences[i] = _VectorLV[i].M() - _NominalMass;
						}
						return true;
					}

				private:

					const std::vector<TLorentzVector>& _VectorLV;
					const double                       _NominalMass;  // constant parameter, needs to be copied
					std::vector<double>&               _MassDifferences;

				};


				class GetRecoilLorentzVec : public Function
				{

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


				class GetECALCorrectedEnergy : public Function
				{

				public:

					GetECALCorrectedEnergy(const std::vector<double>&                      Energies,                 // energies of ECAL clusters to be corrected
					                       const std::vector<int>&                         ECALClusterIndices,       // Cluster Indices of ECAL hits
					                       const int&                                      RunNumber,                // run number of the event
					                       const std::map<int, std::pair<double, double>>& Corrections,              // energy-correction factors for ECAL1 and 2 by run number
					                       std::vector<double>&                            ResultCorrectedEnergies)  // corrected energies of ECAL clusters
						: _Energies               (Energies),
						  _ECALClusterIndices     (ECALClusterIndices),
						  _RunNumber              (RunNumber),
						  _Corrections            (Corrections),
						  _ResultCorrectedEnergies(ResultCorrectedEnergies)
					{ }

					virtual ~GetECALCorrectedEnergy() { }

					bool operator() ()
					{
						const size_t nmbClusters = _Energies.size();
						if (_ECALClusterIndices.size() != nmbClusters) {
							std::cerr << "Input vectors do not have the same size." << std::endl;
							return false;
						}

						_ResultCorrectedEnergies.resize(nmbClusters);
						for (size_t i = 0; i < nmbClusters; ++i) {
							double correction = 1;  // default: no correction
							const auto& it = _Corrections.find(_RunNumber);
							if (it != _Corrections.end()) {
								if        (_ECALClusterIndices[i] == 1) {
									correction = it->second.first;
								} else if (_ECALClusterIndices[i] == 2) {
									correction = it->second.second;
								} else {
									std::cerr << "ECAL index " << _ECALClusterIndices[i] << " is neither 1 nor 2." ;
									return false;
								}
							}
							_ResultCorrectedEnergies[i] = correction * _Energies[i];
						}
						return true;
					}

				private:

					const std::vector<double>&                     _Energies;
					const std::vector<int>&                        _ECALClusterIndices;
					const int&                                     _RunNumber;
					const std::map<int, std::pair<double, double>> _Corrections;  // constant parameter, needs to be copied
					std::vector<double>&                           _ResultCorrectedEnergies;

				};


				// taken from Sebastian Uhl's analysis of pi-pi0pi0
				// http://wwwcompass.cern.ch/compass/publications/theses/2016_phd_uhl.pdf
				// see line 505ff in /nfs/freenas/tuph/e18/project/compass/analysis/suhl/scripts/FinalState_3pi.-00/KinematicPlots.C
				class GetECALCorrectedTiming : public Function
				{

				public:

					GetECALCorrectedTiming(const std::vector<double>&                        Times,                 // times of ECAL clusters to be corrected
					                       const std::vector<double>&                        Energies,              // energies of ECAL clusters
					                       const std::vector<int>&                           ECALClusterIndices,    // ECAL indices of clusters
					                       const std::map<std::string, std::vector<double>>& Calibration,           // calibration coefficients used to correct times
					                       std::vector<double>&                              ResultCorrectedTimes)  // corrected times of ECAL clusters
						: _Times               (Times),
						  _Energies            (Energies),
						  _ECALClusterIndices  (ECALClusterIndices),
						  _Calibration         (Calibration),
						  _ResultCorrectedTimes(ResultCorrectedTimes)
					{ }

					virtual ~GetECALCorrectedTiming() { }

					bool
					operator() ()
					{
						const size_t nmbClusters = _Times.size();
						if ((_Energies.size() != nmbClusters) or (_ECALClusterIndices.size() != nmbClusters)) {
							std::cerr << "Input vectors do not have the same size." << std::endl;
							return false;
						}

						_ResultCorrectedTimes.resize(nmbClusters);
						for (size_t i = 0; i < nmbClusters; ++i) {
							const double energy  = _Energies[i];
							const double energy2 = energy * energy;
							double correction = 0;
							if        (_ECALClusterIndices[i] == 1) {
								const std::vector<double>& coefficients = _Calibration.at("ECAL1");
								correction = coefficients[0] + coefficients[1] / energy  + coefficients[2] * energy
								                             + coefficients[3] / energy2 + coefficients[4] * energy2;
							} else if (_ECALClusterIndices[i] == 2) {
							  const double energy3 = energy * energy2;
								const std::vector<double>& coefficients = _Calibration.at("ECAL2");
								correction = coefficients[0] + coefficients[1] / energy  + coefficients[2] * energy
								                             + coefficients[3] / energy2 + coefficients[4] * energy2
								                             + coefficients[5] / energy3 + coefficients[6] * energy3;
							} else {
								std::cerr << "ECAL index " << _ECALClusterIndices[i] << " is neither 1 nor 2." ;
								return false;
							}
							_ResultCorrectedTimes[i] = _Times[i] - correction;
						}
						return true;
					}

				private:

					const std::vector<double>&                       _Times;
					const std::vector<double>&                       _Energies;
					const std::vector<int>&                          _ECALClusterIndices;
					const std::map<std::string, std::vector<double>> _Calibration;  // constant parameter, needs to be copied
					std::vector<double>&                             _ResultCorrectedTimes;

				};


				class GetCleanedEcalClusters : public Function
				{

				public:

					GetCleanedEcalClusters(const std::vector<TVector3>&                      Positions,                 // positions of ECAL clusters
					                       const std::vector<TVector3>&                      PositionVariances,         // position variances of ECAL clusters
					                       const std::vector<double>&                        Energies,                  // energies of ECAL clusters
					                       const std::vector<double>&                        EnergyVariances,           // energy variances of ECAL clusters
					                       const std::vector<double>&                        Times,                     // times of ECAL clusters
					                       const std::vector<int>&                           ECALClusterIndices,        // ECAL indices of clusters
					                       const double&                                     ECAL1ThresholdEnergy,      // energy threshold applied to ECAL1 clusters
					                       const double&                                     ECAL2ThresholdEnergy,      // energy threshold applied to ECAL2 clusters
					                       const std::map<std::string, std::vector<double>>& ResolutionCoeffs,          // coefficients used to parametrize energy dependence of time resolution
					                       std::vector<TVector3>&                            ResultPositions,           // positions of ECAL clusters
					                       std::vector<TVector3>&                            ResultPositionVariances,   // position variances of ECAL clusters
					                       std::vector<double>&                              ResultEnergies,            // energies of ECAL clusters
					                       std::vector<double>&                              ResultEnergyVariances,     // energy variances of ECAL clusters
					                       std::vector<double>&                              ResultTimes,               // times of ECAL clusters
					                       std::vector<int>&                                 ResultECALClusterIndices)  // indices of the ECAL that measured the photons
						: _Positions               (Positions),
						  _PositionVariances       (PositionVariances),
						  _Energies                (Energies),
						  _EnergyVariances         (EnergyVariances),
						  _Times                   (Times),
						  _ECALClusterIndices      (ECALClusterIndices),
						  _ECAL1ThresholdEnergy    (ECAL1ThresholdEnergy),
						  _ECAL2ThresholdEnergy    (ECAL2ThresholdEnergy),
						  _ResolutionCoeffs        (ResolutionCoeffs),
						  _ResultPositions         (ResultPositions),
						  _ResultPositionVariances (ResultPositionVariances),
						  _ResultEnergies          (ResultEnergies),
						  _ResultEnergyVariances   (ResultEnergyVariances),
						  _ResultTimes             (ResultTimes),
						  _ResultECALClusterIndices(ResultECALClusterIndices)
					{ }

					virtual ~GetCleanedEcalClusters() { }

					bool
					operator() ()
					{
						const size_t nmbClusters = _Positions.size();
						if (   (_PositionVariances.size()  != nmbClusters)
						    or (_Energies.size()           != nmbClusters)
						    or (_EnergyVariances.size()    != nmbClusters)
						    or (_Times.size()              != nmbClusters)
						    or (_ECALClusterIndices.size() != nmbClusters)) {
							std::cerr << "Input vectors do not have the same size." << std::endl;
							return false;
						}

						_ResultPositions.clear();
						_ResultPositionVariances.clear();
						_ResultEnergies.clear();
						_ResultEnergyVariances.clear();
						_ResultTimes.clear();
						_ResultECALClusterIndices.clear();
						_ResultPositions.reserve         (nmbClusters);
						_ResultPositionVariances.reserve (nmbClusters);
						_ResultEnergies.reserve          (nmbClusters);
						_ResultEnergyVariances.reserve   (nmbClusters);
						_ResultTimes.reserve             (nmbClusters);
						_ResultECALClusterIndices.reserve(nmbClusters);
						for (size_t i = 0; i < nmbClusters; ++i) {
							const double energy  = _Energies[i];
							const double energy2 = energy * energy;
							if (_ECALClusterIndices[i] == 1) {
								// apply energy threshold
								if (energy < _ECAL1ThresholdEnergy) {
									continue;
								}
								// apply energy-dependent cut on cluster time
								// taken from Sebastian Uhl's analysis of pi-pi0pi0
								// http://wwwcompass.cern.ch/compass/publications/theses/2016_phd_uhl.pdf
								// see lines 553ff in /nfs/freenas/tuph/e18/project/compass/analysis/suhl/scripts/FinalState_3pi.-00/KinematicPlots.C
								const std::vector<double>& coefficients = _ResolutionCoeffs.at("ECAL1");
								const double sigmaT = std::sqrt(coefficients[0] + coefficients[1] / energy
								                                                + coefficients[2] / energy2);
								if (fabs(_Times[i]) > 3 * sigmaT) {
									continue;
								}
								_ResultECALClusterIndices.push_back(1);
							} else if (_ECALClusterIndices[i] == 2) {
								// apply energy threshold
								if (energy < _ECAL2ThresholdEnergy) {
									continue;
								}
								// apply energy-dependent cut on cluster time; see above
								const std::vector<double>& coefficients = _ResolutionCoeffs.at("ECAL2");
								const double sigmaT = std::sqrt(coefficients[0] + coefficients[1] / energy  + coefficients[2] * energy
								                                                + coefficients[3] / energy2 + coefficients[4] * energy2);
								if (fabs(_Times[i]) > 3 * sigmaT) {
									continue;
								}
								_ResultECALClusterIndices.push_back(2);
							} else {
								std::cerr << "ECAL index " << _ECALClusterIndices[i] << " is neither 1 nor 2." ;
								return false;
							}
							_ResultPositions.push_back        (_Positions[i]);
							_ResultPositionVariances.push_back(_PositionVariances[i]);
							_ResultEnergies.push_back         (energy);
							_ResultEnergyVariances.push_back  (_EnergyVariances[i]);
							_ResultTimes.push_back            (_Times[i]);
						}
						return true;
					}

				private:

					const std::vector<TVector3>&                     _Positions;
					const std::vector<TVector3>&                     _PositionVariances;
					const std::vector<double>&                       _Energies;
					const std::vector<double>&                       _EnergyVariances;
					const std::vector<double>&                       _Times;
					const std::vector<int>&                          _ECALClusterIndices;
					const double                                     _ECAL1ThresholdEnergy;  // constant parameter, needs to be copied
					const double                                     _ECAL2ThresholdEnergy;  // constant parameter, needs to be copied
					const std::map<std::string, std::vector<double>> _ResolutionCoeffs;      // constant parameter, needs to be copied
					std::vector<TVector3>&                           _ResultPositions;
					std::vector<TVector3>&                           _ResultPositionVariances;
					std::vector<double>&                             _ResultEnergies;
					std::vector<double>&                             _ResultEnergyVariances;
					std::vector<double>&                             _ResultTimes;
					std::vector<int>&                                _ResultECALClusterIndices;

				};


				class getECALVariables : public Function
				{

				public:

					getECALVariables(const std::vector<int>&            ECALClusterIndices,        // ECAL indices of cluster
					                 const std::vector<TVector3>&       Positions,                 // positions of ECAL clusters
					                 const std::vector<TVector3>&       PositionVariances,         // variances in position of ECAL clusters
					                 const std::vector<double>&         Energies,                  // energies of ECAL clusters
					                 const std::vector<double>&         EnergyVariances,           // variance in energy of ECAL clusters
					                 const std::vector<double>&         Times,                     // times of ECAL clusters
					                 const int&                         SelectedECALIndex,         // index of selected ECAL
					                 std::vector<int>&                  ResultECALClusterIndices,  // ECAL indices of cluster
					                 std::vector<TVector3>&             ResultPositions,           // positions of ECAL clusters
					                 std::vector<TVector3>&             ResultPositionVariances,   // variances in position of ECAL clusters
					                 std::vector<double>&               ResultEnergies,            // energies of ECAL clusters
					                 std::vector<double>&               ResultEnergyVariances,     // variance in energy of ECAL clusters
					                 std::vector<double>&               ResultTimes)               // times of ECAL clusters
						: _ECALClusterIndices      (ECALClusterIndices),
						  _Positions               (Positions),
						  _PositionVariances       (PositionVariances),
						  _Energies                (Energies),
						  _EnergyVariances         (EnergyVariances),
						  _Times                   (Times),
						  _SelectedECALIndex       (SelectedECALIndex),
						  _ResultECALClusterIndices(ResultECALClusterIndices),
						  _ResultPositions         (ResultPositions),
						  _ResultPositionVariances (ResultPositionVariances),
						  _ResultEnergies          (ResultEnergies),
						  _ResultEnergyVariances   (ResultEnergyVariances),
						  _ResultTimes             (ResultTimes)
					{ }

					virtual ~getECALVariables() { }

					bool
					operator() ()
					{
						// check if all vectors have the same size
						const size_t nmbClusters = _ECALClusterIndices.size();
						if (   (_Positions.size()         != nmbClusters)
						    or (_PositionVariances.size() != nmbClusters)
						    or (_Energies.size()          != nmbClusters)
						    or (_EnergyVariances.size()   != nmbClusters)
						    or (_Times.size()             != nmbClusters)) {
							std::cerr << "Input vectors do not have the same size." << std::endl;
							return false;
						}

						// get number of clusters in given ECAL
						size_t nmbClustersResult = 0;
						for (size_t i = 0; i < nmbClusters; ++i) {
							if (_ECALClusterIndices[i] == _SelectedECALIndex) {
								++nmbClustersResult;
							}
						}
						_ResultECALClusterIndices.resize(nmbClustersResult);
						_ResultPositions.resize         (nmbClustersResult);
						_ResultPositionVariances.resize (nmbClustersResult);
						_ResultEnergies.resize          (nmbClustersResult);
						_ResultEnergyVariances.resize   (nmbClustersResult);
						_ResultTimes.resize             (nmbClustersResult);

						// loop over clusters and save information if cluster is in required ECAL
						int nmbAccepted = 0;
						for (size_t i = 0; i < nmbClusters; ++i) {
							if (_ECALClusterIndices[i] == _SelectedECALIndex) {
								_ResultECALClusterIndices[nmbAccepted] = _ECALClusterIndices[i];
								_ResultPositions         [nmbAccepted] = _Positions         [i];
								_ResultPositionVariances [nmbAccepted] = _PositionVariances [i];
								_ResultEnergies          [nmbAccepted] = _Energies          [i];
								_ResultEnergyVariances   [nmbAccepted] = _EnergyVariances   [i];
								_ResultTimes             [nmbAccepted] = _Times             [i];
								++nmbAccepted;
							}
						}
						return true;
					}

				private:

					const std::vector<int>&      _ECALClusterIndices;
					const std::vector<TVector3>& _Positions;
					const std::vector<TVector3>& _PositionVariances;
					const std::vector<double>&   _Energies;
					const std::vector<double>&   _EnergyVariances;
					const std::vector<double>&   _Times;
					const int                    _SelectedECALIndex;  // const paramerter, needs to be copied
					std::vector<int>&            _ResultECALClusterIndices;
					std::vector<TVector3>&       _ResultPositions;
					std::vector<TVector3>&       _ResultPositionVariances;
					std::vector<double>&         _ResultEnergies;
					std::vector<double>&         _ResultEnergyVariances;
					std::vector<double>&         _ResultTimes;

				};


				class GetPhotonLorentzVecs : public Function
				{

				public:

					GetPhotonLorentzVecs(const std::vector<TVector3>& Positions,           // positions of ECAL clusters
					                     const std::vector<double>&   Energies,            // energies of ECAL clusters
					                     const std::vector<double>&   Times,               // times of ECAL clusters
					                     const std::vector<int>&      ECALClusterIndices,  // ECAL indices of clusters
					                     const TVector3&              VertexPosition,      // position of primary vertex
					                     std::vector<TLorentzVector>& ResultPhotonLVs)     // photon Lorentz vectors
						: _Positions         (Positions),
						  _Energies          (Energies),
						  _Times             (Times),
						  _ECALClusterIndices(ECALClusterIndices),
						  _VertexPosition    (VertexPosition),
						  _ResultPhotonLVs   (ResultPhotonLVs)
					{ }

					virtual ~GetPhotonLorentzVecs() { }

					bool
					operator() ()
					{
						const size_t nmbPhotons = _Positions.size();
						if (   (_Energies.size()           != nmbPhotons)
						    or (_Times.size()              != nmbPhotons)
						    or (_ECALClusterIndices.size() != nmbPhotons)) {
							std::cerr << "Input vectors do not have the same size." << std::endl;
							return false;
						}

						_ResultPhotonLVs.resize(nmbPhotons);
						for (size_t i = 0; i < nmbPhotons; ++i) {
							TVector3 photonMom = _Positions[i] - _VertexPosition;
							photonMom.SetMag(_Energies[i]);
							_ResultPhotonLVs[i] = TLorentzVector(photonMom, _Energies[i]);
						}
						return true;
					}

				private:

					const std::vector<TVector3>& _Positions;
					const std::vector<double>&   _Energies;
					const std::vector<double>&   _Times;
					const std::vector<int>&      _ECALClusterIndices;
					const TVector3&              _VertexPosition;
					std::vector<TLorentzVector>& _ResultPhotonLVs;

				};


				namespace {

					double
					getECALMassWindow(const int    ECALIndex1,
					                  const int    ECALIndex2,
					                  const double ECAL1MassWindow,
					                  const double ECAL2MassWindow,
					                  const double ECALMixedMassWindow)
					{
						if        ((ECALIndex1 == 1) and (ECALIndex2 == 1)) {
							return ECAL1MassWindow;
						} else if ((ECALIndex1 == 2) and (ECALIndex2 == 2)) {
							return ECAL2MassWindow;
						} else if (   ((ECALIndex1 == 1) and (ECALIndex2 == 2))
						           or ((ECALIndex1 == 2) and (ECALIndex2 == 1))) {
							return ECALMixedMassWindow;
						}
						std::stringstream errMsg;
						errMsg << "At least one of the given ECAL indices (" << ECALIndex1 << ", " << ECALIndex2 << ") is unknown. Aborting...";
						throw std::runtime_error(errMsg.str());
					}

				}


				class GetPhotonPairParticles : public Function
				{

				public:

					GetPhotonPairParticles(const std::vector<TLorentzVector>& PhotonLVs,            // Lorentz vectors of photons
					                       const std::vector<int>&            ECALClusterIndices,   // indices of the ECAL that measured the photons
					                       const double&                      NominalMass,          // nominal mass of photon pair
					                       const double&                      ECALMixedMassWindow,  // m(gamma gamma) cut applied around NominalMass when the photon pair is in ECAL1 and 2
					                       const double&                      ECAL1MassWindow,      // m(gamma gamma) cut applied around NominalMass when both photons are in ECAL1
					                       const double&                      ECAL2MassWindow,      // m(gamma gamma) cut applied around NominalMass when both photons are in ECAL2
					                       std::vector<TLorentzVector>&       ResultParticleLVs,    // Lorentz vectors of all particles reconstructed from photon pairs
					                       int&                               ResultHasParticles)   // 1 if at least one photon pair was found; else 0
						: _PhotonLVs          (PhotonLVs),
						  _ECALClusterIndices (ECALClusterIndices),
						  _NominalMass        (NominalMass),
						  _ECALMixedMassWindow(ECALMixedMassWindow),
						  _ECAL1MassWindow    (ECAL1MassWindow),
						  _ECAL2MassWindow    (ECAL2MassWindow),
						  _ResultParticleLVs  (ResultParticleLVs),
						  _ResultHasParticles (ResultHasParticles)
					{ }

					virtual ~GetPhotonPairParticles() { }

					bool
					operator() ()
					{
						const size_t nmbPhotons = _PhotonLVs.size();
						if (_ECALClusterIndices.size() != nmbPhotons) {
							std::cerr << "Input vectors do not have the same size." << std::endl;
							return false;
						}

						_ResultParticleLVs.clear();
						_ResultParticleLVs.reserve(nmbPhotons * nmbPhotons);
						_ResultHasParticles = 0;
						if (nmbPhotons < 2) {
							return true;
						}
						for (size_t i = 0; i < nmbPhotons; ++i) {
							for (size_t j = i + 1; j < nmbPhotons; ++j) {
								const TLorentzVector candidate = _PhotonLVs[i] + _PhotonLVs[j];
								if (std::fabs(candidate.M() - _NominalMass)
								    < getECALMassWindow(_ECALClusterIndices[i], _ECALClusterIndices[j], _ECAL1MassWindow, _ECAL2MassWindow, _ECALMixedMassWindow)) {
									_ResultParticleLVs.push_back(candidate);
								}
							}
						}
						if (_ResultParticleLVs.size() > 1) {
							_ResultHasParticles = 1;
						}
						return true;
					}

				private:

					const std::vector<TLorentzVector>& _PhotonLVs;
					const std::vector<int>&            _ECALClusterIndices;
					const double                       _NominalMass;          // constant parameter, needs to be copied
					const double                       _ECALMixedMassWindow;  // constant parameter, needs to be copied
					const double                       _ECAL1MassWindow;      // constant parameter, needs to be copied
					const double                       _ECAL2MassWindow;      // constant parameter, needs to be copied
					std::vector<TLorentzVector>&       _ResultParticleLVs;
					int&                               _ResultHasParticles;

				};


				class GetPi0Pair : public Function
				{

				public:

					GetPi0Pair(const std::vector<TLorentzVector>& PhotonLVs,                 // Lorentz vectors of photons
					           const std::vector<int>&            ECALClusterIndices,        // indices of the ECAL that measured the photons
					           const double&                      Pi0Mass,                   // nominal pi^0 mass
					           const double&                      ECALMixedMassWindow,       // m(gamma gamma) cut applied around Pi0Mass when the photon pair is in ECAL1 and 2
					           const double&                      ECAL1MassWindow,           // m(gamma gamma) cut applied around Pi0Mass when both photons are in ECAL1
					           const double&                      ECAL2MassWindow,           // m(gamma gamma) cut applied around Pi0Mass when both photons are in ECAL2
					           std::vector<TLorentzVector>&       ResultPi0PairLVs,          // Lorentz vectors of the two pi^0 in the first found pair
					           TLorentzVector&                    ResultPi0LV_0,             // Lorentz vectors of the 1st pi^0 in the first found pair
					           TLorentzVector&                    ResultPi0LV_1,             // Lorentz vectors of the 2nd pi^0 in the first found pair
					           int&                               ResultGoodPi0Pair,         // 1 if exactly one pi^0 pair was found; else 0
					           std::vector<int>&                  ResultECALClusterIndices)  // indices of the ECAL that measured the selected photons
						: _PhotonLVs                 (PhotonLVs),
						  _ECALClusterIndices        (ECALClusterIndices),
						  _Pi0Mass                   (Pi0Mass),
						  _ECALMixedMassWindow       (ECALMixedMassWindow),
						  _ECAL1MassWindow           (ECAL1MassWindow),
						  _ECAL2MassWindow           (ECAL2MassWindow),
						  _ResultPi0PairLVs          (ResultPi0PairLVs),
						  _ResultPi0LV_0             (ResultPi0LV_0),
						  _ResultPi0LV_1             (ResultPi0LV_1),
						  _ResultGoodPi0Pair         (ResultGoodPi0Pair),
						  _ResultECALClusterIndices  (ResultECALClusterIndices)
					{ }

					virtual ~GetPi0Pair() { }

					bool
					operator() ()
					{
						const size_t nmbPhotons = _PhotonLVs.size();
						if (_ECALClusterIndices.size() != nmbPhotons) {
							std::cerr << "Input vectors do not have the same size." << std::endl;
							return false;
						}

						_ResultPi0PairLVs.clear();
						_ResultPi0PairLVs.reserve(2);
						_ResultPi0LV_0.SetPxPyPzE(0, 0, 0, 0);
						_ResultPi0LV_1.SetPxPyPzE(0, 0, 0, 0);
						_ResultGoodPi0Pair = 0;
						_ResultECALClusterIndices = {-1, -1, -1, -1};
						if (nmbPhotons < 4) {
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
								// first photon pair
								const TLorentzVector pi0Candidate0 = _PhotonLVs[i] + _PhotonLVs[j];
								if (std::fabs(pi0Candidate0.M() - _Pi0Mass)
								    > getECALMassWindow(_ECALClusterIndices[i], _ECALClusterIndices[j], _ECAL1MassWindow, _ECAL2MassWindow, _ECALMixedMassWindow)) {
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
										// second photon pair
										if ((m == j) or (n == j)) {
											continue;
										}
										const TLorentzVector pi0Candidate1 = _PhotonLVs[m] + _PhotonLVs[n];
										if (std::fabs(pi0Candidate1.M() - _Pi0Mass)
										    > getECALMassWindow(_ECALClusterIndices[m], _ECALClusterIndices[n], _ECAL1MassWindow, _ECAL2MassWindow, _ECALMixedMassWindow)) {
											continue;
										}
										if (nmbCandidatePairs == 0) {
											_ResultPi0PairLVs.push_back(pi0Candidate0);
											_ResultPi0PairLVs.push_back(pi0Candidate1);
											_ResultPi0LV_0            = pi0Candidate0;
											_ResultPi0LV_1            = pi0Candidate1;
											_ResultECALClusterIndices = {(int)i, (int)j, (int)m, (int)n};
										}
										nmbCandidatePairs++;
									}
								}
							}
						}
						if (nmbCandidatePairs == 1) {
							_ResultGoodPi0Pair = 1;
						}

						return true;
					}

				private:

					const std::vector<TLorentzVector>& _PhotonLVs;
					const std::vector<int>&            _ECALClusterIndices;
					const double                       _Pi0Mass;              // constant parameter, needs to be copied
					const double                       _ECALMixedMassWindow;  // constant parameter, needs to be copied
					const double                       _ECAL1MassWindow;      // constant parameter, needs to be copied
					const double                       _ECAL2MassWindow;      // constant parameter, needs to be copied
					std::vector<TLorentzVector>&       _ResultPi0PairLVs;
					TLorentzVector&                    _ResultPi0LV_0;
					TLorentzVector&                    _ResultPi0LV_1;
					int&                               _ResultGoodPi0Pair;
					std::vector<int>&                  _ResultECALClusterIndices;

				};


				class GetKinematicFittingMass : public Function
				{

				public:

					GetKinematicFittingMass(const TVector3&              VertexPosition,
					                        const std::vector<TVector3>& ClusterPositions,
					                        const std::vector<TVector3>& ClusterPositionVariances,
					                        const std::vector<double>&   ClusterEnergies,
					                        const std::vector<double>&   ClusterEnergyVariances,
					                        const std::vector<int>&      ClusterECALIndices,
					                        const double&                Mass,
					                        const double&                MassLowerLimit,
					                        const double&                MassUpperLimit,
					                        const double&                PrecisionGoal,
					                        const int&                   WhichEnergyVariance,
					                        std::vector<TLorentzVector>& ResultLorentzVectors,
					                        std::vector<double>&         ResultChi2s,
					                        std::vector<double>&         ResultPValues,
					                        std::vector<int>&            ResultNmbIterations,
					                        int&                         ResultSuccess,
					                        std::vector<double>&         ResultPullsX0,
					                        std::vector<double>&         ResultPullsY0,
					                        std::vector<double>&         ResultPullsE0,
					                        std::vector<double>&         ResultPullsX1,
					                        std::vector<double>&         ResultPullsY1,
					                        std::vector<double>&         ResultPullsE1)
						: _VertexPosition          (VertexPosition),
						  _ClusterPositions        (ClusterPositions),
						  _ClusterPositionVariances(ClusterPositionVariances),
						  _ClusterEnergies         (ClusterEnergies),
						  _ClusterEnergyVariances  (ClusterEnergyVariances),
						  _ClusterECALIndices      (ClusterECALIndices),
						  _Mass                    (Mass),
						  _MassLowerLimit          (MassLowerLimit),
						  _MassUpperLimit          (MassUpperLimit),
						  _PrecisionGoal           (PrecisionGoal),
						  _WhichEnergyVariance     (WhichEnergyVariance),
						  _ResultLorentzVectors    (ResultLorentzVectors),
						  _ResultChi2s             (ResultChi2s),
						  _ResultPValues           (ResultPValues),
						  _ResultNmbIterations     (ResultNmbIterations),
						  _ResultSuccess           (ResultSuccess),
						  _ResultPullsX0           (ResultPullsX0),
						  _ResultPullsY0           (ResultPullsY0),
						  _ResultPullsE0           (ResultPullsE0),
						  _ResultPullsX1           (ResultPullsX1),
						  _ResultPullsY1           (ResultPullsY1),
						  _ResultPullsE1           (ResultPullsE1)
					{ }

					virtual ~GetKinematicFittingMass() { }

					bool
					operator() ()
					{
						if (_ClusterECALIndices.size() != 4) {
							std::cerr << "Number of cluster indices " << _ClusterECALIndices.size() << " != 4." << std::endl;
							return false;
						}
						const size_t nmbClusters = _ClusterPositions.size();
						if (   (_ClusterPositionVariances.size() != nmbClusters)
						    or (_ClusterEnergies.size()          != nmbClusters)
						    or (_ClusterEnergyVariances.size()   != nmbClusters)) {
							std::cerr << "Input ECAL vectors do not have the same size." << std::endl;
							return false;
						}

						_ResultLorentzVectors.clear();
						_ResultChi2s.clear();
						_ResultPValues.clear();
						_ResultNmbIterations.clear();
						_ResultPullsX0.clear();
						_ResultPullsY0.clear();
						_ResultPullsE0.clear();
						_ResultPullsX1.clear();
						_ResultPullsY1.clear();
						_ResultPullsE1.clear();
						_ResultLorentzVectors.resize(2);
						_ResultChi2s.resize         (2);
						_ResultPValues.resize       (2);
						_ResultNmbIterations.resize (2);
						_ResultPullsX0.resize       (2);
						_ResultPullsY0.resize       (2);
						_ResultPullsE0.resize       (2);
						_ResultPullsX1.resize       (2);
						_ResultPullsY1.resize       (2);
						_ResultPullsE1.resize       (2);
						_ResultSuccess = 0;

						if (   (_ClusterECALIndices[0] == -1)
						    or (_ClusterECALIndices[1] == -1)
						    or (_ClusterECALIndices[2] == -1)
						    or (_ClusterECALIndices[3] == -1)) {
							return true;
						}
						const double maxPosError = 1e3;
						for (size_t i = 0; i < 4; ++i) {
							if (   (_ClusterPositionVariances[_ClusterECALIndices[i]].X() > maxPosError)
							    or (_ClusterPositionVariances[_ClusterECALIndices[i]].Y() > maxPosError)
							    or (_ClusterPositionVariances[_ClusterECALIndices[i]].Z() > maxPosError)) {
								return true;
							}
						}
						// Fit cluster 0 and 1 and cluster 2 and 3 to Pi0 mass
						std::vector<bool> successes;
						for (size_t i = 0; i < 2; ++i) {
							antok::NeutralFit neutralFit(_VertexPosition,
							                             _ClusterPositions        [_ClusterECALIndices[i*2]],
							                             _ClusterPositions        [_ClusterECALIndices[i*2+1]],
							                             _ClusterPositionVariances[_ClusterECALIndices[i*2]],
							                             _ClusterPositionVariances[_ClusterECALIndices[i*2+1]],
							                             _ClusterEnergies         [_ClusterECALIndices[i*2]],
							                             _ClusterEnergies         [_ClusterECALIndices[i*2+1]],
							                             _ClusterEnergyVariances  [_ClusterECALIndices[i*2]],
							                             _ClusterEnergyVariances  [_ClusterECALIndices[i*2+1]],
							                             _Mass,
							                             _MassLowerLimit,
							                             _MassUpperLimit,
							                             _PrecisionGoal,
							                             _WhichEnergyVariance);
							successes.push_back(neutralFit.doFit());
							if (successes[i]) {
								_ResultLorentzVectors[i] = neutralFit.getImprovedLVSum();
								_ResultPullsX0       [i] = neutralFit.pullValues()[0];
								_ResultPullsY0       [i] = neutralFit.pullValues()[1];
								_ResultPullsE0       [i] = neutralFit.pullValues()[2];
								_ResultPullsX1       [i] = neutralFit.pullValues()[3];
								_ResultPullsY1       [i] = neutralFit.pullValues()[4];
								_ResultPullsE1       [i] = neutralFit.pullValues()[5];
								_ResultChi2s         [i] = neutralFit.chi2Value();
								_ResultPValues       [i] = neutralFit.pValue();
								_ResultNmbIterations [i] = neutralFit.nmbIterations();
							}
						}
						if (successes[0] and successes[1]) {
							_ResultSuccess = 1;
						}
						return true;
					}

				private:

					const TVector3&              _VertexPosition;
					const std::vector<TVector3>& _ClusterPositions;
					const std::vector<TVector3>& _ClusterPositionVariances;
					const std::vector<double>&   _ClusterEnergies;
					const std::vector<double>&   _ClusterEnergyVariances;
					const std::vector<int>&      _ClusterECALIndices;
					const double                 _Mass;                 // constant parameter, needs to be copied
					const double                 _MassLowerLimit;       // constant parameter, needs to be copied
					const double                 _MassUpperLimit;       // constant parameter, needs to be copied
					const double                 _PrecisionGoal;        // constant parameter, needs to be copied
					const int                    _WhichEnergyVariance;  // constant parameter, needs to be copied
					std::vector<TLorentzVector>& _ResultLorentzVectors;
					std::vector<double>&         _ResultChi2s;
					std::vector<double>&         _ResultPValues;
					std::vector<int>&            _ResultNmbIterations;
					int&                         _ResultSuccess;
					std::vector<double>&         _ResultPullsX0;
					std::vector<double>&         _ResultPullsY0;
					std::vector<double>&         _ResultPullsE0;
					std::vector<double>&         _ResultPullsX1;
					std::vector<double>&         _ResultPullsY1;
					std::vector<double>&         _ResultPullsE1;

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
					         const double&         OmegaMassWindow,         // cut around OmegaMass applied on m(pi^- pi^0 pi^+)
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
						  _OmegaMassWindow       (OmegaMassWindow),
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
										if (std::fabs(candidate.M() - _OmegaMass) < _OmegaMassWindow) {
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
					const double          _OmegaMassWindow;  // constant parameter, needs to be copied
					TLorentzVector&       _ResultOmegaLV;
					int&                  _ResultAccepted;
					TLorentzVector&       _ResultNotUsedPi0LV;
					TLorentzVector&       _ResultNotUsedPiMinusLV;

				};


				class GetFittedOmegaMassVsPrecisionGoal : public Function
				{

				public:

					GetFittedOmegaMassVsPrecisionGoal(const TVector3&              VertexPosition,
					                                  const TLorentzVector&        ChargedPartLV_0,  // Lorentz vector of 1st charged particle
					                                  const TLorentzVector&        ChargedPartLV_1,  // Lorentz vector of 2nd charged particle
					                                  const TLorentzVector&        ChargedPartLV_2,  // Lorentz vector of 3rd charged particle
					                                  const int&                   Charge_0,         // charge of 1st charged particle
					                                  const int&                   Charge_1,         // charge of 2nd charged particle
					                                  const int&                   Charge_2,         // charge of 3rd charged particle
					                                  const std::vector<TVector3>& ClusterPositions,
					                                  const std::vector<TVector3>& ClusterPositionVariances,
					                                  const std::vector<double>&   ClusterEnergies,
					                                  const std::vector<double>&   ClusterEnergyVariances,
					                                  const std::vector<int>&      ClusterECALIndices,
					                                  const double&                PiMass,
					                                  const double&                PiMassLowerLimit,
					                                  const double&                PiMassUpperLimit,
					                                  const double&                precisionGoalLowerLimit,
					                                  const double&                precisionGoalUpperLimit,
					                                  const int&                   whichEnergyVariance,
					                                  const double&                OmegaMass,
					                                  const double&                OmegaMasswindow,
					                                  std::vector<double>&         ResultPrecisionGoals,
					                                  std::vector<int>&            ResultAcceptedOmegas,
					                                  std::vector<double>&         ResultOmegaMasses)
						: _VertexPosition          (VertexPosition),
						  _ChargedPartLV_0         (ChargedPartLV_0),
						  _ChargedPartLV_1         (ChargedPartLV_1),
						  _ChargedPartLV_2         (ChargedPartLV_2),
						  _Charge_0                (Charge_0),
						  _Charge_1                (Charge_1),
						  _Charge_2                (Charge_2),
						  _ClusterPositions        (ClusterPositions),
						  _ClusterPositionVariances(ClusterPositionVariances),
						  _ClusterEnergies         (ClusterEnergies),
						  _ClusterEnergyVariances  (ClusterEnergyVariances),
						  _ClusterECALIndices      (ClusterECALIndices),
						  _PiMass                  (PiMass),
						  _PiMassLowerLimit        (PiMassLowerLimit),
						  _PiMassUpperLimit        (PiMassUpperLimit),
						  _precisionGoalLowerLimit (precisionGoalLowerLimit),
						  _precisionGoalUpperLimit (precisionGoalUpperLimit),
						  _whichEnergyVariance     (whichEnergyVariance),
						  _OmegaMass               (OmegaMass),
						  _OmegaMasswindow         (OmegaMasswindow),
						  _ResultPrecisionGoals    (ResultPrecisionGoals),
						  _ResultAcceptedOmegas    (ResultAcceptedOmegas),
						  _ResultOmegaMasses       (ResultOmegaMasses)
					{ }

					virtual ~GetFittedOmegaMassVsPrecisionGoal() { }

					bool
					operator() ()
					{
						const size_t sizeVectors = std::log(_precisionGoalUpperLimit / _precisionGoalLowerLimit);
						_ResultAcceptedOmegas.resize(sizeVectors);
						_ResultPrecisionGoals.resize(sizeVectors);
						_ResultOmegaMasses.resize   (sizeVectors);
						std::vector<TLorentzVector> PiLorentzVectors;
						std::vector<double>         PiChi2s;
						std::vector<double>         PiPullsX0;
						std::vector<double>         PiPullsY0;
						std::vector<double>         PiPullsE0;
						std::vector<double>         PiPullsX1;
						std::vector<double>         PiPullsY1;
						std::vector<double>         PiPullsE1;
						std::vector<double>         PiPValues;
						std::vector<int>            PiNmbIterations;
						int                         PiSuccess = 0;
						TLorentzVector              OmegaLV;
						int                         OmegaAccepted = 0;
						TLorentzVector              FinalPi0LV;
						TLorentzVector              FinalPiMinusLV;
						size_t index = 0;
						if (_precisionGoalLowerLimit > _precisionGoalUpperLimit) {
							std::cerr << "_precisionGoalLowerLimit is larger than _precisionGoalUpperLimit: " << _precisionGoalLowerLimit << ">" << _precisionGoalUpperLimit << ".\n";
						}
						for (double PG = _precisionGoalLowerLimit; PG <= _precisionGoalUpperLimit; PG*= 10 ) {
							GetKinematicFittingMass* PiPairFit = new GetKinematicFittingMass(_VertexPosition, _ClusterPositions, _ClusterPositionVariances, _ClusterEnergies, _ClusterEnergyVariances,
							                                                                 _ClusterECALIndices, _PiMass, _PiMassLowerLimit, _PiMassUpperLimit, PG, _whichEnergyVariance, PiLorentzVectors,
						                                                                  PiChi2s, PiPValues, PiNmbIterations, PiSuccess, PiPullsX0, PiPullsY0, PiPullsE0, PiPullsX1, PiPullsY1, PiPullsE1);
							(*PiPairFit)();
							delete PiPairFit;
							if (not PiSuccess) {
								continue;
							}
							GetOmega* getOmega = new GetOmega(PiLorentzVectors[0], PiLorentzVectors[1], _ChargedPartLV_0, _ChargedPartLV_1, _ChargedPartLV_2, _Charge_0,
							                                  _Charge_1, _Charge_2, _OmegaMass, _OmegaMasswindow, OmegaLV, OmegaAccepted, FinalPi0LV, FinalPiMinusLV);
							(*getOmega)();
							delete getOmega;
							_ResultPrecisionGoals[index] = PG;
							_ResultAcceptedOmegas[index] = OmegaAccepted;
							_ResultOmegaMasses   [index] = OmegaLV.M();
							++index;
						}
						return true;
					}

				private:

					const TVector3&              _VertexPosition;
					const TLorentzVector&        _ChargedPartLV_0;
					const TLorentzVector&        _ChargedPartLV_1;
					const TLorentzVector&        _ChargedPartLV_2;
					const int&                   _Charge_0;
					const int&                   _Charge_1;
					const int&                   _Charge_2;
					const std::vector<TVector3>& _ClusterPositions;
					const std::vector<TVector3>& _ClusterPositionVariances;
					const std::vector<double>&   _ClusterEnergies;
					const std::vector<double>&   _ClusterEnergyVariances;
					const std::vector<int>&      _ClusterECALIndices;
					const double                 _PiMass;                   // constant parameter, needs to be copied
					const double                 _PiMassLowerLimit;         // constant parameter, needs to be copied
					const double                 _PiMassUpperLimit;         // constant parameter, needs to be copied
					const double                 _precisionGoalLowerLimit;  // constant parameter, needs to be copied
					const double                 _precisionGoalUpperLimit;  // constant parameter, needs to be copied
					const int                    _whichEnergyVariance;      // constant parameter, needs to be copied
					const double                 _OmegaMass;                // constant parameter, needs to be copied
					const double                 _OmegaMasswindow;          // constant parameter, needs to be copied
					std::vector<double>&         _ResultPrecisionGoals;
					std::vector<int>&            _ResultAcceptedOmegas;
					std::vector<double>&         _ResultOmegaMasses;

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
						_Result.clear();
						_Result.reserve(4);
						const std::vector<const TLorentzVector*> Pi0LVs         = {&_Pi0LV_0, &_Pi0LV_1};
						const std::vector<const TLorentzVector*> ChargedPartLVs = {&_ChargedPartLV_0, &_ChargedPartLV_1, &_ChargedPartLV_2};
						const std::vector<int>                   Charges        = { _Charge_0,         _Charge_1,         _Charge_2};

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
