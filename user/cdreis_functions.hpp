#ifndef ANTOK_USER_CDREIS_FUNCTIONS_HPP
#define ANTOK_USER_CDREIS_FUNCTIONS_HPP

#include <iostream>
#include <limits>
#include <vector>

#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TRotation.h"

#include "functions.hpp"
#include "neutral_fit.h"


namespace antok {

	namespace user {

		namespace cdreis {

			namespace functions {

				template <typename T>
				class GetSumOverVector : public Function
				{

				public:

					GetSumOverVector(const std::vector<T>& Vector,     // input vector
					                 T&                    ResultSum)  // sum over all elements in vector
						: _Vector   (Vector),
						  _ResultSum(ResultSum)
					{ }

					virtual ~GetSumOverVector() { }

					bool
					operator() ()
					{
						_ResultSum = T();
						for (auto& entry : _Vector) {
							_ResultSum += entry;
						}
						return true;
					}

				private:

					const std::vector<T>& _Vector;
					T&                    _ResultSum;

				};


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

					GetCleanedEcalClusters(const std::vector<TVector3>&                      Positions,                   // positions of ECAL clusters
					                       const std::vector<TVector3>&                      PositionVariances,           // position variances of ECAL clusters
					                       const std::vector<double>&                        Energies,                    // energies of ECAL clusters
					                       const std::vector<double>&                        EnergyVariances,             // energy variances of ECAL clusters
					                       const std::vector<double>&                        Times,                       // times of ECAL clusters
					                       const std::vector<double>&                        DistancesToCharged,          // distances to the next charged track at z Position of the cluster
					                       const std::vector<int>&                           ECALClusterIndices,          // ECAL indices of clusters
					                       const double&                                     ECAL1ThresholdEnergy,        // energy threshold applied to ECAL1 clusters
					                       const double&                                     ECAL2ThresholdEnergy,        // energy threshold applied to ECAL2 clusters
					                       const double&                                     ECAL2YUpperLimit,            // upper limit for the Y position of ECAL2 clusters (HCAL shadow)
					                       const double&                                     ECAL2YLowerLimit,            // lower limit for the Y position of ECAL2 clusters (HCAL shadow)
					                       const double&                                     DistanceToChargedThreshold,  // minimum distance to the next charged track at z Position of the cluster
					                       const int&                                        TimeResolutionMode,          // selects parametrization for time resolution
					                       const std::map<std::string, std::vector<double>>& ResolutionCoeffs,            // coefficients used to parametrize energy dependence of time resolution
					                       std::vector<TVector3>&                            ResultPositions,             // positions of ECAL clusters
					                       std::vector<TVector3>&                            ResultPositionVariances,     // position variances of ECAL clusters
					                       std::vector<double>&                              ResultEnergies,              // energies of ECAL clusters
					                       std::vector<double>&                              ResultEnergyVariances,       // energy variances of ECAL clusters
					                       std::vector<double>&                              ResultTimes,                 // times of ECAL clusters
					                       std::vector<int>&                                 ResultECALClusterIndices)    // indices of the ECAL that measured the photons
						: _Positions                 (Positions),
						  _PositionVariances         (PositionVariances),
						  _Energies                  (Energies),
						  _EnergyVariances           (EnergyVariances),
						  _Times                     (Times),
						  _DistancesToCharged        (DistancesToCharged),
						  _ECALClusterIndices        (ECALClusterIndices),
						  _ECAL1ThresholdEnergy      (ECAL1ThresholdEnergy),
						  _ECAL2ThresholdEnergy      (ECAL2ThresholdEnergy),
						  _ECAL2YUpperLimit          (ECAL2YUpperLimit),
						  _ECAL2YLowerLimit          (ECAL2YLowerLimit),
						  _DistanceToChargedThreshold(DistanceToChargedThreshold),
						  _TimeResolutionMode        (TimeResolutionMode),
						  _ResolutionCoeffs          (ResolutionCoeffs),
						  _ResultPositions           (ResultPositions),
						  _ResultPositionVariances   (ResultPositionVariances),
						  _ResultEnergies            (ResultEnergies),
						  _ResultEnergyVariances     (ResultEnergyVariances),
						  _ResultTimes               (ResultTimes),
						  _ResultECALClusterIndices  (ResultECALClusterIndices)
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
						    or (_DistancesToCharged.size() != nmbClusters)
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
								double sigmaT;
								if        (_TimeResolutionMode == 0) {
									sigmaT = std::sqrt(coefficients[0] + coefficients[1] / energy
								                                       + coefficients[2] / energy2);
								} else if (_TimeResolutionMode == 1) {
									sigmaT = coefficients[0] + coefficients[1] / energy
								                             + coefficients[2] * energy;
								}
								if (fabs(_Times[i]) > 3 * sigmaT) {
									continue;
								}
								// apply distance to next charge threshold
								if (_DistancesToCharged[i] < _DistanceToChargedThreshold) {
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
								double sigmaT;
								if        (_TimeResolutionMode == 0) {
									sigmaT= std::sqrt(coefficients[0] + coefficients[1] / energy  + coefficients[2] * energy
								                                      + coefficients[3] / energy2 + coefficients[4] * energy2);
								} else if (_TimeResolutionMode == 1) {
									sigmaT = coefficients[0] + coefficients[1] / energy
								                             + coefficients[2] * energy;
								}
								if (fabs(_Times[i]) > 3 * sigmaT) {
									continue;
								}
								// apply HCAL shadow veto on Y position of cluster in ECAL2 (see Tobias' PhD thesis p. 62)
								if (not (_ECAL2YUpperLimit == 0 and _ECAL2YLowerLimit == 0)  // check if at least one limit is set
								    and (_Positions[i].Y() >= _ECAL2YUpperLimit or _Positions[i].Y() <= _ECAL2YLowerLimit)) {
									continue;
								}
								// apply distance to next charge threshold
								if (_DistancesToCharged[i] < _DistanceToChargedThreshold) {
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
					const std::vector<double>&                       _DistancesToCharged;
					const std::vector<int>&                          _ECALClusterIndices;
					const double                                     _ECAL1ThresholdEnergy;        // constant parameter, needs to be copied
					const double                                     _ECAL2ThresholdEnergy;        // constant parameter, needs to be copied
					const double                                     _ECAL2YUpperLimit;            // constant parameter, needs to be copied
					const double                                     _ECAL2YLowerLimit;            // constant parameter, needs to be copied
					const double                                     _DistanceToChargedThreshold;  // constant parameter, needs to be copied
					const int                                        _TimeResolutionMode;          // constant parameter, needs to be copied
					const std::map<std::string, std::vector<double>> _ResolutionCoeffs;            // constant parameter, needs to be copied
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
					                 const std::vector<double>&         DistancesToCharged,        // distances to the next charged track
					                 const int&                         SelectedECALIndex,         // index of selected ECAL
					                 std::vector<int>&                  ResultECALClusterIndices,  // ECAL indices of cluster
					                 std::vector<TVector3>&             ResultPositions,           // positions of ECAL clusters
					                 std::vector<TVector3>&             ResultPositionVariances,   // variances in position of ECAL clusters
					                 std::vector<double>&               ResultEnergies,            // energies of ECAL clusters
					                 std::vector<double>&               ResultEnergyVariances,     // variance in energy of ECAL clusters
					                 std::vector<double>&               ResultTimes,               // times of ECAL clusters
					                 std::vector<double>&               ResultDistancesToCharged)  // distances to the next charged track
						: _ECALClusterIndices      (ECALClusterIndices),
						  _Positions               (Positions),
						  _PositionVariances       (PositionVariances),
						  _Energies                (Energies),
						  _EnergyVariances         (EnergyVariances),
						  _Times                   (Times),
						  _DistancesToCharged      (DistancesToCharged),
						  _SelectedECALIndex       (SelectedECALIndex),
						  _ResultECALClusterIndices(ResultECALClusterIndices),
						  _ResultPositions         (ResultPositions),
						  _ResultPositionVariances (ResultPositionVariances),
						  _ResultEnergies          (ResultEnergies),
						  _ResultEnergyVariances   (ResultEnergyVariances),
						  _ResultTimes             (ResultTimes),
						  _ResultDistancesToCharged(ResultDistancesToCharged)
					{ }

					virtual ~getECALVariables() { }

					bool
					operator() ()
					{
						// check if all vectors have the same size
						const size_t nmbClusters = _ECALClusterIndices.size();
						if (   (_Positions.size()          != nmbClusters)
						    or (_PositionVariances.size()  != nmbClusters)
						    or (_Energies.size()           != nmbClusters)
						    or (_EnergyVariances.size()    != nmbClusters)
						    or (_Times.size()              != nmbClusters)
						    or (_DistancesToCharged.size() != nmbClusters)) {
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
						_ResultDistancesToCharged.resize(nmbClustersResult);

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
								_ResultDistancesToCharged[nmbAccepted] = _DistancesToCharged[i];
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
					const std::vector<double>&   _DistancesToCharged;
					const int                    _SelectedECALIndex;  // const paramerter, needs to be copied
					std::vector<int>&            _ResultECALClusterIndices;
					std::vector<TVector3>&       _ResultPositions;
					std::vector<TVector3>&       _ResultPositionVariances;
					std::vector<double>&         _ResultEnergies;
					std::vector<double>&         _ResultEnergyVariances;
					std::vector<double>&         _ResultTimes;
					std::vector<double>&         _ResultDistancesToCharged;

				};


				class GetPhotonLorentzVecs : public Function
				{

				public:

					GetPhotonLorentzVecs(const std::vector<TVector3>& Positions,           // positions of ECAL clusters
					                     const std::vector<double>&   Energies,            // energies of ECAL clusters
					                     const TVector3&              VertexPosition,      // position of primary vertex
					                     std::vector<TLorentzVector>& ResultPhotonLVs)     // photon Lorentz vectors
						: _Positions         (Positions),
						  _Energies          (Energies),
						  _VertexPosition    (VertexPosition),
						  _ResultPhotonLVs   (ResultPhotonLVs)
					{ }

					virtual ~GetPhotonLorentzVecs() { }

					bool
					operator() ()
					{
						const size_t nmbPhotons = _Positions.size();
						if (   (_Energies.size()           != nmbPhotons)) {
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
					const TVector3&              _VertexPosition;
					std::vector<TLorentzVector>& _ResultPhotonLVs;

				};


				class GetPhotonPairParticles : public Function
				{

				public:

					GetPhotonPairParticles(const std::vector<TLorentzVector>& PhotonLVs,               // Lorentz vectors of photons
					                       const std::vector<int>&            ECALClusterIndices,      // indices of the ECAL that measured the photons
					                       const int&                         SelectionMode,           // nominal mass of photon pair
					                       std::vector<TLorentzVector>&       ResultPhotonPairLVs,     // Lorentz vectors of all particles reconstructed from photon pairs
					                       std::vector<TLorentzVector>&       ResultPhoton0inPairLVs,  // Lorentz vectors of all first photons used in photon pairs above
					                       std::vector<TLorentzVector>&       ResultPhoton1inPairLVs,  // Lorentz vectors of all second photons used in photon pairs above
					                       std::vector<TLorentzVector>&       ResultPhotonPairsLVs_0,  // Lorentz vectors of first particle reconstructed from photon pair
					                       std::vector<TLorentzVector>&       ResultPhotonPairsLVs_1)  // Lorentz vectors of second particle reconstructed from photon pair
						: _PhotonLVs             (PhotonLVs),
						  _ECALClusterIndices    (ECALClusterIndices),
						  _SelectionMode         (SelectionMode),
						  _ResultPhotonPairLVs   (ResultPhotonPairLVs),
						  _ResultPhoton0inPairLVs(ResultPhoton0inPairLVs),
						  _ResultPhoton1inPairLVs(ResultPhoton1inPairLVs),
						  _ResultPhotonPairsLVs_0(ResultPhotonPairsLVs_0),
						  _ResultPhotonPairsLVs_1(ResultPhotonPairsLVs_1)
					{ }

					virtual ~GetPhotonPairParticles() { }

					bool
					operator() ()
					{
						enum selectionType {noSelection = -1, allSelection = 0, ECAL1Selection = 1, ECAL2Selection = 2, mixedSelection = 3};

						const size_t nmbPhotons = _PhotonLVs.size();
						if (nmbPhotons != 4) {
							return true;
						}
						if (_ECALClusterIndices.size() != nmbPhotons) {
							std::cerr << "Input vectors do not have the same size." << std::endl;
							return false;
						}
						// set selection mode to 0 if its not in {1, 2, 3}
						selectionType internSelectionMode = allSelection;
						if (   _SelectionMode == 1
							  or _SelectionMode == 2
							  or _SelectionMode == 3) {
							internSelectionMode = (selectionType)_SelectionMode;
						}

						_ResultPhotonPairLVs.clear();
						_ResultPhoton0inPairLVs.clear();
						_ResultPhoton1inPairLVs.clear();
						_ResultPhotonPairsLVs_0.clear();
						_ResultPhotonPairsLVs_1.clear();

						selectionType PhotonPairMode_0;
						selectionType PhotonPairMode_1;
						for (size_t i = 0; i < nmbPhotons - 3; ++i) {
							for (size_t j = i + 1; j < nmbPhotons; ++j) {
								// select only photon pairs with both photons in ECAL1
								if      (internSelectionMode == ECAL1Selection and (_ECALClusterIndices[i] == 1 and _ECALClusterIndices[j] == 1)) {
									PhotonPairMode_0 = ECAL1Selection;
									_ResultPhotonPairLVs.push_back(_PhotonLVs[i] + _PhotonLVs[j]);
									_ResultPhoton0inPairLVs.push_back(_PhotonLVs[i]);
									_ResultPhoton1inPairLVs.push_back(_PhotonLVs[j]);
								}
								// select only photon pairs with both photons in ECAL2
								else if (internSelectionMode == ECAL2Selection and (_ECALClusterIndices[i] == 2 or _ECALClusterIndices[j] == 2)) {
									PhotonPairMode_0 = ECAL2Selection;
									_ResultPhotonPairLVs.push_back(_PhotonLVs[i] + _PhotonLVs[j]);
									_ResultPhoton0inPairLVs.push_back(_PhotonLVs[i]);
									_ResultPhoton1inPairLVs.push_back(_PhotonLVs[j]);
								}
								// select only photon pairs with one photon in ECAL1 and one in ECAL2
								else if (internSelectionMode == mixedSelection
								         and (   (_ECALClusterIndices[i] == 1 and _ECALClusterIndices[j] == 2)
								              or (_ECALClusterIndices[i] == 2 and _ECALClusterIndices[j] == 1))) {
									PhotonPairMode_0 = mixedSelection;
									_ResultPhotonPairLVs.push_back(_PhotonLVs[i] + _PhotonLVs[j]);
									_ResultPhoton0inPairLVs.push_back(_PhotonLVs[i]);
									_ResultPhoton1inPairLVs.push_back(_PhotonLVs[j]);
								}
								// select all photon pairs
								else if (internSelectionMode == allSelection) {
									PhotonPairMode_0 = allSelection;
									_ResultPhotonPairLVs.push_back(_PhotonLVs[i] + _PhotonLVs[j]);
									_ResultPhoton0inPairLVs.push_back(_PhotonLVs[i]);
									_ResultPhoton1inPairLVs.push_back(_PhotonLVs[j]);
								}
								else {
									PhotonPairMode_0 = noSelection;
								}
								for (size_t k = i + 1; k < nmbPhotons-1; ++k) {
									if (k == j) {
										continue;
									}
									for (size_t l = k + 1; l < nmbPhotons; ++l) {
										if (l == j) {
											continue;
										}
										// select only photon pairs with both photons in ECAL1
										if      (internSelectionMode == ECAL1Selection and (_ECALClusterIndices[k] == 1 and _ECALClusterIndices[l] == 1)) {
											PhotonPairMode_1 = ECAL1Selection;
											_ResultPhotonPairLVs.push_back(_PhotonLVs[k] + _PhotonLVs[l]);
											_ResultPhoton0inPairLVs.push_back(_PhotonLVs[k]);
											_ResultPhoton1inPairLVs.push_back(_PhotonLVs[l]);
										}
										// select only photon pairs with both photons in ECAL2
										else if (internSelectionMode == ECAL2Selection and (_ECALClusterIndices[k] == 2 or _ECALClusterIndices[l] == 2)) {
											PhotonPairMode_1 = ECAL2Selection;
											_ResultPhotonPairLVs.push_back(_PhotonLVs[k] + _PhotonLVs[l]);
											_ResultPhoton0inPairLVs.push_back(_PhotonLVs[k]);
											_ResultPhoton1inPairLVs.push_back(_PhotonLVs[l]);
										}
										// select only photon pairs with one photon in ECAL1 and one in ECAL2
										else if (internSelectionMode == mixedSelection
										         and (   (_ECALClusterIndices[k] == 1 and _ECALClusterIndices[l] == 2)
										              or (_ECALClusterIndices[k] == 2 and _ECALClusterIndices[l] == 1))) {
											PhotonPairMode_1 = mixedSelection;
											_ResultPhotonPairLVs.push_back(_PhotonLVs[k] + _PhotonLVs[l]);
											_ResultPhoton0inPairLVs.push_back(_PhotonLVs[k]);
											_ResultPhoton1inPairLVs.push_back(_PhotonLVs[l]);
										}
										// select all photon pairs
										else if (internSelectionMode == allSelection) {
											PhotonPairMode_1 = allSelection;
											_ResultPhotonPairLVs.push_back(_PhotonLVs[k] + _PhotonLVs[l]);
											_ResultPhoton0inPairLVs.push_back(_PhotonLVs[k]);
											_ResultPhoton1inPairLVs.push_back(_PhotonLVs[l]);
										}
										else {
											PhotonPairMode_1 = noSelection;
										}
										// push back LV for two pairs if selected mode is true for both pairs
										// randomise order of photon pairs to undo energy ordering by PHAST
										if (    PhotonPairMode_0 == internSelectionMode
										    and PhotonPairMode_1 == internSelectionMode) {
											size_t vecSize = _ResultPhotonPairLVs.size();
											// generate random int in {0, 1} to decide if values get swapped
											const int swapMode = gRandom->Integer(2);
											// get next-to-last element if swapMode is 1, else last element
											_ResultPhotonPairsLVs_0.push_back(_ResultPhotonPairLVs[vecSize - (1 + swapMode)]);
											// get last element if swapMode is 1, else next-to-last element using modulo
											_ResultPhotonPairsLVs_1.push_back(_ResultPhotonPairLVs[vecSize - (1 + (1 + swapMode) % 2)]);
										}
									}
								}
							}
						}
						return true;
					}

				private:

					const std::vector<TLorentzVector>& _PhotonLVs;
					const std::vector<int>&            _ECALClusterIndices;
					const int                          _SelectionMode;  // constant parameter, needs to be copied
					std::vector<TLorentzVector>&       _ResultPhotonPairLVs;
					std::vector<TLorentzVector>&       _ResultPhoton0inPairLVs;
					std::vector<TLorentzVector>&       _ResultPhoton1inPairLVs;
					std::vector<TLorentzVector>&       _ResultPhotonPairsLVs_0;
					std::vector<TLorentzVector>&       _ResultPhotonPairsLVs_1;

				};


				namespace {

					void
					getECALMassParamters(const int    ECALIndex1,
					                     const int    ECALIndex2,
										 const double ECAL1Mass,
					                     const double ECAL1MassWindow,
										 const double ECAL2Mass,
					                     const double ECAL2MassWindow,
										 const double ECALMixedMass,
					                     const double ECALMixedMassWindow,
									     double &outMass,
									     double &outMassWindow)
					{
						if        ((ECALIndex1 == 1) and (ECALIndex2 == 1)) {
							outMass = ECAL1Mass;
							outMassWindow = ECAL1MassWindow;
							return;
						} else if ((ECALIndex1 == 2) and (ECALIndex2 == 2)) {
							outMass = ECAL2Mass;
							outMassWindow = ECAL2MassWindow;
							return;
						} else if (   ((ECALIndex1 == 1) and (ECALIndex2 == 2))
						           or ((ECALIndex1 == 2) and (ECALIndex2 == 1))) {
							outMass = ECALMixedMass;
							outMassWindow = ECALMixedMassWindow;
							return;
						}
						std::stringstream errMsg;
						errMsg << "At least one of the given ECAL indices (" << ECALIndex1 << ", " << ECALIndex2 << ") is unknown. Aborting...";
						throw std::runtime_error(errMsg.str());
					}

				}


				class GetPi0Pair : public Function
				{

				public:

					GetPi0Pair(const std::vector<TLorentzVector>& PhotonLVs,                    // Lorentz vectors of photons
					           const std::vector<int>&            ECALClusterIndices,           // indices of the ECAL that measured the photons
					           const double&                      ECALMixedMass,                // pi^0 mass when the photon pair is in ECAL1 and 2
					           const double&                      ECALMixedMassWindow,          // m(gamma gamma) cut applied around Pi0Mass when the photon pair is in ECAL1 and 2
					           const double&                      ECAL1Mass,                    // pi^0 mass when both photons are in ECAL1
					           const double&                      ECAL1MassWindow,              // m(gamma gamma) cut applied around Pi0Mass when both photons are in ECAL1
					           const double&                      ECAL2Mass,                    // pi^0 mass when both photons are in ECAL2
					           const double&                      ECAL2MassWindow,              // m(gamma gamma) cut applied around Pi0Mass when both photons are in ECAL2
					           std::vector<TLorentzVector>&       ResultPi0PairLVs,             // Lorentz vectors of the two pi^0 in the first found pair
							   std::vector<int>&                  ResultPi0CombTypes,           // vector of combination types of first found pair: both photons ecal1 = 1, both photons ecal2 = 2, mixed photons = 3
					           int&                               ResultNmbGoodPi0Pairs,        // 1 if exactly one pi^0 pair was found; else 0
					           std::vector<int>&                  ResultSelectedClusterIndices, // indices of the selected clusters
					           std::vector<TLorentzVector>&       ResultGammaLVsForPi0_0,       // Lorentz vectors of the two gammas in the first pi0_0
					           std::vector<TLorentzVector>&       ResultGammaLVsForPi0_1)       // Lorentz vectors of the two gammas in the first pi0_1
						: _PhotonLVs                    (PhotonLVs),
						  _ECALClusterIndices           (ECALClusterIndices),
						  _ECALMixedMass                (ECALMixedMass),
						  _ECALMixedMassWindow          (ECALMixedMassWindow),
						  _ECAL1Mass                    (ECAL1Mass),
						  _ECAL1MassWindow              (ECAL1MassWindow),
						  _ECAL2Mass                    (ECAL2Mass),
						  _ECAL2MassWindow              (ECAL2MassWindow),
						  _ResultPi0PairLVs             (ResultPi0PairLVs),
						  _ResultPi0CombTypes           (ResultPi0CombTypes),
						  _ResultNmbGoodPi0Pairs        (ResultNmbGoodPi0Pairs),
						  _ResultSelectedClusterIndices (ResultSelectedClusterIndices),
						  _ResultGammaLVsForPi0_0       (ResultGammaLVsForPi0_0),
						  _ResultGammaLVsForPi0_1       (ResultGammaLVsForPi0_1)
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
						_ResultPi0CombTypes.clear();
						_ResultPi0PairLVs.reserve(2);
						_ResultPi0CombTypes.reserve(2);
						_ResultNmbGoodPi0Pairs = 0;
						_ResultSelectedClusterIndices = {-1, -1, -1, -1};
						if (nmbPhotons < 4) {
							return true;
						}

						// searches for pairs of pi^0 candidates
						// search is aborted if more than one pair is found
						// always the first found pair is returned
						for (size_t i = 0; i < nmbPhotons - 1; ++i) {
							for (size_t j = i + 1; j < nmbPhotons; ++j) {
								// photon pair 0
								const TLorentzVector pi0Candidate0 = _PhotonLVs[i] + _PhotonLVs[j];
								const std::vector<TLorentzVector> gammasForpi0Candidate0 = {_PhotonLVs[i], _PhotonLVs[j]};
								double mass0 = 0.0;
								double massWindow0 = 0.0;
								getECALMassParamters(_ECALClusterIndices[i], _ECALClusterIndices[j], _ECAL1Mass, _ECAL1MassWindow, _ECAL2Mass, _ECAL2MassWindow, _ECALMixedMass, _ECALMixedMassWindow, mass0, massWindow0);
								const double massDiff0 = std::fabs(pi0Candidate0.M() - mass0);
								if (massDiff0 > massWindow0) {
									continue;
								}
								for (size_t m = i + 1; m < nmbPhotons - 1; ++m) {
									for (size_t n = m + 1; n < nmbPhotons; ++n) {
										// exclude photon pairs that share photon(s) with pair 0
										if (m == j or n == j) {
											continue;
										}
										// photon pair 1
										const TLorentzVector pi0Candidate1 = _PhotonLVs[m] + _PhotonLVs[n];
										const std::vector<TLorentzVector> gammasForpi0Candidate1 = {_PhotonLVs[m], _PhotonLVs[n]};
										double mass1 = 0.0;
										double massWindow1 = 0.0;
										getECALMassParamters(_ECALClusterIndices[m], _ECALClusterIndices[n], _ECAL1Mass, _ECAL1MassWindow, _ECAL2Mass, _ECAL2MassWindow, _ECALMixedMass, _ECALMixedMassWindow, mass1, massWindow1);
										const double massDiff1 = std::fabs(pi0Candidate0.M() - mass1);
										// elliptic cut in mass vs mass plane
										if (massDiff0 * massDiff0 / (massWindow0 * massWindow0) + massDiff1 * massDiff1 / (massWindow1 * massWindow1) > 1) {
											continue;
										}
										if (_ResultNmbGoodPi0Pairs == 0) {
											_ResultPi0PairLVs.push_back(pi0Candidate0);
											_ResultGammaLVsForPi0_0 = gammasForpi0Candidate0;
											_ResultPi0PairLVs.push_back(pi0Candidate1);
											_ResultGammaLVsForPi0_1 = gammasForpi0Candidate1;
											_ResultSelectedClusterIndices = {(int)i, (int)j, (int)m, (int)n};
											// find combination type of both photon pairs: 1 if both photons in ECAL1, 2 if both photons in ECAL2, 3 if mixed case
											for (size_t indexEntry = 0; indexEntry < 4; indexEntry = indexEntry + 2) {
												int index0 = _ECALClusterIndices[_ResultSelectedClusterIndices[indexEntry]];
												int index1 = _ECALClusterIndices[_ResultSelectedClusterIndices[indexEntry+1]];
												switch (index0 + index1) {
													case 2: { // both photons in ECAL1
														_ResultPi0CombTypes.push_back(1);
														break;
													}
													case 4: { // both photons in ECAL2
														_ResultPi0CombTypes.push_back(2);
														break;
													}
													case 3: { // one photon in ECAL1, one in ECAL2
														_ResultPi0CombTypes.push_back(3);
														break;
													}
													default:
														abort();
												}
											}
										}
										_ResultNmbGoodPi0Pairs++;
									}
								}
							}
						}

						return true;
					}

				private:

					const std::vector<TLorentzVector>& _PhotonLVs;
					const std::vector<int>&            _ECALClusterIndices;
					const double                       _ECALMixedMass;        // constant parameter, needs to be copied
					const double                       _ECALMixedMassWindow;  // constant parameter, needs to be copied
					const double                       _ECAL1Mass;            // constant parameter, needs to be copied
					const double                       _ECAL1MassWindow;      // constant parameter, needs to be copied
					const double                       _ECAL2Mass;            // constant parameter, needs to be copied
					const double                       _ECAL2MassWindow;      // constant parameter, needs to be copied
					std::vector<TLorentzVector>&       _ResultPi0PairLVs;
					std::vector<int>&                  _ResultPi0CombTypes;
					int&                               _ResultNmbGoodPi0Pairs;
					std::vector<int>&                  _ResultSelectedClusterIndices;
					std::vector<TLorentzVector>&       _ResultGammaLVsForPi0_0;
					std::vector<TLorentzVector>&       _ResultGammaLVsForPi0_1;

				};


				// perform kinematic fit of two photon pairs to given mass
				class GetKinematicFittingMass : public Function
				{

				public:

					GetKinematicFittingMass(const TVector3&              VertexPosition,            // position of primary vertex
					                        const std::vector<TVector3>& ClusterPositions,          // positions of ECAL clusters
					                        const std::vector<TVector3>& ClusterPositionVariances,  // variances in position of ECAL clusters
					                        const std::vector<double>&   ClusterEnergies,           // energies of ECAL clusters
					                        const std::vector<double>&   ClusterEnergyVariances,    // variance in energy of ECAL clusters
					                        const std::vector<int>&      ClusterECALIndices,        // indices in the arrays above of ECAL clusters in the two photon pairs
					                        const double&                Mass,                      // mass that each of the two photon pairs is constrained to
					                        const double&                PrecisionGoal,             // defines convergence criterion for fit
					                        const int&                   WhichEnergyVariance,       // defines how variance of ECAL energies is calculated; see src/neutral_fit.h
					                        std::vector<TLorentzVector>& ResultLorentzVectors,      // Lorentz vectors of photon pairs after fit
					                        std::vector<double>&         ResultChi2s,               // chi^2 values of fits
					                        std::vector<double>&         ResultPValues,             // P-values of fits
					                        std::vector<int>&            ResultNmbIterations,       // number of iterations required to reach PrecisionGoal
					                        int&                         ResultSuccess,             // indicates whether fit was successful
					                        std::vector<double>&         ResultPullsX0,             // pulls for x direction of first photon in pairs
					                        std::vector<double>&         ResultPullsY0,             // pulls for y direction of first photon in pairs
					                        std::vector<double>&         ResultPullsE0,             // pulls for energie of first photon in pairs
					                        std::vector<double>&         ResultPullsX1,             // pulls for x direction of second photon in pairs
					                        std::vector<double>&         ResultPullsY1,             // pulls for y direction of second photon in pairs
					                        std::vector<double>&         ResultPullsE1)             // pulls for energy of second photon in pairs
						: _VertexPosition          (VertexPosition),
						  _ClusterPositions        (ClusterPositions),
						  _ClusterPositionVariances(ClusterPositionVariances),
						  _ClusterEnergies         (ClusterEnergies),
						  _ClusterEnergyVariances  (ClusterEnergyVariances),
						  _ClusterECALIndices      (ClusterECALIndices),
						  _Mass                    (Mass),
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

						// Fit clusters corresponding to first and second pair of indices in _ClusterECALIndices to given _Mass
						std::vector<bool> successes = {false, false};
						for (size_t i = 0; i < 2; ++i) {
							const int clusterIndexA = _ClusterECALIndices[2 * i];
							const int clusterIndexB = _ClusterECALIndices[2 * i + 1];
							antok::NeutralFit neutralFit(
								_VertexPosition,
								_ClusterPositions        [clusterIndexA],
								_ClusterPositions        [clusterIndexB],
								_ClusterPositionVariances[clusterIndexA],
								_ClusterPositionVariances[clusterIndexB],
								_ClusterEnergies         [clusterIndexA],
								_ClusterEnergies         [clusterIndexB],
								_ClusterEnergyVariances  [clusterIndexA],
								_ClusterEnergyVariances  [clusterIndexB],
								_Mass,
								0, // MassLowerLimit for neutralFit::massIsInWindow(), which is not used
								0, // MassUpperLimit for neutralFit::massIsInWindow(), which is not used
								_PrecisionGoal,
								_WhichEnergyVariance);
							successes[i] = neutralFit.doFit();
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
					         int&                  ResultNmbOmegas,         // 1 if there is exactly one omega candidate, 0 otherwise
					         TLorentzVector&       ResultNotUsedPi0LV,      // Lorentz vector of the pi^0 that is not part of the omega
					         TLorentzVector&       ResultNotUsedPiMinusLV,  // Lorentz vector of the pi^- that is not part of the omega
					         TLorentzVector&       ResultPiMinusInOmegaLV,  // Lorentz vector of the pi^- that is part of the omega
					         TLorentzVector&       ResultPi0InOmegaLV,      // Lorentz vector of the pi^0 that is part of the omega
					         TLorentzVector&       ResultPiPlusInOmegaLV)   // Lorentz vector of the pi^+ that is part of the omega
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
						  _ResultNmbOmegas       (ResultNmbOmegas),
						  _ResultNotUsedPi0LV    (ResultNotUsedPi0LV),
						  _ResultNotUsedPiMinusLV(ResultNotUsedPiMinusLV),
						  _ResultPiMinusInOmegaLV(ResultPiMinusInOmegaLV),
						  _ResultPi0InOmegaLV    (ResultPi0InOmegaLV),
						  _ResultPiPlusInOmegaLV (ResultPiPlusInOmegaLV)
					{ }

					virtual ~GetOmega() { }

					bool
					operator() ()
					{
						// searches for omega(782) -> pi^- pi^0 pi^+ candidates in pi^- pi^- pi^0 pi^0 pi^+ final state
						// if several omega candidates are found, a random omega is returned
						const std::vector<const TLorentzVector*> pi0s       = {&_Pi0LV_O, &_Pi0LV_1};
						const std::vector<const TLorentzVector*> chargedLVs = {&_ChargedPartLV_0, &_ChargedPartLV_1, &_ChargedPartLV_2};
						const std::vector<const int*>            charges    = {&_Charge_0,        &_Charge_1,        &_Charge_2};
						_ResultNmbOmegas = 0;
						std::vector<TLorentzVector> _OmegaLVs = {};
						std::vector<TLorentzVector> _Pi0InOmegaLVs = {};
						std::vector<TLorentzVector> _NotUsedPi0LVs = {};
						std::vector<TLorentzVector> _PiMinusInOmegaLVs = {};
						std::vector<TLorentzVector> _NotUsedPiMinusLVs = {};
						std::vector<TLorentzVector> _PiPlusInOmegaLVs = {};
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
											_OmegaLVs.push_back(candidate);
											_ResultNmbOmegas++;  // Count omega candidates
											// find pi^0 that is not part of the omega candidate
											for (size_t l = 0; l < pi0s.size(); l++) {
												if (i != l) {
													_NotUsedPi0LVs.push_back(*pi0s[l]);
													_Pi0InOmegaLVs.push_back(*pi0s[i]);
												}
											}
											// find pi^- that is not part of the omega candidate
											for (size_t m = 0; m < chargedLVs.size(); ++m) {
												if ((m != j) and (m != k)) {
													_NotUsedPiMinusLVs.push_back(*chargedLVs[m]);
													if (*charges[j] == +1) {
														_PiMinusInOmegaLVs.push_back(*chargedLVs[k]);
														_PiPlusInOmegaLVs.push_back(*chargedLVs[j]);
													} else {
														_PiMinusInOmegaLVs.push_back(*chargedLVs[j]);
														_PiPlusInOmegaLVs.push_back(*chargedLVs[k]);
													}
												}
											}
										}
									}
								}
							}
						}
						if (_ResultNmbOmegas == 0) return true;
						if (_ResultNmbOmegas == 1) {
							_ResultOmegaLV          = _OmegaLVs[0];
							_ResultNotUsedPi0LV     = _NotUsedPi0LVs[0];
							_ResultNotUsedPiMinusLV = _NotUsedPiMinusLVs[0];
							_ResultPiMinusInOmegaLV = _PiMinusInOmegaLVs[0];
							_ResultPi0InOmegaLV     = _Pi0InOmegaLVs[0];
							_ResultPiPlusInOmegaLV  = _PiPlusInOmegaLVs[0]; 
						} else {
							// randomly choose one of the events with a valid omega candidate
							TRandom* randomSelector = new TRandom();
							const size_t selectedEvent = randomSelector->Integer(_ResultNmbOmegas);
							delete randomSelector;
							_ResultOmegaLV          = _OmegaLVs[selectedEvent];
							_ResultNotUsedPi0LV     = _NotUsedPi0LVs[selectedEvent];
							_ResultNotUsedPiMinusLV = _NotUsedPiMinusLVs[selectedEvent];
							_ResultPiMinusInOmegaLV = _PiMinusInOmegaLVs[selectedEvent];
							_ResultPi0InOmegaLV     = _Pi0InOmegaLVs[selectedEvent];
							_ResultPiPlusInOmegaLV  = _PiPlusInOmegaLVs[selectedEvent]; 
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
					int&                  _ResultNmbOmegas;
					TLorentzVector&       _ResultNotUsedPi0LV;
					TLorentzVector&       _ResultNotUsedPiMinusLV;
					TLorentzVector&       _ResultPiMinusInOmegaLV;
					TLorentzVector&       _ResultPi0InOmegaLV;
					TLorentzVector&       _ResultPiPlusInOmegaLV;

				};


				class GetPionLVs : public Function
				{

				public:

					GetPionLVs(const TLorentzVector&        Pi0LV_0,          // Lorentz vector of 1st pi^0
					           const TLorentzVector&        Pi0LV_1,          // Lorentz vector of 2nd pi^0
					           const TLorentzVector&        ChargedPartLV_0,  // Lorentz vector of 1st charged particle
					           const TLorentzVector&        ChargedPartLV_1,  // Lorentz vector of 2nd charged particle
					           const TLorentzVector&        ChargedPartLV_2,  // Lorentz vector of 3rd charged particle
					           const int&                   Charge_0,         // charge of 1st charged particle
					           const int&                   Charge_1,         // charge of 2nd charged particle
					           const int&                   Charge_2,         // charge of 3rd charged particle
					           const int&                   SelectedCharge,
					           std::vector<TLorentzVector>& Result)           // result LVs
						: _Pi0LV_0        (Pi0LV_0),
						  _Pi0LV_1        (Pi0LV_1),
						  _ChargedPartLV_0(ChargedPartLV_0),
						  _ChargedPartLV_1(ChargedPartLV_1),
						  _ChargedPartLV_2(ChargedPartLV_2),
						  _Charge_0       (Charge_0),
						  _Charge_1       (Charge_1),
						  _Charge_2       (Charge_2),
						  _SelectedCharge (SelectedCharge),
						  _Result         (Result)
					{ }

					virtual ~GetPionLVs() { }

					bool
					operator() ()
					{
						_Result.clear();
						_Result.reserve(2);
						const std::vector<const TLorentzVector*> ChargedPartLVs = {&_ChargedPartLV_0, &_ChargedPartLV_1, &_ChargedPartLV_2};
						const std::vector<int>                   Charges        = { _Charge_0,         _Charge_1,         _Charge_2};

						if (_SelectedCharge == 0) {
							_Result = {_Pi0LV_0, _Pi0LV_1};
						} else {
							for (size_t i = 0; i < ChargedPartLVs.size(); ++i) {
								if (Charges[i] == _SelectedCharge) {
									_Result.push_back(*(ChargedPartLVs[i]));
								}
							}
						}

						return true;
					}

				private:

					const TLorentzVector&        _Pi0LV_0;
					const TLorentzVector&        _Pi0LV_1;
					const TLorentzVector&        _ChargedPartLV_0;
					const TLorentzVector&        _ChargedPartLV_1;
					const TLorentzVector&        _ChargedPartLV_2;
					const int&                   _Charge_0;
					const int&                   _Charge_1;
					const int&                   _Charge_2;
					const int                    _SelectedCharge;
					std::vector<TLorentzVector>& _Result;

				};


				class GetTwoPionCombinationLV : public Function
				{

				public:

					GetTwoPionCombinationLV(const TLorentzVector&        Pi0LV_0,          // Lorentz vector of 1st pi^0
					                        const TLorentzVector&        Pi0LV_1,          // Lorentz vector of 2nd pi^0
					                        const TLorentzVector&        ChargedPartLV_0,  // Lorentz vector of 1st charged particle
					                        const TLorentzVector&        ChargedPartLV_1,  // Lorentz vector of 2nd charged particle
					                        const TLorentzVector&        ChargedPartLV_2,  // Lorentz vector of 3rd charged particle
					                        const int&                   Charge_0,         // charge of 1st charged particle
					                        const int&                   Charge_1,         // charge of 2nd charged particle
					                        const int&                   Charge_2,         // charge of 3rd charged particle
					                        const int&                   CombinationMode,
					                        std::vector<TLorentzVector>& Result)           // result LVs
						: _Pi0LV_0        (Pi0LV_0),
						  _Pi0LV_1        (Pi0LV_1),
						  _ChargedPartLV_0(ChargedPartLV_0),
						  _ChargedPartLV_1(ChargedPartLV_1),
						  _ChargedPartLV_2(ChargedPartLV_2),
						  _Charge_0       (Charge_0),
						  _Charge_1       (Charge_1),
						  _Charge_2       (Charge_2),
						  _CombinationMode(CombinationMode),
						  _Result         (Result)
					{ }

					virtual ~GetTwoPionCombinationLV() { }

					bool
					operator() ()
					{
						enum combMode {Pi0PiMinusCombinations = 0, Pi0PiPlusCombinations = 1, PiMinusPiPlusCombinations = 2};

						_Result.clear();
						_Result.reserve(4);
						const std::vector<const TLorentzVector*> Pi0LVs         = {&_Pi0LV_0, &_Pi0LV_1};
						const std::vector<const TLorentzVector*> ChargedPartLVs = {&_ChargedPartLV_0, &_ChargedPartLV_1, &_ChargedPartLV_2};
						const std::vector<int>                   Charges        = { _Charge_0,         _Charge_1,         _Charge_2};

						std::vector<const TLorentzVector*> PiMinusLVs;
						std::vector<const TLorentzVector*> PiPlusLVs;
						for (size_t i = 0; i < ChargedPartLVs.size(); ++i) {
							if (Charges[i] == -1) {
								PiMinusLVs.push_back(ChargedPartLVs[i]);
							} else if (Charges[i] == 1) {
								PiPlusLVs.push_back(ChargedPartLVs[i]);
							}
						}

						//get requested Pi LVs
						std::vector<const TLorentzVector*> LVs_0, LVs_1;
						switch ((combMode)_CombinationMode) {
							case Pi0PiMinusCombinations: {
								LVs_0 = Pi0LVs;
								LVs_1 = PiMinusLVs;
								break;
							}
							case Pi0PiPlusCombinations: {
								LVs_0 = Pi0LVs;
								LVs_1 = PiPlusLVs;
								break;
							}
							case PiMinusPiPlusCombinations: {
								LVs_0 = PiMinusLVs;
								LVs_1 = PiPlusLVs;
								break;
							}
							default:
								abort();
						}

						// get all Pi combinations
						for (const TLorentzVector* LV_0 : LVs_0) {
							for (const TLorentzVector* LV_1 : LVs_1) {
								_Result.push_back(*LV_0 + *LV_1);
							}
						}

						return true;
					}

				private:

					const TLorentzVector&        _Pi0LV_0;
					const TLorentzVector&        _Pi0LV_1;
					const TLorentzVector&        _ChargedPartLV_0;
					const TLorentzVector&        _ChargedPartLV_1;
					const TLorentzVector&        _ChargedPartLV_2;
					const int&                   _Charge_0;
					const int&                   _Charge_1;
					const int&                   _Charge_2;
					const int                    _CombinationMode;
					std::vector<TLorentzVector>& _Result;

				};


				class GetThreePionCombinationLV : public Function
				{

				public:

					GetThreePionCombinationLV(const TLorentzVector&        Pi0LV_0,           // Lorentz vector of 1st pi^0
					                          const TLorentzVector&        Pi0LV_1,           // Lorentz vector of 2nd pi^0
					                          const TLorentzVector&        ChargedPartLV_0,   // Lorentz vector of 1st charged particle
					                          const TLorentzVector&        ChargedPartLV_1,   // Lorentz vector of 2nd charged particle
					                          const TLorentzVector&        ChargedPartLV_2,   // Lorentz vector of 3rd charged particle
					                          const int&                   Charge_0,          // charge of 1st charged particle
					                          const int&                   Charge_1,          // charge of 2nd charged particle
					                          const int&                   Charge_2,          // charge of 3rd charged particle
					                          std::vector<TLorentzVector>& Result3PiLVs,      // result 3Pi LVs
											  std::vector<TLorentzVector>& ResultPi0In3PiLVs) // result Pi0 in 3Pi comb LVs
						: _Pi0LV_0          (Pi0LV_0),
						  _Pi0LV_1          (Pi0LV_1),
						  _ChargedPartLV_0  (ChargedPartLV_0),
						  _ChargedPartLV_1  (ChargedPartLV_1),
						  _ChargedPartLV_2  (ChargedPartLV_2),
						  _Charge_0         (Charge_0),
						  _Charge_1         (Charge_1),
						  _Charge_2         (Charge_2),
						  _Result3PiLVs     (Result3PiLVs),
						  _ResultPi0In3PiLVs(ResultPi0In3PiLVs)

					{ }

					virtual ~GetThreePionCombinationLV() { }

					bool
					operator() ()
					{
						_Result3PiLVs.clear();
						_Result3PiLVs.reserve(4);
						_ResultPi0In3PiLVs.clear();
						_ResultPi0In3PiLVs.reserve(4);
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
										_Result3PiLVs.push_back(*Pi0LV + *chargedFirst + *chargedSecond);
										_ResultPi0In3PiLVs.push_back(*Pi0LV);
									}
								}
							}
						}

						return true;
					}

				private:

					const TLorentzVector&        _Pi0LV_0;
					const TLorentzVector&        _Pi0LV_1;
					const TLorentzVector&        _ChargedPartLV_0;
					const TLorentzVector&        _ChargedPartLV_1;
					const TLorentzVector&        _ChargedPartLV_2;
					const int&                   _Charge_0;
					const int&                   _Charge_1;
					const int&                   _Charge_2;
					std::vector<TLorentzVector>& _Result3PiLVs;
					std::vector<TLorentzVector>& _ResultPi0In3PiLVs;

				};


				class GetFourPionCombinationLV : public Function
				{

				public:

					GetFourPionCombinationLV(const TLorentzVector&        Pi0LV_0,          // Lorentz vector of 1st pi^0
					                         const TLorentzVector&        Pi0LV_1,          // Lorentz vector of 2nd pi^0
					                         const TLorentzVector&        ChargedPartLV_0,  // Lorentz vector of 1st charged particle
					                         const TLorentzVector&        ChargedPartLV_1,  // Lorentz vector of 2nd charged particle
					                         const TLorentzVector&        ChargedPartLV_2,  // Lorentz vector of 3rd charged particle
					                         const int&                   Charge_0,         // charge of 1st charged particle
					                         const int&                   Charge_1,         // charge of 2nd charged particle
					                         const int&                   Charge_2,         // charge of 3rd charged particle
					                         const int&                   CombinationMode,
					                         std::vector<TLorentzVector>& Result)           // result LVs
						: _Pi0LV_0        (Pi0LV_0),
						  _Pi0LV_1        (Pi0LV_1),
						  _ChargedPartLV_0(ChargedPartLV_0),
						  _ChargedPartLV_1(ChargedPartLV_1),
						  _ChargedPartLV_2(ChargedPartLV_2),
						  _Charge_0       (Charge_0),
						  _Charge_1       (Charge_1),
						  _Charge_2       (Charge_2),
						  _CombinationMode(CombinationMode),
						  _Result         (Result)
					{ }

					virtual ~GetFourPionCombinationLV() { }

					bool
					operator() ()
					{
						enum combMode {Pi00MinusPlusCombinations = -1, Pi0MinusMinusPlusCombinations = 0, Pi00MinusMinusCombinations = 1};

						_Result.clear();
						_Result.reserve(2);
						const std::vector<const TLorentzVector*> Pi0LVs         = {&_Pi0LV_0, &_Pi0LV_1};
						const std::vector<const TLorentzVector*> ChargedPartLVs = {&_ChargedPartLV_0, &_ChargedPartLV_1, &_ChargedPartLV_2};
						const std::vector<int>                   Charges        = { _Charge_0,         _Charge_1,         _Charge_2};

						std::vector<const TLorentzVector*> PiMinusLVs;
						std::vector<const TLorentzVector*> PiPlusLVs;
						for (size_t i = 0; i < ChargedPartLVs.size(); ++i) {
							if (Charges[i] == -1) {
								PiMinusLVs.push_back(ChargedPartLVs[i]);
							} else if (Charges[i] == 1) {
								PiPlusLVs.push_back(ChargedPartLVs[i]);
							}
						}

						//get requested Pi LVs
						switch ((combMode)_CombinationMode) {
							case Pi00MinusPlusCombinations: {
								for (auto& selectedLV: PiMinusLVs) {
									_Result.push_back(*selectedLV + *(PiPlusLVs[0]) + _Pi0LV_0 + _Pi0LV_1);
								}
								break;
							}
							case Pi0MinusMinusPlusCombinations: {
								for (auto& selectedLV: Pi0LVs) {
									_Result.push_back(*selectedLV + *(PiPlusLVs[0]) + *(PiMinusLVs[0]) + *(PiMinusLVs[1]));
								}
								break;
							}
							case Pi00MinusMinusCombinations: {
								_Result.push_back(*(PiMinusLVs[0]) + *(PiMinusLVs[1]) + _Pi0LV_0 + _Pi0LV_1);
								break;
							}
							default:
								abort();
						}

						return true;
					}

				private:

					const TLorentzVector&        _Pi0LV_0;
					const TLorentzVector&        _Pi0LV_1;
					const TLorentzVector&        _ChargedPartLV_0;
					const TLorentzVector&        _ChargedPartLV_1;
					const TLorentzVector&        _ChargedPartLV_2;
					const int&                   _Charge_0;
					const int&                   _Charge_1;
					const int&                   _Charge_2;
					const int                    _CombinationMode;
					std::vector<TLorentzVector>& _Result;

				};


			}  // functions namespace

		}  // cdreis namespace

	}  // user namespace

}  // antok namespace

#endif  // ANTOK_USER_CDREIS_FUNCTIONS_HPP
