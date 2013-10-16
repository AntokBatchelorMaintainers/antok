#ifndef ANTOK_USER_KBICKER_FUNCTIONS_HPP
#define ANTOK_USER_KBICKER_FUNCTIONS_HPP

#include<iostream>
#include<vector>

#include<TLorentzVector.h>

#include<basic_calcs.h>

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
								antok::getRPDDeltaPhiResProjection((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr));
								return true;
							case 1:
								if(_protonPhiAddr == 0 and _xPhiAddr == 0) {
									antok::getRPDDeltaPhiResRotation((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr));
								} else {
									antok::getRPDDeltaPhiResRotation((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr), (*_protonPhiAddr), (*_xPhiAddr));
								}
								return true;
							case 2:
								if(not _vertex) {
									std::cerr<<"Found null pointer for vertex in GetRpdPhi while using prediction method."<<std::endl;
									return false;
								}
								if(_protonPhiAddr == 0 and _xPhiAddr == 0) {
									antok::getRPDDeltaPhiResPrediction((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_vertex), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr));
								} else {
									antok::getRPDDeltaPhiResPrediction((*_beamLorentzVecAddr), (*_rpdProtonLorentzVecAddr), (*_xLorentzVecAddr), (*_vertex), (*_rpdDeltaPhiAddr), (*_rpdDeltaPhiResAddr), (*_protonPhiAddr), (*_xPhiAddr));
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
						antok::getRPDExpectedHitsParameters(*_pBeamAddr,
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

			}

		}

	}

}

#endif
