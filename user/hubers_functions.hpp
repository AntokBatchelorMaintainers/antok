#ifndef ANTOK_USER_HUBERS_FUNCTIONS_HPP
#define ANTOK_USER_HUBERS_FUNCTIONS_HPP

#include<iostream>
#include<vector>
#include<NNpoly.h>

#include<TLorentzVector.h>

namespace antok {

	namespace user {

		namespace hubers {

			namespace functions {

				class Sqrt : public Function
				{
					public:
						Sqrt(double* inAddr, double* outAddr)
							: _inAddr(inAddr),
							_outAddr(outAddr) { }

						virtual ~Sqrt() { }

						bool operator() () {
							if(*_inAddr>0)
								(*_outAddr) = std::sqrt(*_inAddr);
							else
								(*_outAddr) = -std::sqrt(-*_inAddr);
							return true;
						}

					private:
						double* _inAddr;
						double* _outAddr;
				};

				class Frac : public Function
				{
					public:
						Frac(double* inAddr1, double* inAddr2, double* outAddr)
							: _inAddr1(inAddr1),
							_inAddr2(inAddr2),
							_outAddr(outAddr) { }

						virtual ~Frac() { }

						bool operator() () {
							if(*_inAddr2!=0)
								(*_outAddr) = (*_inAddr1) / (*_inAddr2);
							else
								(*_outAddr) = (1<<31);
							return true;
						}

					private:
						double* _inAddr1;
						double* _inAddr2;
						double* _outAddr;
				};



				class GetPt : public Function
				{
					public:
						GetPt(TLorentzVector* pLorentzVec, TLorentzVector* beamLorentzVec, double *pTAddr)
							: _pLorentzVec(pLorentzVec),
							_beamLorentzVec(beamLorentzVec),
							_pTAddr(pTAddr){}

						virtual ~GetPt() { }

						bool operator() () {
							TVector3 beamVec = _beamLorentzVec->Vect();
							TVector3 pVec = _pLorentzVec->Vect();
							(*_pTAddr) = beamVec.Cross( pVec ) .Cross( beamVec ) .Unit() .Dot( pVec );
							return true;
						}

					private:
						TLorentzVector* _pLorentzVec;
						TLorentzVector* _beamLorentzVec;
						double* _pTAddr;
				};

				template<typename T>
				class EnforceEConservation : public Function
				{
					public:
						EnforceEConservation(T* beamAddr, T* pionAddr, T* gammaAddr,  T* outAddr) {
							_beamAddr = beamAddr;
							_pionAddr = pionAddr;
							_gammaAddr = gammaAddr;
							_outAddr = outAddr;
							_mode=0;

						};

						virtual ~EnforceEConservation() { }

						bool operator() () {
							if(_mode==0){
								const double E = _beamAddr->E()  - _pionAddr->E();
								TVector3 g3( _gammaAddr->Vect() );
								g3.SetMag(E);
								_outAddr->SetVect(g3);
								_outAddr->SetE(E);
							}
							else if(_mode==1){
								const double E = _beamAddr->E()  - _gammaAddr->E();
								TVector3 pi3( _pionAddr->Vect() );
								pi3.SetMag(sqrt(E*E- 0.13957018*0.13957018));
								_outAddr->SetVect(pi3);
								_outAddr->SetE(E);
							}
							return true;
						};


					private:

						T* _beamAddr;
						T* _pionAddr;
						T* _gammaAddr;
						T* _outAddr;
						int _mode;
				};

				//***********************************
				//Calculates the beam energy using a
				//polynomial represention from a
				//neuronal network
				//calibrated for primakoff 2009
				//***********************************
				class GetNeuronalBeam : public Function
				{
					public:
						GetNeuronalBeam(double* xAddr, double* yAddr, double* dxAddr,
						                double* dyAddr, double* eAddr, TLorentzVector* LVAddr)
						               :_xAddr(xAddr), _yAddr(yAddr), _dxAddr(dxAddr), _dyAddr(dyAddr),
						                _eAddr(eAddr), _LVAddr(LVAddr){}

						virtual ~GetNeuronalBeam() { }

						bool operator() () {
							double xarr[4]={*_xAddr,*_yAddr,*_dxAddr,*_dyAddr};
							if ( std::fabs(*_xAddr) > 1.8 || std::fabs(*_yAddr) > 1.8 || std::fabs(*_dxAddr) > 5e-4 || std::fabs(*_dyAddr+3e-4) > 5e-4 )
								*_eAddr = 190;
							else
								*_eAddr = NNpoly::Ebeam(xarr,NNpoly::getParams2009());
							TVector3 v3(*_dxAddr, *_dyAddr, std::sqrt( 1 - sqr(*_dxAddr) - sqr(*_dyAddr) ));
							v3.SetMag(*_eAddr);
							_LVAddr->SetXYZT(v3.X(),v3.Y(), v3.Z(), std::sqrt( sqr(*_eAddr) + sqr(antok::Constants::chargedPionMass())) );
							return true;
						}

					private:
						double* _xAddr;
						double* _yAddr;
						double* _dxAddr;
						double* _dyAddr;
						double* _eAddr;
						TLorentzVector* _LVAddr;
				};



			}

		}

	}

}
#endif
