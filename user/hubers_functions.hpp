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

				class GetNeuronalBeamEnergy : public Function
				{
					public:
						GetNeuronalBeamEnergy(double* xAddr,
								double* yAddr,
								double* dxAddr,
								double* dyAddr,
								double* eAddr)
							: _xAddr(xAddr),
							_yAddr(yAddr),
							_dxAddr(dxAddr),
							_dyAddr(dyAddr),
							_eAddr(eAddr) { }

						virtual ~GetNeuronalBeamEnergy() { }

						bool operator() () {
							double xarr[4]={*_xAddr,*_yAddr,*_dxAddr,*_dyAddr};
							*_eAddr=					antok::user::hubers::NNpoly::Ebeam(xarr,NNpoly::params);
							return true;
						}

					private:
						double* _xAddr;
						double* _yAddr;
						double* _dxAddr;
						double* _dyAddr;
						double* _eAddr;
				};



			}

		}

	}

}


#endif
