#ifndef ANTOK_FUNCTIONS_H
#define ANTOK_FUNCTIONS_H

#include<iostream>
#include<vector>

#include<TLorentzVector.h>

namespace antok {

	class Function
	{

	  public:

		virtual bool operator() () = 0;

	};

	template<typename T>
	class Sum : public Function
	{

	  public:

	  	Sum(std::vector<T*> inputAddrs, T* outAddr) {
			_inputAddrs = inputAddrs;
			_outAddr = outAddr;
			if(_inputAddrs.size() < 1) {
				std::cerr<<"Got empty address vector as input for a sum."<<std::endl;
				throw 1;
			}
		};
		
		bool operator() () {
			(*_outAddr) = *(_inputAddrs[0]);
			for(unsigned int i = 1; i < _inputAddrs.size(); ++i) {
				*_outAddr = *_outAddr + *(_inputAddrs.at(i));
			}
			return true;
		};


	  private:

		std::vector<T*> _inputAddrs;
		T* _outAddr;

	};

	class Mass: public Function
	{

	  public:
		
		Mass(TLorentzVector* inputAddr, double* outAddr) {
			_inputAddr = inputAddr;
			_outAddr = outAddr;
		}

		bool operator() () {
			(*_outAddr) = _inputAddr->M();
			return true;
		}

	  private:

		TLorentzVector* _inputAddr;
		double* _outAddr;

	};

	class GetLorentzVec: public Function
	{

	  public:

		GetLorentzVec(double* xAddr, double* yAddr, double* zAddr, double* mAddr, TLorentzVector* outAddr) {
			_xAddr = xAddr;
			_yAddr = yAddr;
			_zAddr = zAddr;
			_mAddr = mAddr;
			_outAddr = outAddr;
		}

		bool operator() () {
			_outAddr->SetXYZM(*_xAddr, *_yAddr, *_zAddr, *_mAddr);
//			_outAddr->Print();
			return true;
		}

	  private:

		double* _xAddr;
		double* _yAddr;
		double* _zAddr;
		double* _mAddr;
		TLorentzVector* _outAddr;

	};

}

#endif

